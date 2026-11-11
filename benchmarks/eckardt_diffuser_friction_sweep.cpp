#include <turbocompgpu/complete_stage_solver.hpp>
#include <turbocompgpu/geometry_coupled_choking.hpp>
#include <turbocompgpu/choking_mass_flux.hpp>
#include <turbocompgpu/throat_geometry.hpp>
#include <turbocompgpu/vaneless_diffuser_variable_width.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/caloric.hpp>
#include <thermomodel/peng_robinson.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numbers>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using namespace turbocompgpu;

namespace {

struct MapPoint {
    double rpm;
    double mass_flow;
    double measured_pressure_ratio;
};

struct EfficiencyPoint {
    double rpm;
    double mass_flow;
    double measured_efficiency;
};

std::vector<EfficiencyPoint> read_efficiency_map(
    const std::filesystem::path& path)
{
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open efficiency data: " + path.string());
    std::vector<EfficiencyPoint> points;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line.front() == '#')
            continue;
        std::istringstream row(line);
        EfficiencyPoint point{};
        char comma1 = '\0';
        char comma2 = '\0';
        if (!(row >> point.rpm >> comma1 >> point.mass_flow >> comma2
                  >> point.measured_efficiency) ||
            comma1 != ',' || comma2 != ',')
            throw std::runtime_error("Invalid efficiency row: " + line);
        points.push_back(point);
    }
    if (points.empty())
        throw std::runtime_error("Efficiency map contains no points");
    return points;
}

std::vector<MapPoint> read_map(const std::filesystem::path& path)
{
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open map data: " + path.string());

    std::vector<MapPoint> points;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line.front() == '#')
            continue;
        std::istringstream row(line);
        MapPoint point{};
        char comma1 = '\0';
        char comma2 = '\0';
        if (!(row >> point.rpm >> comma1 >> point.mass_flow >> comma2
                  >> point.measured_pressure_ratio) ||
            comma1 != ',' || comma2 != ',') {
            throw std::runtime_error("Invalid map row: " + line);
        }
        points.push_back(point);
    }
    return points;
}

ImpellerSystemGeometry geometry()
{
    // Eckardt backswept Impeller A. Do not compare this geometry directly
    // against the radial-discharge Impeller O traverse dataset.
    constexpr double exit_diameter_m = 0.400;
    constexpr double exit_width_m = 0.026;
    constexpr double blade_count = 20.0;
    constexpr double exit_blade_thickness_m = 0.00108;
    constexpr double exit_blockage = blade_count * exit_blade_thickness_m /
        (2.0 * std::numbers::pi * (exit_diameter_m / 2.0));

    return {
        {0.120, 0.280, 0.0},
        {exit_diameter_m, exit_width_m, exit_blockage},
        {exit_diameter_m, 60.0,
         2.0 * std::sqrt(0.5 * (0.140 * 0.140 + 0.060 * 0.060))}
    };
}

NonlinearResidualScales scales()
{
    return {50000.0, 100.0, 300000.0};
}

OperatingPointRequest make_request(
    const MapPoint& point,
    const thermomodel::Model& air,
    const thermomodel::DifferentiableModel& derivatives,
    const ImpellerSystemGeometry& g)
{
    constexpr double T01 = 288.15;
    constexpr double p01 = 101000.0;
    constexpr double beta_author_deg = 30.0;
    constexpr double blades = 20.0;

    const auto inlet_total = evaluate(air, T01, p01);
    const double area1 = effective_flow_area(
        inlet_annulus_area(g.inlet), g.inlet.blockage);
    const double maximum_velocity = 0.999 *
        evaluate_isentropic_sound_speed(derivatives, T01, p01).sound_speed;
    const auto mass_residual = [&](double meridional_velocity) {
        const auto state = state_from_hs(
            air, inlet_total.h - 0.5 * meridional_velocity * meridional_velocity,
            inlet_total.s, T01, p01);
        return state.rho * area1 * meridional_velocity - point.mass_flow;
    };
    if (!(maximum_velocity > 0.0 && mass_residual(maximum_velocity) >= 0.0))
        throw std::runtime_error("Inlet flow exceeds the subsonic mass-flow bracket");

    double velocity_lo = 0.0;
    double velocity_hi = maximum_velocity;
    for (int i = 0; i < 100; ++i) {
        const double velocity = 0.5 * (velocity_lo + velocity_hi);
        if (mass_residual(velocity) < 0.0)
            velocity_lo = velocity;
        else
            velocity_hi = velocity;
    }
    const double vm1 = 0.5 * (velocity_lo + velocity_hi);
    const auto inlet_static = state_from_hs(
        air, inlet_total.h - 0.5 * vm1 * vm1,
        inlet_total.s, T01, p01);
    const double slip = 1.0 -
        std::sqrt(std::cos(beta_author_deg * std::numbers::pi / 180.0)) /
        std::pow(blades, 0.7);

    return {
        {inlet_static.T, inlet_static.p, point.mass_flow,
         point.rpm, slip, 0.0, 0.0},
        {inlet_static.T * 1.15,
         std::log(1.5 * inlet_static.p / scales().pressure)}
    };
}

struct DiffuserTraverseState {
    double impeller_exit_static_pressure_pa;
    double radius_ratio;
    double radius_m;
    double width_m;
    double bulk_static_pressure_pa;
    double total_pressure_pa;
    double total_temperature_k;
    double flow_angle_from_meridional_deg;
    double radial_velocity_m_s;
    double tangential_velocity_m_s;
};

enum class PublishedLossSet {
    CurrentScreen,
    EckardtA2042,
    CrossGeometry2137
};

struct Prediction {
    double stage_pr;
    double impeller_pr;
    double diffuser_pr;
    double diffuser_loss_percent;
    double darcy_factor;
    double inlet_reynolds;
    double total_to_total_efficiency;
    double impeller_skin_friction_loss_j_kg;
    double impeller_blade_loading_loss_j_kg;
    double impeller_skin_friction_factor;
    double impeller_diffusion_factor;
    double impeller_incidence_loss_j_kg;
    double inlet_relative_angle_deg;
    double optimum_inlet_angle_deg;
    double estimated_throat_area_m2;
    double critical_throat_flow_kg_s;
    double throat_capacity_utilization;
    double throat_capacity_margin_percent;
    bool choke_screen_flag;
    double impeller_specific_work_j_kg;
    double outlet_total_temperature_k;
    double impeller_reynolds;
    double impeller_clearance_loss_j_kg;
    double impeller_mixing_loss_j_kg;
    double impeller_disk_friction_loss_j_kg;
    double impeller_hydraulic_diameter_m;
    double internal_loss_entropy_generation_j_kg_k;
    double impeller_blade_path_length_m;
    double impeller_mean_relative_velocity_m_s;
    double impeller_mean_density_kg_m3;
    double impeller_mean_temperature_k;
    double impeller_skin_fanning_factor;
    double impeller_mixing_chi;
    double impeller_mixing_bstar;
    double impeller_correlation_exit_angle_deg;
    double impeller_inlet_mean_relative_velocity_m_s;
    double impeller_inlet_shroud_relative_velocity_m_s;
    double impeller_exit_relative_velocity_m_s;
    double impeller_throat_relative_velocity_m_s;
    double impeller_throat_density_kg_m3;
    double impeller_throat_enthalpy_j_kg;
    double impeller_throat_area_m2;
    double impeller_throat_blade_speed_m_s;
    double impeller_throat_rothalpy_j_kg;
    double impeller_throat_mass_residual_kg_s;
    double impeller_throat_rothalpy_residual_j_kg;
    double impeller_throat_entropy_residual_j_kg_k;
    double impeller_throat_mach;
    DiffuserTraverseState near_traverse;
    DiffuserTraverseState far_traverse;
    double impeller_choke_loss_j_kg;
    double impeller_entrance_diffusion_loss_j_kg;
    double impeller_recirculation_loss_j_kg;
    double stage_total_enthalpy_rise_j_kg;
    double impeller_euler_work_j_kg;
    double external_loss_enthalpy_addition_j_kg;
    double stage_energy_residual_vs_euler_j_kg;
    double stage_energy_residual_vs_euler_disk_recirc_j_kg;
    double stage_energy_residual_vs_euler_external_addition_j_kg;
};

ThroatGeometry eckardt_meanline_throat(const ImpellerSystemGeometry& g)
{
    constexpr int blade_count = 20;
    constexpr double inlet_blade_thickness_m = 0.00211;
    const double mean_inlet_radius = g.impeller.inlet_diameter / 2.0;
    const double inlet_span =
        (g.inlet.shroud_diameter - g.inlet.hub_diameter) / 2.0;
    // Eckardt A's inlet metal angles are reported from the meridional
    // direction. ThroatGeometry expects its angle from the tangential
    // direction, so convert the mean angle before projecting the pitch onto
    // the blade normal. This remains a local meanline approximation; it does
    // not locate the 3D minimum section or reconstruct spanwise camber.
    constexpr double beta_mean_meridional_deg = 0.5 * (63.0 + 40.0);
    constexpr double beta_mean_from_tangential_deg =
        90.0 - beta_mean_meridional_deg;
    return {blade_count, mean_inlet_radius, inlet_span,
            beta_mean_from_tangential_deg,
            inlet_blade_thickness_m, 0.0};
}

struct ThroatCapacityEstimate {
    double geometric_area_m2;
    double effective_area_m2;
    double critical_mass_flow_kg_s;
};

struct OperatingThroatEstimate {
    double relative_speed_m_s;
    double density_kg_m3;
    double enthalpy_j_kg;
    double effective_area_m2;
    double blade_speed_m_s;
    double inlet_rothalpy_j_kg;
    double mass_residual_kg_s;
    double rothalpy_residual_j_kg;
    double entropy_residual_j_kg_k;
    double mach;
};

using ThroatCapacityCache =
    std::map<std::tuple<double, double, double>, ThroatCapacityEstimate>;

ThroatCapacityEstimate estimate_throat_capacity(
    const thermomodel::Model& air,
    const thermomodel::DifferentiableModel& derivatives,
    const ImpellerSystemGeometry& g,
    const OperatingPointRequest& request,
    ThroatCapacityCache& cache)
{
    const auto& operating = request.operating;
    const auto key = std::make_tuple(
        operating.rotational_speed_rpm,
        operating.inlet_temperature,
        operating.inlet_pressure);
    if (const auto found = cache.find(key); found != cache.end())
        return found->second;

    const auto throat_geometry = eckardt_meanline_throat(g);
    const auto throat_area = evaluate_throat_geometry(throat_geometry);
    const auto choke = evaluate_geometry_coupled_choking(
        air, derivatives, g,
        {operating.inlet_temperature,
         operating.inlet_pressure,
         operating.rotational_speed_rpm,
         operating.inlet_prewhirl,
         throat_geometry});
    const ThroatCapacityEstimate estimate{
        throat_area.geometric_area,
        throat_area.effective_area,
        choke.coupled.critical_mass_flow};
    if (!(std::isfinite(estimate.critical_mass_flow_kg_s) &&
          estimate.critical_mass_flow_kg_s > 0.0))
        throw std::runtime_error("Invalid estimated impeller throat capacity");
    cache.emplace(key, estimate);
    return estimate;
}

// Solve the subcritical, blade-aligned operating throat state. In the rotor
// frame, entropy and rothalpy are conserved and mdot = rho*A_eff*W_n. The
// existing choke solver supplies the sonic upper bracket; the operating
// solution is the subsonic root below it.
OperatingThroatEstimate estimate_operating_throat(
    const thermomodel::Model& air,
    const thermomodel::DifferentiableModel& derivatives,
    const ImpellerSystemGeometry& g,
    const OperatingPointRequest& request)
{
    const auto& op = request.operating;
    const auto inlet = evaluate(air, op.inlet_temperature, op.inlet_pressure);
    const auto geometry = eckardt_meanline_throat(g);
    const auto area = evaluate_throat_geometry(geometry);
    const double inlet_area = effective_flow_area(
        inlet_annulus_area(g.inlet), g.inlet.blockage);
    const double vm1 = op.mass_flow / (inlet.rho * inlet_area);
    const double omega = 2.0 * std::numbers::pi * op.rotational_speed_rpm / 60.0;
    const double u1 = omega * (g.impeller.inlet_diameter / 2.0);
    const double ut = omega * geometry.radius;
    const double wtheta1 = op.inlet_prewhirl - u1;
    const double inlet_rothalpy = inlet.h +
        0.5 * (vm1 * vm1 + wtheta1 * wtheta1) - 0.5 * u1 * u1;
    const double effective_h0 = inlet_rothalpy + 0.5 * ut * ut;
    const auto effective = state_from_hs(
        air, effective_h0, inlet.s, inlet.T, inlet.p);
    const auto sonic = evaluate_choking_mass_flux(
        air, derivatives, effective.T, effective.p);

    const auto state_at_speed = [&](double speed) {
        return state_from_hs(
            air, effective_h0 - 0.5 * speed * speed,
            inlet.s, effective.T, effective.p);
    };
    const auto mass_residual = [&](double speed) {
        return area.effective_area * state_at_speed(speed).rho * speed -
            op.mass_flow;
    };

    double lo = 0.0;
    double hi = sonic.sound_speed;
    if (!(std::isfinite(hi) && hi > 0.0 &&
          mass_residual(hi) >= 0.0))
        throw std::runtime_error(
            "Operating throat flow exceeds its rothalpy-based sonic capacity");
    for (int iteration = 0; iteration < 100; ++iteration) {
        const double mid = 0.5 * (lo + hi);
        if (mass_residual(mid) < 0.0)
            lo = mid;
        else
            hi = mid;
        if (hi - lo <= 1.0e-12 * std::max(1.0, mid))
            break;
    }
    const double speed = 0.5 * (lo + hi);
    const auto throat = state_at_speed(speed);
    const double rothalpy_residual = throat.h + 0.5 * speed * speed -
        0.5 * ut * ut - inlet_rothalpy;
    const double continuity_residual = area.effective_area * throat.rho * speed -
        op.mass_flow;
    const double mach = speed /
        evaluate_isentropic_sound_speed(
            derivatives, throat.T, throat.p).sound_speed;
    return {speed, throat.rho, throat.h, area.effective_area, ut,
            inlet_rothalpy, continuity_residual, rothalpy_residual,
            throat.s - inlet.s, mach};
}

double air_viscosity_sutherland(double temperature_k)
{
    constexpr double reference_temperature_k = 273.15;
    constexpr double reference_viscosity_pa_s = 1.716e-5;
    constexpr double sutherland_temperature_k = 110.4;
    return reference_viscosity_pa_s *
        std::pow(temperature_k / reference_temperature_k, 1.5) *
        (reference_temperature_k + sutherland_temperature_k) /
        (temperature_k + sutherland_temperature_k);
}

double smooth_wall_darcy_factor(double reynolds)
{
    if (!(std::isfinite(reynolds) && reynolds > 0.0))
        throw std::invalid_argument("Invalid diffuser Reynolds number");
    if (reynolds < 2300.0)
        return 64.0 / reynolds;
    // Petukhov smooth-pipe estimate for turbulent flow. This is a
    // reference estimate, not a diffuser-specific calibrated correlation.
    const double denominator = 0.79 * std::log(reynolds) - 1.64;
    return 1.0 / (denominator * denominator);
}

struct ImpellerLossEstimate {
    double skin_friction_j_kg;
    double blade_loading_j_kg;
    double friction_factor;
    double diffusion_factor;
    double incidence_j_kg;
    double inlet_relative_angle_deg;
    double optimum_inlet_angle_deg;
    double entropy_generation_j_kg_k;
    double clearance_j_kg;
    double mixing_j_kg;
    double disk_friction_j_kg;
    double hydraulic_diameter_m;
    double reynolds;
    double blade_path_length_m;
    double mean_relative_velocity_m_s;
    double mean_density_kg_m3;
    double mean_temperature_k;
    double fanning_factor;
    double mixing_chi;
    double mixing_bstar;
    double correlation_exit_angle_deg;
    double inlet_mean_relative_velocity_m_s;
    double inlet_shroud_relative_velocity_m_s;
    double exit_relative_velocity_m_s;
    double throat_relative_velocity_m_s;
    double throat_density_kg_m3;
    double throat_enthalpy_j_kg;
    double throat_area_m2;
    double throat_blade_speed_m_s;
    double throat_rothalpy_j_kg;
    double throat_mass_residual_kg_s;
    double throat_rothalpy_residual_j_kg;
    double throat_entropy_residual_j_kg_k;
    double throat_mach;
    double choke_j_kg;
    double entrance_diffusion_j_kg;
    double recirculation_j_kg;
};

enum class LossClosure {
    MeanTemperatureEntropy,
    EnthalpyEquivalentEntropy
};

enum class HydraulicDiameterModel {
    Galvas,
    RectangularEndpointApproximation
};

struct LossModelOptions {
    bool skin_friction = true;
    bool blade_loading = true;
    bool incidence = true;
    bool clearance = true;
    bool mixing = true;
    bool disk_friction = true;
    LossClosure closure = LossClosure::EnthalpyEquivalentEntropy;
    HydraulicDiameterModel hydraulic_diameter_model =
        HydraulicDiameterModel::Galvas;
    double mixing_chi = 0.15;
    double mixing_bstar = 1.0;
    double skin_friction_multiplier = 1.0;
    double blade_loading_multiplier = 1.0;
    double incidence_multiplier = 1.0;
    double clearance_multiplier = 1.0;
    double mixing_multiplier = 1.0;
    double disk_friction_multiplier = 1.0;
    PublishedLossSet published_set = PublishedLossSet::CurrentScreen;
    bool use_shroud_inlet_speed_proxy = false;
    // Equation-form alternatives for auditing the empirical loss closures.
    // Defaults retain the existing baseline; each alternative is isolated in
    // the ablation table below.
    bool blade_loading_use_shroud_speed_consistently = false;
    double blade_loading_loading_coefficient = 0.75;
    bool incidence_use_mean_inlet_speed = false;
    bool clearance_use_exit_meridional_speed = false;
    bool disk_friction_use_diameter_clearance_ratio = false;
    // Diagnostic extension: integrate the local smooth-wall Cf(Re) relation
    // along a one-dimensional passage instead of evaluating it at mean state.
    bool distributed_skin_friction = false;
    // Broader passage diagnostic: interpolate endpoint rectangular Dh and
    // area, then obtain local meridional speed from continuity.
    bool geometry_continuity_skin_friction = false;
    // The source dataset gives beta2b=-30 deg; loss correlations use its
    // 30-deg magnitude. Keep this separate from the solver's complementary
    // 60-deg velocity-triangle angle.
    double correlation_exit_angle_deg = 30.0;
};

double schlichting_fanning_factor(double reynolds)
{
    if (!(std::isfinite(reynolds) && reynolds > 0.0))
        throw std::invalid_argument("Invalid impeller Reynolds number");
    // Schlichting correlation as reproduced by Meroni et al. The transition
    // band is linearly blended to avoid a discontinuity in this screening run.
    const double low = 2.67 * std::pow(reynolds, -0.5);
    const double high = 0.0622 * std::pow(reynolds, -0.2);
    if (reynolds <= 2.0e5)
        return low;
    if (reynolds >= 3.0e5)
        return high;
    const double weight = (reynolds - 2.0e5) / 1.0e5;
    return low * (1.0 - weight) + high * weight;
}

double distributed_skin_friction_loss(
    double path_length_m,
    double hydraulic_diameter_m,
    double inlet_relative_speed_m_s,
    double exit_relative_speed_m_s,
    double inlet_density_kg_m3,
    double exit_density_kg_m3,
    double inlet_temperature_k,
    double exit_temperature_k)
{
    // First-order passage screen: W, rho and T vary linearly from the
    // mean inlet state to the exit state; Dh remains the existing Galvas
    // equivalent diameter. Integrate the published differential form
    // d(h_sf)=2 Cf(Re) W^2 ds/Dh using composite Simpson quadrature.
    constexpr int intervals = 64;
    const double dx = 1.0 / intervals;
    double weighted_integrand_sum = 0.0;
    for (int i = 0; i <= intervals; ++i) {
        const double x = static_cast<double>(i) / intervals;
        const double speed = inlet_relative_speed_m_s + x *
            (exit_relative_speed_m_s - inlet_relative_speed_m_s);
        const double density = inlet_density_kg_m3 + x *
            (exit_density_kg_m3 - inlet_density_kg_m3);
        const double temperature = inlet_temperature_k + x *
            (exit_temperature_k - inlet_temperature_k);
        const double reynolds = density * speed * hydraulic_diameter_m /
            air_viscosity_sutherland(temperature);
        const double integrand = schlichting_fanning_factor(reynolds) *
            speed * speed;
        const double weight = (i == 0 || i == intervals)
            ? 1.0 : (i % 2 == 0 ? 2.0 : 4.0);
        weighted_integrand_sum += weight * integrand;
    }
    const double path_average = dx * weighted_integrand_sum / 3.0;
    return 2.0 * (path_length_m / hydraulic_diameter_m) * path_average;
}

double geometry_continuity_skin_friction_loss(
    double path_length_m,
    double inlet_hydraulic_diameter_m,
    double exit_hydraulic_diameter_m,
    double inlet_area_m2,
    double exit_area_m2,
    double mass_flow_kg_s,
    double inlet_density_kg_m3,
    double exit_density_kg_m3,
    double inlet_temperature_k,
    double exit_temperature_k,
    double inlet_radius_m,
    double exit_radius_m,
    double angular_speed_rad_s,
    double inlet_swirl_velocity_m_s,
    double exit_swirl_velocity_m_s)
{
    // Screening path: linearly interpolate the effective flow area,
    // rectangular-endpoint Dh, radius, thermodynamic state, and absolute
    // swirl. Continuity supplies Vm(x); the velocity triangle supplies W(x).
    // This is not an energy-closed throughflow solution.
    constexpr int intervals = 64;
    const double dx = 1.0 / intervals;
    double weighted_integrand_sum = 0.0;
    for (int i = 0; i <= intervals; ++i) {
        const double x = static_cast<double>(i) / intervals;
        const double area = inlet_area_m2 + x *
            (exit_area_m2 - inlet_area_m2);
        const double hydraulic_diameter = inlet_hydraulic_diameter_m + x *
            (exit_hydraulic_diameter_m - inlet_hydraulic_diameter_m);
        const double density = inlet_density_kg_m3 + x *
            (exit_density_kg_m3 - inlet_density_kg_m3);
        const double temperature = inlet_temperature_k + x *
            (exit_temperature_k - inlet_temperature_k);
        const double radius = inlet_radius_m + x *
            (exit_radius_m - inlet_radius_m);
        const double swirl = inlet_swirl_velocity_m_s + x *
            (exit_swirl_velocity_m_s - inlet_swirl_velocity_m_s);
        const double meridional_speed = mass_flow_kg_s / (density * area);
        const double blade_speed = angular_speed_rad_s * radius;
        const double relative_speed = std::hypot(
            meridional_speed, swirl - blade_speed);
        const double reynolds = density * relative_speed * hydraulic_diameter /
            air_viscosity_sutherland(temperature);
        const double integrand = schlichting_fanning_factor(reynolds) *
            relative_speed * relative_speed / hydraulic_diameter;
        const double weight = (i == 0 || i == intervals)
            ? 1.0 : (i % 2 == 0 ? 2.0 : 4.0);
        weighted_integrand_sum += weight * integrand;
    }
    const double path_average = dx * weighted_integrand_sum / 3.0;
    return 2.0 * path_length_m * path_average;
}

ImpellerLossEstimate estimate_impeller_losses(
    const MapPoint& point,
    const thermomodel::Model& air,
    const thermomodel::DifferentiableModel& derivatives,
    const ImpellerSystemGeometry& g,
    const OperatingPointRequest& request,
    const OperatingPointResult& trial,
    const OperatingThroatEstimate& throat,
    double throat_capacity_utilization,
    const LossModelOptions& options)
{
    constexpr double blades = 20.0;
    constexpr double impeller_axial_length_m = 0.130;
    constexpr double beta1_shroud_deg = 63.0;
    constexpr double beta1_hub_deg = 40.0;

    const double r1s = g.inlet.shroud_diameter / 2.0;
    const double r1h = g.inlet.hub_diameter / 2.0;
    const double r1m = g.impeller.inlet_diameter / 2.0;
    const double r2 = g.impeller.exit_diameter / 2.0;
    const double b2 = g.exit.width;
    const double pi = std::numbers::pi;
    const double omega = 2.0 * pi * point.rpm / 60.0;
    const double U1s = omega * r1s;
    const double U1h = omega * r1h;
    const double U2 = omega * r2;
    const double beta2 = g.impeller.exit_blade_angle_deg * pi / 180.0;
    const double beta2_correlation =
        options.correlation_exit_angle_deg * pi / 180.0;

    const auto inlet_state = evaluate(
        air, request.operating.inlet_temperature,
        request.operating.inlet_pressure);
    const auto outlet_state = evaluate(
        air, trial.outlet_temperature, trial.outlet_pressure);
    const double inlet_area = effective_flow_area(
        inlet_annulus_area(g.inlet), g.inlet.blockage);
    const double exit_area = effective_flow_area(
        impeller_exit_area(g.exit), g.exit.blockage);
    const double Vm1 = point.mass_flow / (inlet_state.rho * inlet_area);
    const double Vm2 = point.mass_flow / (outlet_state.rho * exit_area);
    const double W1s = std::hypot(Vm1, request.operating.inlet_prewhirl - U1s);
    const double W1h = std::hypot(Vm1, request.operating.inlet_prewhirl - U1h);
    const double W1 = std::sqrt(0.5 * (W1s * W1s + W1h * W1h));
    const auto throat_geometry = eckardt_meanline_throat(g);
    const auto throat_area = evaluate_throat_geometry(throat_geometry);
    const double inlet_area_ratio = inlet_area / throat_area.effective_area;
    // Incidence correlation uses the source's meridional angle convention;
    // the throat-area geometry above stores the complementary tangential angle.
    const double beta1_mean_meridional_deg =
        options.published_set == PublishedLossSet::CurrentScreen
            ? 0.5 * (63.0 + 40.0) : 56.0;
    const double beta1_blade = beta1_mean_meridional_deg * pi / 180.0;
    double beta1_opt = std::atan(
        inlet_area_ratio * std::tan(beta1_blade));
    if (options.published_set != PublishedLossSet::CurrentScreen) {
        // Conrad's inlet-metal-angle correction, using the documented mean
        // inlet blade thickness and the source's meridional angle convention.
        constexpr double inlet_blade_thickness_m = 0.00211;
        const double D1 = 2.0 * r1m;
        const double open_pitch_circumference =
            pi * D1 - blades * inlet_blade_thickness_m;
        if (!(open_pitch_circumference > 0.0))
            throw std::runtime_error("Invalid Conrad inlet pitch");
        beta1_opt = std::atan(
            (pi * D1 / open_pitch_circumference) * std::tan(beta1_blade));
    }
    // Meroni's inlet metal angles are referenced to the meridional axis.
    // Convert the relative-velocity triangle to that same convention.
    const double U1m = omega * r1m;
    const double W1m = std::hypot(
        Vm1, request.operating.inlet_prewhirl - U1m);
    const double Wu1m = request.operating.inlet_prewhirl - U1m;
    const double beta1_flow = std::atan2(
        U1m - request.operating.inlet_prewhirl, Vm1);
    const double incidence_angle = beta1_flow - beta1_opt;
    const double incidence_speed =
        options.published_set != PublishedLossSet::CurrentScreen
            ? W1
            : options.incidence_use_mean_inlet_speed ? W1m : W1s;
    const double incidence_sigma =
        options.published_set == PublishedLossSet::CurrentScreen ? 1.0 : 0.6;
    const double incidence_loss = incidence_sigma * 0.5 *
        incidence_speed * incidence_speed *
        std::pow(std::sin(incidence_angle), 2.0);
    const double Ctheta2 = request.operating.slip_factor *
        (U2 - Vm2 / std::tan(beta2));
    const double W2 = std::hypot(Vm2, Ctheta2 - U2);
    // Jansen's mean relative velocity used by the published skin-friction
    // correlation. Preserve the prior RMS/shroud proxy only as a diagnostic
    // alternative; the baseline follows the source definition.
    const double Wbar_jansen = (2.0 * W2 + W1h + W1s) / 4.0;
    const double Wbar_from_mean_inlet =
        std::sqrt(0.5 * (W1 * W1 + W2 * W2));
    const double Wbar_from_shroud =
        std::sqrt(0.5 * (W1s * W1s + W2 * W2));
    const double Wbar_aungier = std::sqrt(0.5 * (W1 * W1 + W2 * W2));
    const double Wbar = options.published_set == PublishedLossSet::EckardtA2042
        ? (2.0 * W2 + W1s - W1h) / 4.0
        : options.published_set == PublishedLossSet::CrossGeometry2137
        ? Wbar_aungier
        : options.use_shroud_inlet_speed_proxy
        ? std::max(Wbar_from_mean_inlet, Wbar_from_shroud)
        : Wbar_jansen;

    // Galvas hydraulic diameter: the sum of the outlet and inducer
    // contributions, using the published blade-angle and annulus geometry.
    const double beta1s = beta1_shroud_deg * pi / 180.0;
    const double beta1h = beta1_hub_deg * pi / 180.0;
    const double beta1_mean_cos = 0.5 * (std::cos(beta1s) + std::cos(beta1h));
    const double d1s = 2.0 * r1s;
    const double d1h = 2.0 * r1h;
    const double d2 = 2.0 * r2;
    const double dh_exit = d2 * std::cos(beta2_correlation) /
        (blades / pi + d2 * std::cos(beta2_correlation) / b2);
    const double dh_inlet = 0.5 * (d1s + d1h) * beta1_mean_cos /
        (blades / pi + (d1s + d1h) / (d1s - d1h) * beta1_mean_cos);
    const double galvas_hydraulic_diameter = dh_exit + dh_inlet;
    constexpr double inlet_blade_thickness_m = 0.00211;
    const double inlet_pitch = 2.0 * pi * r1m / blades -
        inlet_blade_thickness_m;
    const double exit_pitch = 2.0 * pi * r2 / blades *
        (1.0 - g.exit.blockage);
    if (!(inlet_pitch > 0.0 && exit_pitch > 0.0))
        throw std::runtime_error("Invalid impeller blade passage pitch");
    const double rectangular_dh1 = 2.0 * (r1s - r1h) * inlet_pitch /
        ((r1s - r1h) + inlet_pitch);
    const double rectangular_dh2 = 2.0 * b2 * exit_pitch /
        (b2 + exit_pitch);
    const double rectangular_hydraulic_diameter =
        0.5 * (rectangular_dh1 + rectangular_dh2);
    const double hydraulic_diameter =
        options.hydraulic_diameter_model == HydraulicDiameterModel::Galvas
            ? galvas_hydraulic_diameter
            : rectangular_hydraulic_diameter;

    // Jansen's impeller path-length approximation, using Eckardt A's
    // measured radii, axial length, and inlet/exit blade angles.
    const auto cos_deg = [](double degrees) {
        return std::cos(degrees * std::numbers::pi / 180.0);
    };
    const double angle_denominator =
        0.5 * (cos_deg(beta1_shroud_deg) + cos_deg(beta1_hub_deg)) +
        cos_deg(options.correlation_exit_angle_deg);
    const double blade_length = pi / 8.0 *
        (2.0 * r2 - (r1s + r1h) - b2 + 2.0 * impeller_axial_length_m) *
        (2.0 / angle_denominator);

    const double mean_temperature = 0.5 *
        (request.operating.inlet_temperature + trial.outlet_temperature);
    const double mean_density = 0.5 * (inlet_state.rho + outlet_state.rho);
    const double reynolds = mean_density * Wbar * hydraulic_diameter /
        air_viscosity_sutherland(mean_temperature);
    const double cf = schlichting_fanning_factor(reynolds);
    const double friction_factor = 4.0 * cf;
    // Published Jansen relation: dh_sf = 2 Cf (L_H / d_H) Wbar^2.
    // friction_factor remains the corresponding Darcy factor (4 Cf).
    const double skin_friction_loss = options.geometry_continuity_skin_friction
        ? geometry_continuity_skin_friction_loss(
              blade_length, rectangular_dh1, rectangular_dh2,
              inlet_area, exit_area, point.mass_flow,
              inlet_state.rho, outlet_state.rho,
              request.operating.inlet_temperature, trial.outlet_temperature,
              r1m, r2, omega, request.operating.inlet_prewhirl, Ctheta2)
        : options.distributed_skin_friction
        ? distributed_skin_friction_loss(
              blade_length, hydraulic_diameter,
              0.5 * (W1h + W1s), W2,
              inlet_state.rho, outlet_state.rho,
              request.operating.inlet_temperature, trial.outlet_temperature)
        : 2.0 * cf * (blade_length / hydraulic_diameter) * Wbar * Wbar;

    // Coppage diffusion-factor correlation for an impeller without splitters.
    const double diffusion_geometry =
        blades / pi * (1.0 - r1s / r2) + 2.0 * r1s / r2;
    const double work = std::max(0.0, trial.specific_work);
    const double W1_blade_loading =
        options.blade_loading_use_shroud_speed_consistently ? W1s : W1;
    const double diffusion_factor = 1.0 - W2 / W1_blade_loading +
        options.blade_loading_loading_coefficient *
            (work / (U2 * U2)) * (W2 / W1_blade_loading) /
            diffusion_geometry;
    const double delta_w_aungier =
        2.0 * pi * g.impeller.exit_diameter * Ctheta2 /
        (blades * blade_length);
    const double diffusion_factor_published = 1.0 - W2 / W1s +
        0.75 * (work / (U2 * U2)) * (W2 / W1s) /
            diffusion_geometry;
    const double blade_loading_loss =
        options.published_set == PublishedLossSet::CurrentScreen
            ? 0.05 * diffusion_factor * diffusion_factor * U2 * U2
            : delta_w_aungier * delta_w_aungier / 48.0;

    // Jansen tip-clearance leakage, with Eckardt A clearances from Meroni
    // Table 5: axial 0.235 mm and radial 0.190 mm.
    const double mean_clearance_m = 0.5 * (0.000235 + 0.000190);
    const double C1m = options.clearance_use_exit_meridional_speed
        ? Vm2 : Vm1;
    const double C2t = Ctheta2;
    const double jansen_clearance_loss = 0.6 * mean_clearance_m / b2 *
        std::abs(C2t) * std::sqrt(
            4.0 * pi / (b2 * blades) *
            ((r1s * r1s - r1h * r1h) /
             ((r2 - r1s) * (1.0 + outlet_state.rho / inlet_state.rho))) *
            std::abs(C2t) * C1m);
    const double clearance_loss =
        options.published_set == PublishedLossSet::EckardtA2042
            ? 0.1 * (0.0005 / b2) * U2 * U2
            : options.published_set == PublishedLossSet::CrossGeometry2137
            ? 0.6 * (0.0005 / b2) * std::abs(C2t) * std::sqrt(
                2.0 * pi * (std::pow(2.0 * r1s, 2.0) -
                            std::pow(2.0 * r1h, 2.0)) /
                (b2 * blades * (2.0 * r2 - 2.0 * r1s) *
                 (1.0 + outlet_state.rho / inlet_state.rho)) *
                std::abs(C2t) * Vm1)
            : jansen_clearance_loss;

    // Johnston-Dean jet/wake mixing. Meroni uses chi=0.15 and b*=1 for this
    // screening model. Here b* is diffuser-inlet depth / impeller-exit depth;
    // the sudden-expansion expression is stated for b*>1. Invert
    // chi=0.93 eps_w^2+0.07 eps_w to obtain the wake fraction.
    const double secondary_flow_fraction = options.mixing_chi;
    const double diffuser_to_impeller_width_ratio = options.mixing_bstar;
    const double wake_fraction =
        (-0.07 + std::sqrt(0.07 * 0.07 +
                           4.0 * 0.93 * secondary_flow_fraction)) /
        (2.0 * 0.93);
    const double Ctheta = Ctheta2;
    const double C2 = std::hypot(Vm2, Ctheta);
    const double alpha2 = std::atan2(Vm2, Ctheta);
    const double aungier_delta_w =
        2.0 * pi * g.impeller.exit_diameter * Ctheta2 /
        (blades * blade_length);
    // Aungier's mixing correlation uses mean-streamline W1xi to form Wmax,
    // then the inlet relative tangential component Wu1xi in Wout (Kovar et
    // al., 2021, Eqs. 43-45). Keep the shroud speed only for correlations
    // that explicitly call for W1s.
    const double w_max = 0.5 * (W1m + W2 + aungier_delta_w);
    const double deq = w_max / W2;
    const double w_separation = deq <= 2.0 ? W2 : W2 * deq / 2.0;
    const double aungier_w_out = std::hypot(
        Vm2 * exit_area / (pi * g.impeller.exit_diameter * b2), Wu1m);
    const double aungier_mixing_loss =
        0.5 * std::pow(w_separation - aungier_w_out, 2.0);
    const double mixing_loss =
        options.published_set == PublishedLossSet::CurrentScreen
            ? 0.5 * C2 * C2 *
                std::pow((1.0 - wake_fraction - diffuser_to_impeller_width_ratio) /
                         (1.0 - wake_fraction), 2.0) * std::cos(alpha2) *
                std::cos(alpha2)
            : aungier_mixing_loss;

    // Daily-Nece disk friction, using Eckardt A backface clearance 0.235 mm.
    const double mu2 = air_viscosity_sutherland(outlet_state.T);
    const double Re2 = outlet_state.rho * U2 * r2 / mu2;
    const double epsilon_backface_ratio =
        options.disk_friction_use_diameter_clearance_ratio
            ? 2.0 * 0.000235 / (2.0 * r2)
            : 0.000235 / b2;
    const double Kf = Re2 < 3.0e5
        ? 3.7 * std::pow(epsilon_backface_ratio, 0.1) /
              std::sqrt(Re2)
        : 0.102 * std::pow(epsilon_backface_ratio, 0.1) /
              std::pow(Re2, 0.2);
    const double shepherd_reynolds = outlet_state.rho * U2 * (2.0 * r2) / mu2;
    const double disk_loss =
        options.published_set == PublishedLossSet::CurrentScreen
            ? 0.25 * (0.5 * (inlet_state.rho + outlet_state.rho)) *
                U2 * U2 * U2 * r2 * r2 * Kf / point.mass_flow
            : 0.01356 * outlet_state.rho * U2 * U2 * U2 *
                (2.0 * r2) * (2.0 * r2) /
                (point.mass_flow * std::pow(shepherd_reynolds, 0.2));

    double choke_loss = 0.0;
    double entrance_diffusion_loss = 0.0;
    double recirculation_loss = 0.0;
    if (options.published_set != PublishedLossSet::CurrentScreen) {
        // Aungier choke closure. The source's corrected throat-area ratio is
        // unavailable in this mean-line reconstruction, so capacity
        // utilization is used as a documented proxy for A_th*/(Cr A_th).
        const double sonic_area_ratio_proxy = 1.0 /
            std::max(1.0e-12, throat_capacity_utilization);
        const double x_choke = std::max(
            0.0, 10.0 * (1.1 - sonic_area_ratio_proxy));
        choke_loss = 0.5 * W1 * W1 *
            (0.05 * x_choke + std::pow(x_choke, 7.0));
        // Aungier entrance diffusion/stall closure, using the coupled
        // operating-throat relative speed as W_th.
        const double base_entrance = std::max(
            0.0, 0.4 * std::pow(W1 - throat.relative_speed_m_s, 2.0) -
                     incidence_loss);
        entrance_diffusion_loss = base_entrance;
        if (throat.relative_speed_m_s > 0.0 &&
            W1s / throat.relative_speed_m_s > 1.75) {
            entrance_diffusion_loss = std::min(
                base_entrance,
                std::max(0.0, 0.5 * std::pow(
                    W1s - 1.75 * throat.relative_speed_m_s, 2.0) -
                    incidence_loss));
        }
        recirculation_loss = 0.02 *
            std::sqrt(std::max(0.0, std::tan(alpha2))) *
            diffusion_factor_published * diffusion_factor_published * U2 * U2;
    }

    const double internal_loss =
        (options.skin_friction
             ? options.skin_friction_multiplier * skin_friction_loss : 0.0) +
        (options.blade_loading
             ? options.blade_loading_multiplier * blade_loading_loss : 0.0) +
        (options.incidence
             ? options.incidence_multiplier * incidence_loss : 0.0) +
        (options.clearance
             ? options.clearance_multiplier * clearance_loss : 0.0) +
        choke_loss + entrance_diffusion_loss +
        (options.published_set == PublishedLossSet::CurrentScreen
             ? 0.0 : options.mixing_multiplier * mixing_loss);
    double entropy_generation = internal_loss / mean_temperature;
    if (options.closure == LossClosure::EnthalpyEquivalentEntropy) {
        // Convert the loss enthalpy to entropy at the trial outlet total
        // pressure using the selected property model, instead of Δh/T_mean.
        const FlowStation inlet_station{
            inlet_state, {Vm1, request.operating.inlet_prewhirl}};
        const FlowStation outlet_station{
            outlet_state, {Vm2, Ctheta2}};
        const auto total_in = evaluate_stagnation_state(air, inlet_station)
                                  .stagnation_state;
        const auto total_out = evaluate_stagnation_state(air, outlet_station)
                                   .stagnation_state;
        const auto isentropic_out = invert_temperature(
            air, total_out.p, total_in.s,
            Property::Entropy, total_in.T).state;
        if (!(std::isfinite(isentropic_out.h) && internal_loss >= 0.0))
            throw std::runtime_error("Invalid enthalpy-equivalent loss state");
        const auto loss_state = invert_temperature(
            air, total_out.p, isentropic_out.h + internal_loss,
            Property::Enthalpy, isentropic_out.T).state;
        entropy_generation = internal_loss == 0.0
            ? 0.0 : loss_state.s - total_in.s;
    }
    if (!(std::isfinite(blade_length) && blade_length > 0.0 &&
          std::isfinite(hydraulic_diameter) && hydraulic_diameter > 0.0 &&
          std::isfinite(skin_friction_loss) && skin_friction_loss >= 0.0 &&
          std::isfinite(blade_loading_loss) && blade_loading_loss >= 0.0 &&
          std::isfinite(incidence_loss) && incidence_loss >= 0.0 &&
          std::isfinite(clearance_loss) && clearance_loss >= 0.0 &&
          std::isfinite(mixing_loss) && mixing_loss >= 0.0 &&
          std::isfinite(disk_loss) && disk_loss >= 0.0 &&
          std::isfinite(choke_loss) && choke_loss >= 0.0 &&
          std::isfinite(entrance_diffusion_loss) &&
              entrance_diffusion_loss >= 0.0 &&
          std::isfinite(recirculation_loss) && recirculation_loss >= 0.0 &&
          std::isfinite(entropy_generation) && entropy_generation >= 0.0))
        throw std::runtime_error("Invalid impeller loss estimate");

    return {skin_friction_loss, blade_loading_loss, friction_factor,
            diffusion_factor, incidence_loss,
            beta1_flow * 180.0 / pi, beta1_opt * 180.0 / pi,
            entropy_generation, clearance_loss, mixing_loss, disk_loss,
            hydraulic_diameter, reynolds, blade_length, Wbar,
            mean_density, mean_temperature,
            cf, secondary_flow_fraction,
            diffuser_to_impeller_width_ratio,
            options.correlation_exit_angle_deg,
            W1, W1s, W2, throat.relative_speed_m_s,
            throat.density_kg_m3, throat.enthalpy_j_kg,
            throat.effective_area_m2,
            throat.blade_speed_m_s, throat.inlet_rothalpy_j_kg,
            throat.mass_residual_kg_s, throat.rothalpy_residual_j_kg,
            throat.entropy_residual_j_kg_k, throat.mach,
            choke_loss, entrance_diffusion_loss, recirculation_loss};
}

Prediction predict(
    const MapPoint& point,
    const thermomodel::Model& air,
    const thermomodel::DifferentiableModel& derivatives,
    const ImpellerSystemGeometry& g,
    double darcy_factor,
    bool apply_galvas_incidence,
    ThroatCapacityCache& throat_capacity_cache,
    bool include_skin_friction = true,
    bool include_blade_loading = true,
    LossModelOptions loss_options = {},
    double diffuser_friction_multiplier = 1.0)
{
    loss_options.skin_friction = include_skin_friction;
    loss_options.blade_loading = include_blade_loading;
    loss_options.incidence = apply_galvas_incidence;
    auto request = make_request(point, air, derivatives, g);
    const auto throat_capacity = estimate_throat_capacity(
        air, derivatives, g, request, throat_capacity_cache);
    const auto operating_throat = estimate_operating_throat(
        air, derivatives, g, request);
    NonlinearSolverOptions solver_options;
    solver_options.residual_tolerance = 1.0e-12;

    auto impeller = evaluate_operating_point(
        air, g, request, scales(), solver_options);
    if (!(impeller.converged && impeller.independently_validated &&
          impeller.basic_physics_valid)) {
        throw std::runtime_error("Impeller solve failed: " + impeller.message);
    }

    ImpellerLossEstimate impeller_losses{};
    for (int iteration = 0; iteration < 20; ++iteration) {
        impeller_losses = estimate_impeller_losses(
            point, air, derivatives, g, request, impeller,
            operating_throat,
            point.mass_flow / throat_capacity.critical_mass_flow_kg_s,
            loss_options);
        request.operating.entropy_generation =
            impeller_losses.entropy_generation_j_kg_k;
        const auto corrected = evaluate_operating_point(
            air, g, request, scales(), solver_options);
        if (!(corrected.converged && corrected.independently_validated &&
              corrected.basic_physics_valid)) {
            throw std::runtime_error(
                "Loss-corrected impeller solve failed: " + corrected.message);
        }
        const double relative_change = std::abs(
            corrected.outlet_pressure - impeller.outlet_pressure) /
            std::max(1.0, std::abs(corrected.outlet_pressure));
        impeller = corrected;
        if (relative_change < 1.0e-10)
            break;
    }
    // Keep diagnostics consistent with the final loss-corrected state.
    impeller_losses = estimate_impeller_losses(
        point, air, derivatives, g, request, impeller, operating_throat,
        point.mass_flow / throat_capacity.critical_mass_flow_kg_s,
        loss_options);

    // Mixing and backface disk friction are external losses in the cited
    // mean-line formulation: apply their equivalent enthalpy rise at the
    // impeller exit pressure before solving the vaneless diffuser.
    const double external_loss =
        (loss_options.published_set == PublishedLossSet::CurrentScreen &&
         loss_options.mixing ? loss_options.mixing_multiplier *
             impeller_losses.mixing_j_kg : 0.0) +
        (loss_options.disk_friction ? loss_options.disk_friction_multiplier *
             impeller_losses.disk_friction_j_kg : 0.0) +
        impeller_losses.recirculation_j_kg;
    const auto impeller_exit_static = evaluate(
        air, impeller.outlet_temperature, impeller.outlet_pressure);
    const double outlet_static_temperature = invert_temperature(
        air, impeller.outlet_pressure,
        impeller_exit_static.h + external_loss,
        Property::Enthalpy, impeller.outlet_temperature).state.T;

    const double mass = request.operating.mass_flow;
    const double r2 = g.exit.diameter / 2.0;
    const double r3 = 2.69 * r2;
    const double b2 = g.exit.width;
    const double b3 = 0.51 * b2;
    const double area2 = effective_flow_area(
        impeller_exit_area(g.exit), g.exit.blockage);
    const double width2 = area2 / (2.0 * std::numbers::pi * r2);

    const auto state1 = evaluate(
        air, request.operating.inlet_temperature,
        request.operating.inlet_pressure);
    const auto state2 = evaluate(
        air, outlet_static_temperature, impeller.outlet_pressure);
    const double vm2 = mass / (state2.rho * area2);
    const double omega = 2.0 * std::numbers::pi *
        request.operating.rotational_speed_rpm / 60.0;
    const double U2 = omega * g.impeller.exit_diameter / 2.0;
    const double beta = g.impeller.exit_blade_angle_deg *
        std::numbers::pi / 180.0;
    const double whirl2 = request.operating.slip_factor *
        (U2 - vm2 / std::tan(beta));

    const double inlet_speed = std::hypot(vm2, whirl2);
    const double hydraulic_diameter = 2.0 * width2;
    const double reynolds = state2.rho * inlet_speed *
        hydraulic_diameter / air_viscosity_sutherland(state2.T);
    double used_darcy_factor = (darcy_factor < 0.0
        ? smooth_wall_darcy_factor(reynolds)
        : darcy_factor) * diffuser_friction_multiplier;

    const VanelessDiffuserInlet inlet{state2, vm2, whirl2, mass};
    const VariableWidthDiffuserGeometry diffuser_geometry{
        r2,
        r3,
        [=](double r) {
            return width2 + (b3 - width2) * (r - r2) / (r3 - r2);
        },
        [=](double) { return (b3 - width2) / (r3 - r2); }
    };
    FrictionalDiffuserOptions diffuser_options;
    diffuser_options.darcy_friction_factor = used_darcy_factor;
    diffuser_options.radial_steps = 1024;

    if (loss_options.published_set != PublishedLossSet::CurrentScreen &&
        darcy_factor < 0.0) {
        // Map the literature's vaneless-diffuser loss to an equivalent Darcy
        // factor while retaining the existing variable-width integration.
        // Coppage gives an explicit loss. Stanitz is evaluated on the current
        // smooth-wall reference solution, then its entropy target is matched.
        const FlowStation diffuser_inlet_station{state2, {vm2, whirl2}};
        const auto total_in = evaluate_stagnation_state(
            air, diffuser_inlet_station).stagnation_state;
        const double C2 = std::hypot(vm2, whirl2);
        const double alpha2 = std::atan2(vm2, whirl2);
        double target_loss = 0.0;
        if (loss_options.published_set == PublishedLossSet::EckardtA2042) {
            const double cf_diffuser = schlichting_fanning_factor(reynolds);
            target_loss = (2.0 / 3.0) * std::pow(r2 / r3, 1.5) *
                (cf_diffuser * (2.0 * r2) / (8.0 * width2)) *
                C2 * C2 / std::pow(std::cos(alpha2), 2.0);
        } else {
            const auto reference_diffuser =
                evaluate_variable_width_vaneless_diffuser(
                    air, derivatives, diffuser_geometry, inlet, diffuser_options);
            const FlowStation outlet_station{
                reference_diffuser.outlet_static,
                {reference_diffuser.outlet_radial_velocity,
                 reference_diffuser.outlet_tangential_velocity}};
            const auto total_out = evaluate_stagnation_state(
                air, outlet_station).stagnation_state;
            const double gamma = 1004.5 / (1004.5 -
                8.31446261815324 / 0.0289652);
            const double exponent = (gamma - 1.0) / gamma;
            target_loss = 1004.5 * total_in.T *
                (std::pow(reference_diffuser.outlet_static.p / total_out.p,
                          exponent) -
                 std::pow(reference_diffuser.outlet_static.p / total_in.p,
                          exponent));
            target_loss = std::max(0.0, target_loss);
        }
        const double target_entropy = target_loss /
            std::max(1.0, total_in.T);
        auto entropy_residual = [&](double factor) {
            auto trial_options = diffuser_options;
            trial_options.darcy_friction_factor = factor;
            const auto trial_diffuser = evaluate_variable_width_vaneless_diffuser(
                air, derivatives, diffuser_geometry, inlet, trial_options);
            return trial_diffuser.integrated_entropy_production - target_entropy;
        };
        double lo = 0.0;
        double hi = 0.2;
        while (entropy_residual(hi) < 0.0 && hi < 10.0)
            hi *= 2.0;
        for (int i = 0; i < 32; ++i) {
            const double mid = 0.5 * (lo + hi);
            if (entropy_residual(mid) < 0.0)
                lo = mid;
            else
                hi = mid;
        }
        used_darcy_factor = 0.5 * (lo + hi);
        diffuser_options.darcy_friction_factor = used_darcy_factor;
    }
    const auto diffuser = evaluate_variable_width_vaneless_diffuser(
        air, derivatives, diffuser_geometry, inlet, diffuser_options);

    // Japikse reports diffuser traverses at four radii, including r/r2=1.075
    // and 1.6875. Integrate the model to these planes for future comparison
    // with digitized measurements. The static pressure below is a bulk station
    // value; the model does not predict the corresponding wall-tap pressure.
    const auto evaluate_traverse = [&](double radius_ratio) {
        const double radius = radius_ratio * r2;
        const double width_at_radius = width2 +
            (b3 - width2) * (radius - r2) / (r3 - r2);
        const VariableWidthDiffuserGeometry station_geometry{
            r2,
            radius,
            [=](double r) {
                return width2 + (b3 - width2) * (r - r2) / (r3 - r2);
            },
            [=](double) { return (b3 - width2) / (r3 - r2); }
        };
        const auto station = evaluate_variable_width_vaneless_diffuser(
            air, derivatives, station_geometry, inlet, diffuser_options);
        const FlowStation flow_station{
            station.outlet_static,
            {station.outlet_radial_velocity,
             station.outlet_tangential_velocity}
        };
        const auto total = evaluate_stagnation_state(air, flow_station)
            .stagnation_state;
        return DiffuserTraverseState{
            impeller.outlet_pressure,
            radius_ratio,
            radius,
            width_at_radius,
            station.outlet_static.p,
            total.p,
            total.T,
            std::atan2(station.outlet_tangential_velocity,
                       station.outlet_radial_velocity) *
                180.0 / std::numbers::pi,
            station.outlet_radial_velocity,
            station.outlet_tangential_velocity
        };
    };
    const auto near_traverse = evaluate_traverse(1.075);
    const auto far_traverse = evaluate_traverse(1.6875);

    const double vm1 = mass /
        (state1.rho * effective_flow_area(
            inlet_annulus_area(g.inlet), g.inlet.blockage));
    const FlowStation station1{state1, {vm1, request.operating.inlet_prewhirl}};
    const FlowStation station2{state2, {vm2, whirl2}};
    const FlowStation station3{
        diffuser.outlet_static,
        {diffuser.outlet_radial_velocity,
         diffuser.outlet_tangential_velocity}
    };

    const double p01 = evaluate_stagnation_state(air, station1)
        .stagnation_state.p;
    const double p02 = evaluate_stagnation_state(air, station2)
        .stagnation_state.p;
    const double p03 = evaluate_stagnation_state(air, station3)
        .stagnation_state.p;
    const auto total1 = evaluate_stagnation_state(air, station1).stagnation_state;
    const auto total3 = evaluate_stagnation_state(air, station3).stagnation_state;
    const auto total3_isentropic = invert_temperature(
        air, p03, total1.s, Property::Entropy, total1.T).state;
    const double eta_tt = (total3_isentropic.h - total1.h) /
        (total3.h - total1.h);

    const double capacity_utilization = point.mass_flow /
        throat_capacity.critical_mass_flow_kg_s;
    const double capacity_margin_percent =
        100.0 * (1.0 - capacity_utilization);

    return {p03 / p01, p02 / p01, p03 / p02,
            100.0 * (1.0 - p03 / p02), used_darcy_factor, reynolds,
            eta_tt, impeller_losses.skin_friction_j_kg,
            impeller_losses.blade_loading_j_kg,
            impeller_losses.friction_factor,
            impeller_losses.diffusion_factor,
            impeller_losses.incidence_j_kg,
            impeller_losses.inlet_relative_angle_deg,
            impeller_losses.optimum_inlet_angle_deg,
            throat_capacity.effective_area_m2,
            throat_capacity.critical_mass_flow_kg_s,
            capacity_utilization,
            capacity_margin_percent,
            capacity_margin_percent <= 0.0,
            impeller.specific_work +
                (loss_options.disk_friction
                     ? impeller_losses.disk_friction_j_kg : 0.0) +
                impeller_losses.recirculation_j_kg,
            total3.T,
            impeller_losses.reynolds,
            impeller_losses.clearance_j_kg,
            impeller_losses.mixing_j_kg,
            impeller_losses.disk_friction_j_kg,
            impeller_losses.hydraulic_diameter_m,
            impeller_losses.entropy_generation_j_kg_k,
            impeller_losses.blade_path_length_m,
            impeller_losses.mean_relative_velocity_m_s,
            impeller_losses.mean_density_kg_m3,
            impeller_losses.mean_temperature_k,
            impeller_losses.fanning_factor,
            impeller_losses.mixing_chi,
            impeller_losses.mixing_bstar,
            impeller_losses.correlation_exit_angle_deg,
            impeller_losses.inlet_mean_relative_velocity_m_s,
            impeller_losses.inlet_shroud_relative_velocity_m_s,
            impeller_losses.exit_relative_velocity_m_s,
            impeller_losses.throat_relative_velocity_m_s,
            impeller_losses.throat_density_kg_m3,
            impeller_losses.throat_enthalpy_j_kg,
            impeller_losses.throat_area_m2,
            impeller_losses.throat_blade_speed_m_s,
            impeller_losses.throat_rothalpy_j_kg,
            impeller_losses.throat_mass_residual_kg_s,
            impeller_losses.throat_rothalpy_residual_j_kg,
            impeller_losses.throat_entropy_residual_j_kg_k,
            impeller_losses.throat_mach,
            near_traverse,
            far_traverse,
            impeller_losses.choke_j_kg,
            impeller_losses.entrance_diffusion_j_kg,
            impeller_losses.recirculation_j_kg,
            total3.h - total1.h,
            impeller.specific_work,
            external_loss,
            (total3.h - total1.h) - impeller.specific_work,
            (total3.h - total1.h) - impeller.specific_work -
                (loss_options.disk_friction
                     ? impeller_losses.disk_friction_j_kg : 0.0) -
                impeller_losses.recirculation_j_kg,
            (total3.h - total1.h) - impeller.specific_work - external_loss};
}

} // namespace

#ifndef TURBOCOMP_ECKARDT_NO_MAIN
int main(int argc, char** argv)
{
    try {
        const std::filesystem::path input = argc > 1
            ? argv[1] : "data/validation/eckardt_impeller_a.csv";
        const std::filesystem::path output = argc > 2
            ? argv[2]
            : "benchmarks/results-eckardt-impeller-a/diffuser_friction_sweep.csv";
        const std::optional<std::filesystem::path> efficiency_input = argc > 3
            ? std::optional<std::filesystem::path>{argv[3]} : std::nullopt;
        const std::optional<std::filesystem::path> efficiency_output = argc > 4
            ? std::optional<std::filesystem::path>{argv[4]} : std::nullopt;
        const std::string property_model = argc > 5 ? argv[5] : "ideal-gas";
        if (property_model != "ideal-gas" && property_model != "pr-air")
            throw std::invalid_argument(
                "Thermodynamic model must be ideal-gas or pr-air");
        const auto points = read_map(input);
        const auto g = geometry();
        const thermomodel::IdealGas ideal_air;
        thermomodel::PRParameters pseudo_air_parameters;
        pseudo_air_parameters.molar_mass = 0.0289652;
        pseudo_air_parameters.Tc = 132.53;
        pseudo_air_parameters.pc = 3.786e6;
        pseudo_air_parameters.acentric = 0.035;
        pseudo_air_parameters.Tref = 298.15;
        pseudo_air_parameters.pref = 101325.0;
        const auto air_caloric = std::make_shared<const thermomodel::ConstantCp>(
            1004.5);
        const thermomodel::PengRobinson pr_air(
            air_caloric, pseudo_air_parameters);
        const thermomodel::Model& air = property_model == "pr-air"
            ? static_cast<const thermomodel::Model&>(pr_air)
            : static_cast<const thermomodel::Model&>(ideal_air);
        const thermomodel::DifferentiableModel& derivatives =
            property_model == "pr-air"
                ? static_cast<const thermomodel::DifferentiableModel&>(pr_air)
                : static_cast<const thermomodel::DifferentiableModel&>(ideal_air);
        ThroatCapacityCache throat_capacity_cache;
        std::filesystem::create_directories(output.parent_path());
        std::ofstream csv(output);
        if (!csv)
            throw std::runtime_error("Cannot write: " + output.string());
        std::ofstream model_record(output.parent_path() /
            "thermodynamic_model.txt");
        if (!model_record)
            throw std::runtime_error("Cannot write thermodynamic model record");
        model_record << "property_model=" << property_model << '\n'
                     << "published_loss_sets=diagnostic scenarios in loss_term_ablation.csv\n";
        if (property_model == "pr-air")
            model_record << "eos=Peng-Robinson pseudo-pure dry air\n"
                            "molar_mass_kg_mol=0.0289652\n"
                            "critical_temperature_K=132.53\n"
                            "critical_pressure_Pa=3786000\n"
                            "acentric_factor=0.035\n"
                            "ideal_cp_J_kg_K=1004.5\n";
        std::cout << "Thermodynamic model: " << property_model << '\n';
        csv << std::setprecision(12)
            << "model,darcy_factor,inlet_reynolds,rpm,mass_flow_kg_s,"
               "measured_total_pressure_ratio,predicted_total_pressure_ratio,"
               "predicted_impeller_total_pressure_ratio,"
               "diffuser_total_pressure_ratio,diffuser_loss_percent,"
               "relative_error_percent,predicted_eta_tt,"
               "impeller_skin_friction_loss_j_kg,"
               "impeller_blade_loading_loss_j_kg,"
               "impeller_friction_factor,impeller_diffusion_factor,"
               "impeller_incidence_loss_j_kg,inlet_relative_angle_from_meridional_deg,"
               "optimum_inlet_angle_deg,estimated_throat_area_m2,"
               "critical_throat_flow_kg_s,throat_capacity_utilization,"
               "throat_capacity_margin_percent,choke_screen_flag\n";
        const std::vector<double> friction_factors{
            0.0, 0.005, 0.010, 0.015, 0.020, 0.030, 0.040
        };

        std::cout << "Linear diffuser-width taper: b2_eff="
                  << std::setprecision(6)
                  << effective_flow_area(impeller_exit_area(g.exit),
                                         g.exit.blockage) /
                         (std::numbers::pi * g.exit.diameter)
                  << " m to b3=" << 0.51 * g.exit.width << " m; r3/r2=2.69\n"
                  << "Darcy factor, complete, bias_pct, MAPE_pct, max_abs_pct, "
                     "mean_diffuser_loss_pct\n";

        for (const double f : friction_factors) {
            double bias = 0.0;
            double mape = 0.0;
            double max_abs = 0.0;
            double mean_loss = 0.0;
            int complete = 0;
            for (const auto& point : points) {
                const auto result = predict(
                    point, air, derivatives, g, f, true, throat_capacity_cache);
                const double error = 100.0 *
                    (result.stage_pr - point.measured_pressure_ratio) /
                    point.measured_pressure_ratio;
                csv << "fixed_galvas_incidence," << f << ',' << result.inlet_reynolds << ','
                    << point.rpm << ',' << point.mass_flow << ','
                    << point.measured_pressure_ratio << ',' << result.stage_pr
                    << ',' << result.impeller_pr << ',' << result.diffuser_pr
                    << ',' << result.diffuser_loss_percent << ',' << error
                    << ',' << result.total_to_total_efficiency
                    << ',' << result.impeller_skin_friction_loss_j_kg
                    << ',' << result.impeller_blade_loading_loss_j_kg
                    << ',' << result.impeller_skin_friction_factor
                    << ',' << result.impeller_diffusion_factor
                    << ',' << result.impeller_incidence_loss_j_kg
                    << ',' << result.inlet_relative_angle_deg
                    << ',' << result.optimum_inlet_angle_deg
                    << ',' << result.estimated_throat_area_m2
                    << ',' << result.critical_throat_flow_kg_s
                    << ',' << result.throat_capacity_utilization
                    << ',' << result.throat_capacity_margin_percent
                    << ',' << (result.choke_screen_flag ? 1 : 0)
                    << '\n';
                bias += error;
                mape += std::abs(error);
                max_abs = std::max(max_abs, std::abs(error));
                mean_loss += result.diffuser_loss_percent;
                ++complete;
            }
            const double n = static_cast<double>(complete);
            std::cout << f << ", " << complete << ", " << bias / n
                      << ", " << mape / n << ", " << max_abs
                      << ", " << mean_loss / n << '\n';
        }

        double smooth_bias = 0.0;
        double smooth_mape = 0.0;
        double smooth_max_abs = 0.0;
        double smooth_mean_f = 0.0;
        double smooth_min_f = std::numeric_limits<double>::infinity();
        double smooth_max_f = 0.0;
        double smooth_mean_re = 0.0;
        for (const auto& point : points) {
            const auto result = predict(
                point, air, derivatives, g, -1.0, true, throat_capacity_cache);
            const double error = 100.0 *
                (result.stage_pr - point.measured_pressure_ratio) /
                point.measured_pressure_ratio;
            csv << "smooth_wall_inlet_Re_galvas_incidence," << result.darcy_factor << ','
                << result.inlet_reynolds << ',' << point.rpm << ','
                << point.mass_flow << ',' << point.measured_pressure_ratio
                << ',' << result.stage_pr << ',' << result.impeller_pr << ','
                << result.diffuser_pr << ',' << result.diffuser_loss_percent
                << ',' << error << ',' << result.total_to_total_efficiency
                << ',' << result.impeller_skin_friction_loss_j_kg
                << ',' << result.impeller_blade_loading_loss_j_kg
                << ',' << result.impeller_skin_friction_factor
                << ',' << result.impeller_diffusion_factor
                << ',' << result.impeller_incidence_loss_j_kg
                << ',' << result.inlet_relative_angle_deg
                << ',' << result.optimum_inlet_angle_deg
                << ',' << result.estimated_throat_area_m2
                << ',' << result.critical_throat_flow_kg_s
                << ',' << result.throat_capacity_utilization
                << ',' << result.throat_capacity_margin_percent
                << ',' << (result.choke_screen_flag ? 1 : 0) << '\n';
            smooth_bias += error;
            smooth_mape += std::abs(error);
            smooth_max_abs = std::max(smooth_max_abs, std::abs(error));
            smooth_mean_f += result.darcy_factor;
            smooth_min_f = std::min(smooth_min_f, result.darcy_factor);
            smooth_max_f = std::max(smooth_max_f, result.darcy_factor);
            smooth_mean_re += result.inlet_reynolds;
        }
        double no_incidence_bias = 0.0;
        double no_incidence_mape = 0.0;
        double no_incidence_max_abs = 0.0;
        for (const auto& point : points) {
            const auto result = predict(
                point, air, derivatives, g, -1.0, false, throat_capacity_cache);
            const double error = 100.0 *
                (result.stage_pr - point.measured_pressure_ratio) /
                point.measured_pressure_ratio;
            csv << "smooth_wall_inlet_Re_no_incidence," << result.darcy_factor
                << ',' << result.inlet_reynolds << ',' << point.rpm << ','
                << point.mass_flow << ',' << point.measured_pressure_ratio
                << ',' << result.stage_pr << ',' << result.impeller_pr << ','
                << result.diffuser_pr << ',' << result.diffuser_loss_percent
                << ',' << error << ',' << result.total_to_total_efficiency
                << ',' << result.impeller_skin_friction_loss_j_kg
                << ',' << result.impeller_blade_loading_loss_j_kg
                << ',' << result.impeller_skin_friction_factor
                << ',' << result.impeller_diffusion_factor
                << ',' << result.impeller_incidence_loss_j_kg
                << ',' << result.inlet_relative_angle_deg
                << ',' << result.optimum_inlet_angle_deg
                << ',' << result.estimated_throat_area_m2
                << ',' << result.critical_throat_flow_kg_s
                << ',' << result.throat_capacity_utilization
                << ',' << result.throat_capacity_margin_percent
                << ',' << (result.choke_screen_flag ? 1 : 0) << '\n';
            no_incidence_bias += error;
            no_incidence_mape += std::abs(error);
            no_incidence_max_abs = std::max(
                no_incidence_max_abs, std::abs(error));
        }
        const double count = static_cast<double>(points.size());
        std::cout << "smooth-wall inlet-Re estimate: mean_Re="
                  << smooth_mean_re / count
                  << ", Darcy_factor_mean=" << smooth_mean_f / count
                  << ", range=[" << smooth_min_f << ", " << smooth_max_f
                  << "], bias_pct=" << smooth_bias / count
                  << ", MAPE_pct=" << smooth_mape / count
                  << ", max_abs_pct=" << smooth_max_abs << '\n';
        std::cout << "smooth-wall without Galvas incidence: bias_pct="
                  << no_incidence_bias / count
                  << ", MAPE_pct=" << no_incidence_mape / count
                  << ", max_abs_pct=" << no_incidence_max_abs << '\n';
        std::cout << "Per-point CSV: " << output << '\n';
        const auto traverse_output = output.parent_path() /
            "impeller_a_diffuser_traverse_predictions.csv";
        std::ofstream traverse_csv(traverse_output);
        if (!traverse_csv)
            throw std::runtime_error(
                "Cannot write: " + traverse_output.string());
        traverse_csv << std::setprecision(12)
            << "model,model_geometry,rpm,mass_flow_kg_s,station_id,r_over_r2,radius_m,"
               "model_width_m,predicted_impeller_exit_static_pressure_pa,"
               "predicted_bulk_static_pressure_pa,"
               "predicted_total_pressure_pa,predicted_total_temperature_k,"
               "predicted_flow_angle_from_meridional_deg,"
               "predicted_radial_velocity_m_s,"
               "predicted_tangential_velocity_m_s\n";
        for (const auto& point : points) {
            const auto result = predict(
                point, air, derivatives, g, -1.0, true, throat_capacity_cache);
            const auto write_traverse = [&](
                const char* station_id,
                const DiffuserTraverseState& station) {
                traverse_csv << "smooth_wall_inlet_Re_galvas_incidence,"
                    << "eckardt_impeller_a_backswept,"
                    << point.rpm << ',' << point.mass_flow << ',' << station_id
                    << ',' << station.radius_ratio << ',' << station.radius_m
                    << ',' << station.width_m
                    << ',' << station.impeller_exit_static_pressure_pa
                    << ',' << station.bulk_static_pressure_pa
                    << ',' << station.total_pressure_pa
                    << ',' << station.total_temperature_k
                    << ',' << station.flow_angle_from_meridional_deg
                    << ',' << station.radial_velocity_m_s
                    << ',' << station.tangential_velocity_m_s << '\n';
            };
            write_traverse("near_1.075", result.near_traverse);
            write_traverse("far_1.6875", result.far_traverse);
        }
        std::cout << "Diffuser traverse predictions CSV: "
                  << traverse_output << '\n';
        int pressure_choke_flags = 0;
        double minimum_pressure_capacity_margin =
            std::numeric_limits<double>::infinity();
        for (const auto& point : points) {
            const auto result = predict(
                point, air, derivatives, g, -1.0, true, throat_capacity_cache);
            pressure_choke_flags += result.choke_screen_flag ? 1 : 0;
            minimum_pressure_capacity_margin = std::min(
                minimum_pressure_capacity_margin,
                result.throat_capacity_margin_percent);
        }
        std::cout << "Estimated meanline inlet-throat screen: "
                  << pressure_choke_flags << '/' << points.size()
                  << " pressure points exceed capacity; minimum margin_pct="
                  << minimum_pressure_capacity_margin << '\n';

        if (efficiency_input) {
            if (!efficiency_output)
                throw std::invalid_argument(
                    "Efficiency CSV output path is required when efficiency input is given");
            const auto efficiency_points = read_efficiency_map(*efficiency_input);
            std::ofstream efficiency_csv(*efficiency_output);
            if (!efficiency_csv)
                throw std::runtime_error("Cannot write: " + efficiency_output->string());
            efficiency_csv << std::setprecision(12)
                << "model,darcy_factor,inlet_reynolds,rpm,mass_flow_kg_s,"
                   "measured_eta_tt,predicted_eta_tt,error_percentage_points,"
                   "predicted_total_pressure_ratio,diffuser_loss_percent,"
                   "impeller_skin_friction_loss_j_kg,"
                   "impeller_blade_loading_loss_j_kg,"
                   "impeller_friction_factor,impeller_diffusion_factor,"
                   "impeller_incidence_loss_j_kg,inlet_relative_angle_from_meridional_deg,"
                   "optimum_inlet_angle_deg,estimated_throat_area_m2,"
                   "critical_throat_flow_kg_s,throat_capacity_utilization,"
                   "throat_capacity_margin_percent,choke_screen_flag\n";
            std::cout << "\nEfficiency digitization comparison (percentage points):\n"
                         "Darcy factor, complete, bias_pp, MAE_pp, max_abs_pp\n";
            for (const double f : friction_factors) {
                double bias = 0.0;
                double mae = 0.0;
                double max_abs = 0.0;
                for (const auto& point : efficiency_points) {
                    const MapPoint request{point.rpm, point.mass_flow, 1.0};
                    const auto result = predict(
                        request, air, derivatives, g, f, true, throat_capacity_cache);
                    const double error_pp = 100.0 *
                        (result.total_to_total_efficiency - point.measured_efficiency);
                    efficiency_csv << "fixed_galvas_incidence," << f << ',' << result.inlet_reynolds
                        << ',' << point.rpm << ',' << point.mass_flow << ','
                        << point.measured_efficiency << ','
                        << result.total_to_total_efficiency << ',' << error_pp << ','
                        << result.stage_pr << ',' << result.diffuser_loss_percent << ','
                        << result.impeller_skin_friction_loss_j_kg << ','
                        << result.impeller_blade_loading_loss_j_kg << ','
                        << result.impeller_skin_friction_factor << ','
                        << result.impeller_diffusion_factor << ','
                        << result.impeller_incidence_loss_j_kg << ','
                        << result.inlet_relative_angle_deg << ','
                        << result.optimum_inlet_angle_deg << ','
                        << result.estimated_throat_area_m2 << ','
                        << result.critical_throat_flow_kg_s << ','
                        << result.throat_capacity_utilization << ','
                        << result.throat_capacity_margin_percent << ','
                        << (result.choke_screen_flag ? 1 : 0) << '\n';
                    bias += error_pp;
                    mae += std::abs(error_pp);
                    max_abs = std::max(max_abs, std::abs(error_pp));
                }
                const double n = static_cast<double>(efficiency_points.size());
                std::cout << f << ", " << n << ", " << bias / n << ", "
                          << mae / n << ", " << max_abs << '\n';
            }
            double bias = 0.0;
            double mae = 0.0;
            double max_abs = 0.0;
            for (const auto& point : efficiency_points) {
                const MapPoint request{point.rpm, point.mass_flow, 1.0};
                const auto result = predict(
                    request, air, derivatives, g, -1.0, true, throat_capacity_cache);
                const double error_pp = 100.0 *
                    (result.total_to_total_efficiency - point.measured_efficiency);
                efficiency_csv << "smooth_wall_inlet_Re_galvas_incidence," << result.darcy_factor
                    << ',' << result.inlet_reynolds << ',' << point.rpm << ','
                    << point.mass_flow << ',' << point.measured_efficiency << ','
                    << result.total_to_total_efficiency << ',' << error_pp << ','
                    << result.stage_pr << ',' << result.diffuser_loss_percent << ','
                    << result.impeller_skin_friction_loss_j_kg << ','
                    << result.impeller_blade_loading_loss_j_kg << ','
                    << result.impeller_skin_friction_factor << ','
                    << result.impeller_diffusion_factor << ','
                    << result.impeller_incidence_loss_j_kg << ','
                    << result.inlet_relative_angle_deg << ','
                    << result.optimum_inlet_angle_deg << ','
                    << result.estimated_throat_area_m2 << ','
                    << result.critical_throat_flow_kg_s << ','
                    << result.throat_capacity_utilization << ','
                    << result.throat_capacity_margin_percent << ','
                    << (result.choke_screen_flag ? 1 : 0) << '\n';
                bias += error_pp;
                mae += std::abs(error_pp);
                max_abs = std::max(max_abs, std::abs(error_pp));
            }
            const double n = static_cast<double>(efficiency_points.size());
            std::cout << "smooth_wall_inlet_Re_galvas_incidence, " << n << ", " << bias / n
                      << ", " << mae / n << ", " << max_abs << '\n'
                      << "Efficiency comparison CSV: " << *efficiency_output << '\n';

            double no_incidence_bias = 0.0;
            double no_incidence_mae = 0.0;
            double no_incidence_max_abs = 0.0;
            for (const auto& point : efficiency_points) {
                const MapPoint request{point.rpm, point.mass_flow, 1.0};
                const auto result = predict(
                    request, air, derivatives, g, -1.0, false, throat_capacity_cache);
                const double error_pp = 100.0 *
                    (result.total_to_total_efficiency - point.measured_efficiency);
                efficiency_csv << "smooth_wall_inlet_Re_no_incidence,"
                    << result.darcy_factor << ',' << result.inlet_reynolds
                    << ',' << point.rpm << ',' << point.mass_flow << ','
                    << point.measured_efficiency << ','
                    << result.total_to_total_efficiency << ',' << error_pp << ','
                    << result.stage_pr << ',' << result.diffuser_loss_percent << ','
                    << result.impeller_skin_friction_loss_j_kg << ','
                    << result.impeller_blade_loading_loss_j_kg << ','
                    << result.impeller_skin_friction_factor << ','
                    << result.impeller_diffusion_factor << ','
                    << result.impeller_incidence_loss_j_kg << ','
                    << result.inlet_relative_angle_deg << ','
                    << result.optimum_inlet_angle_deg << ','
                    << result.estimated_throat_area_m2 << ','
                    << result.critical_throat_flow_kg_s << ','
                    << result.throat_capacity_utilization << ','
                    << result.throat_capacity_margin_percent << ','
                    << (result.choke_screen_flag ? 1 : 0) << '\n';
                no_incidence_bias += error_pp;
                no_incidence_mae += std::abs(error_pp);
                no_incidence_max_abs = std::max(
                    no_incidence_max_abs, std::abs(error_pp));
            }
            std::cout << "smooth_wall_inlet_Re_no_incidence, " << n << ", "
                      << no_incidence_bias / n << ", "
                      << no_incidence_mae / n << ", "
                      << no_incidence_max_abs << '\n';

            const auto ablation_output = efficiency_output->parent_path() /
                "loss_term_ablation.csv";
            std::ofstream ablation_csv(ablation_output);
            if (!ablation_csv)
                throw std::runtime_error(
                    "Cannot write: " + ablation_output.string());
            ablation_csv << std::setprecision(12)
                << "scenario,rpm,mass_flow_kg_s,measured_eta_tt,predicted_eta_tt,"
                   "eta_error_percentage_points,predicted_total_pressure_ratio,"
                   "predicted_impeller_pressure_ratio,diffuser_pressure_ratio,"
                   "darcy_factor,impeller_skin_friction_loss_j_kg,"
                   "impeller_blade_loading_loss_j_kg,impeller_incidence_loss_j_kg,"
                   "impeller_clearance_loss_j_kg,impeller_mixing_loss_j_kg,"
                   "impeller_disk_friction_loss_j_kg,impeller_hydraulic_diameter_m,"
                   "impeller_reynolds,impeller_darcy_friction_factor,"
                   "internal_loss_entropy_generation_j_kg_k,"
                   "impeller_specific_work_j_kg,outlet_total_temperature_k,"
                   "diffuser_loss_percent,impeller_blade_path_length_m,"
                   "impeller_mean_relative_velocity_m_s,impeller_mean_density_kg_m3,"
                   "impeller_mean_temperature_k,impeller_skin_fanning_factor,"
                   "impeller_mixing_chi,impeller_mixing_bstar,"
                   "impeller_correlation_exit_angle_deg,"
                   "impeller_inlet_mean_relative_velocity_m_s,"
                   "impeller_inlet_shroud_relative_velocity_m_s,"
                   "impeller_exit_relative_velocity_m_s,"
                   "impeller_throat_relative_velocity_m_s,"
                   "impeller_throat_density_kg_m3,impeller_throat_enthalpy_j_kg,"
                   "impeller_throat_area_m2,impeller_throat_blade_speed_m_s,"
                   "impeller_throat_rothalpy_j_kg,"
                   "impeller_throat_mass_residual_kg_s,"
                   "impeller_throat_rothalpy_residual_j_kg,"
                   "impeller_throat_entropy_residual_j_kg_k,"
                   "impeller_throat_mach,predicted_near_total_pressure_pa,"
                   "predicted_far_total_temperature_k,impeller_choke_loss_j_kg,"
                   "impeller_entrance_diffusion_loss_j_kg,"
                   "impeller_recirculation_loss_j_kg,"
                   "stage_total_enthalpy_rise_j_kg,impeller_euler_work_j_kg,"
                   "external_loss_enthalpy_addition_j_kg,"
                   "stage_energy_residual_vs_euler_j_kg,"
                   "stage_energy_residual_vs_euler_disk_recirc_j_kg,"
                   "stage_energy_residual_vs_euler_external_addition_j_kg\n";
            struct AblationScenario {
                std::string name;
                LossModelOptions losses;
                double diffuser_darcy_factor;
            };
            LossModelOptions full_losses;
            LossModelOptions published_2042 = full_losses;
            published_2042.published_set = PublishedLossSet::EckardtA2042;
            LossModelOptions published_2137 = full_losses;
            published_2137.published_set = PublishedLossSet::CrossGeometry2137;
            LossModelOptions shroud_proxy = full_losses;
            shroud_proxy.use_shroud_inlet_speed_proxy = true;
            LossModelOptions distributed_skin = full_losses;
            distributed_skin.distributed_skin_friction = true;
            LossModelOptions geometry_continuity_skin = full_losses;
            geometry_continuity_skin.geometry_continuity_skin_friction = true;
            LossModelOptions no_skin = full_losses;
            no_skin.skin_friction = false;
            LossModelOptions no_blade = full_losses;
            no_blade.blade_loading = false;
            LossModelOptions no_incidence = full_losses;
            no_incidence.incidence = false;
            LossModelOptions no_impeller = full_losses;
            no_impeller.skin_friction = false;
            no_impeller.blade_loading = false;
            no_impeller.incidence = false;
            no_impeller.clearance = false;
            LossModelOptions no_diffuser = full_losses;
            LossModelOptions mean_temperature_closure = full_losses;
            mean_temperature_closure.closure =
                LossClosure::MeanTemperatureEntropy;
            LossModelOptions no_clearance = full_losses;
            no_clearance.clearance = false;
            LossModelOptions no_mixing = full_losses;
            no_mixing.mixing = false;
            LossModelOptions no_disk = full_losses;
            no_disk.disk_friction = false;
            LossModelOptions no_mixing_disk = full_losses;
            no_mixing_disk.mixing = false;
            no_mixing_disk.disk_friction = false;
            LossModelOptions rectangular_diameter = full_losses;
            rectangular_diameter.hydraulic_diameter_model =
                HydraulicDiameterModel::RectangularEndpointApproximation;
            LossModelOptions blade_loading_shroud_consistent = full_losses;
            blade_loading_shroud_consistent
                .blade_loading_use_shroud_speed_consistently = true;
            LossModelOptions blade_loading_06 = blade_loading_shroud_consistent;
            blade_loading_06.blade_loading_loading_coefficient = 0.6;
            LossModelOptions incidence_mean_speed = full_losses;
            incidence_mean_speed.incidence_use_mean_inlet_speed = true;
            LossModelOptions clearance_exit_speed = full_losses;
            clearance_exit_speed.clearance_use_exit_meridional_speed = true;
            LossModelOptions disk_diameter_gap = full_losses;
            disk_diameter_gap.disk_friction_use_diameter_clearance_ratio = true;
            LossModelOptions source_consistent_losses = blade_loading_06;
            source_consistent_losses.incidence_use_mean_inlet_speed = true;
            source_consistent_losses.clearance_use_exit_meridional_speed = true;
            source_consistent_losses.disk_friction_use_diameter_clearance_ratio = true;
            std::vector<AblationScenario> ablations{
                {"all_losses_enthalpy_equivalent", full_losses, -1.0},
                {"published_eckardt_a_set_2042", published_2042, -1.0},
                {"published_cross_geometry_set_2137", published_2137, -1.0},
                {"distributed_skin_friction_screen", distributed_skin, -1.0},
                {"geometry_continuity_skin_friction_screen",
                 geometry_continuity_skin, -1.0},
                {"shroud_inlet_speed_proxy", shroud_proxy, -1.0},
                {"legacy_mean_temperature_entropy", mean_temperature_closure, -1.0},
                {"rectangular_endpoint_hydraulic_diameter", rectangular_diameter, -1.0},
                {"without_skin_friction", no_skin, -1.0},
                {"without_blade_loading", no_blade, -1.0},
                {"without_incidence", no_incidence, -1.0},
                {"without_clearance", no_clearance, -1.0},
                {"without_mixing", no_mixing, -1.0},
                {"without_disk_friction", no_disk, -1.0},
                {"without_mixing_and_disk_friction", no_mixing_disk, -1.0},
                {"without_impeller_losses", no_impeller, -1.0},
                {"without_diffuser_friction", no_diffuser, 0.0},
                {"audit_blade_loading_shroud_speed", blade_loading_shroud_consistent, -1.0},
                {"audit_blade_loading_06_shroud_speed", blade_loading_06, -1.0},
                {"audit_incidence_mean_speed", incidence_mean_speed, -1.0},
                {"audit_clearance_exit_meridional_speed", clearance_exit_speed, -1.0},
                {"audit_disk_diameter_gap_ratio", disk_diameter_gap, -1.0},
                {"audit_combined_alternative_forms", source_consistent_losses, -1.0},
            };

            // b* is diffuser-inlet depth / impeller-exit depth, not the
            // diffuser outlet taper ratio b3/b2=0.51. The geometric ratio is
            // approximately one here. Johnston-Dean documents its sudden-
            // expansion expression for b*>1, so b*=1 is a boundary/default
            // case and values just above one show the in-domain response.
            for (double chi : {0.05, 0.10, 0.20, 0.25, 0.30}) {
                auto losses = full_losses;
                losses.mixing_chi = chi;
                std::ostringstream name;
                name << "mixing_chi_" << std::fixed << std::setprecision(2)
                     << chi << "_bstar_1.00";
                ablations.push_back({name.str(), losses, -1.0});
            }
            for (double bstar : {1.02, 1.10}) {
                auto losses = full_losses;
                losses.mixing_bstar = bstar;
                std::ostringstream name;
                name << "mixing_chi_0.15_bstar_" << std::fixed
                     << std::setprecision(2) << bstar;
                ablations.push_back({name.str(), losses, -1.0});
            }
            {
                auto losses = full_losses;
                losses.correlation_exit_angle_deg = 60.0;
                ablations.push_back({"correlation_exit_angle_60deg_screen",
                                     losses, -1.0});
            }
            std::cout << "\nOne-at-a-time loss ablation and bounded sensitivity (efficiency points):\n"
                         "scenario,complete,bias_pp,MAE_pp,max_abs_pp,"
                         "mean_work_J_kg,mean_pressure_ratio\n";
            for (const auto& scenario : ablations) {
                double scenario_bias = 0.0;
                double scenario_mae = 0.0;
                double scenario_max_abs = 0.0;
                double scenario_work = 0.0;
                double scenario_pr = 0.0;
                for (const auto& point : efficiency_points) {
                    const MapPoint request{point.rpm, point.mass_flow, 1.0};
                    const auto result = [&] {
                        try {
                            return predict(
                                request, air, derivatives, g,
                                scenario.diffuser_darcy_factor,
                                scenario.losses.incidence,
                                throat_capacity_cache,
                                scenario.losses.skin_friction,
                                scenario.losses.blade_loading,
                                scenario.losses);
                        }
                        catch (const std::exception& error) {
                            throw std::runtime_error(
                                "Ablation scenario '" + scenario.name +
                                "' failed at " + std::to_string(point.rpm) +
                                " rpm, " + std::to_string(point.mass_flow) +
                                " kg/s: " + error.what());
                        }
                    }();
                    const double error_pp = 100.0 *
                        (result.total_to_total_efficiency -
                         point.measured_efficiency);
                    ablation_csv << scenario.name << ',' << point.rpm << ','
                        << point.mass_flow << ',' << point.measured_efficiency
                        << ',' << result.total_to_total_efficiency << ','
                        << error_pp << ',' << result.stage_pr << ','
                        << result.impeller_pr << ',' << result.diffuser_pr << ','
                        << result.darcy_factor << ','
                        << result.impeller_skin_friction_loss_j_kg << ','
                        << result.impeller_blade_loading_loss_j_kg << ','
                        << result.impeller_incidence_loss_j_kg << ','
                        << result.impeller_clearance_loss_j_kg << ','
                        << result.impeller_mixing_loss_j_kg << ','
                        << result.impeller_disk_friction_loss_j_kg << ','
                        << result.impeller_hydraulic_diameter_m << ','
                        << result.impeller_reynolds << ','
                        << result.impeller_skin_friction_factor << ','
                        << result.internal_loss_entropy_generation_j_kg_k << ','
                        << result.impeller_specific_work_j_kg << ','
                        << result.outlet_total_temperature_k << ','
                        << result.diffuser_loss_percent << ','
                        << result.impeller_blade_path_length_m << ','
                        << result.impeller_mean_relative_velocity_m_s << ','
                        << result.impeller_mean_density_kg_m3 << ','
                        << result.impeller_mean_temperature_k << ','
                        << result.impeller_skin_fanning_factor << ','
                        << result.impeller_mixing_chi << ','
                        << result.impeller_mixing_bstar << ','
                        << result.impeller_correlation_exit_angle_deg << ','
                        << result.impeller_inlet_mean_relative_velocity_m_s << ','
                        << result.impeller_inlet_shroud_relative_velocity_m_s << ','
                        << result.impeller_exit_relative_velocity_m_s << ','
                        << result.impeller_throat_relative_velocity_m_s << ','
                        << result.impeller_throat_density_kg_m3 << ','
                        << result.impeller_throat_enthalpy_j_kg << ','
                        << result.impeller_throat_area_m2 << ','
                        << result.impeller_throat_blade_speed_m_s << ','
                        << result.impeller_throat_rothalpy_j_kg << ','
                        << result.impeller_throat_mass_residual_kg_s << ','
                        << result.impeller_throat_rothalpy_residual_j_kg << ','
                        << result.impeller_throat_entropy_residual_j_kg_k << ','
                        << result.impeller_throat_mach << ','
                        << result.near_traverse.total_pressure_pa << ','
                        << result.far_traverse.total_temperature_k << ','
                        << result.impeller_choke_loss_j_kg << ','
                        << result.impeller_entrance_diffusion_loss_j_kg << ','
                        << result.impeller_recirculation_loss_j_kg << ','
                        << result.stage_total_enthalpy_rise_j_kg << ','
                        << result.impeller_euler_work_j_kg << ','
                        << result.external_loss_enthalpy_addition_j_kg << ','
                        << result.stage_energy_residual_vs_euler_j_kg << ','
                        << result.stage_energy_residual_vs_euler_disk_recirc_j_kg << ','
                        << result.stage_energy_residual_vs_euler_external_addition_j_kg
                        << '\n';
                    scenario_bias += error_pp;
                    scenario_mae += std::abs(error_pp);
                    scenario_max_abs = std::max(
                        scenario_max_abs, std::abs(error_pp));
                    scenario_work += result.impeller_specific_work_j_kg;
                    scenario_pr += result.stage_pr;
                }
                const double n_points =
                    static_cast<double>(efficiency_points.size());
                std::cout << scenario.name << ',' << n_points << ','
                    << scenario_bias / n_points << ','
                    << scenario_mae / n_points << ','
                    << scenario_max_abs << ','
                    << scenario_work / n_points << ','
                    << scenario_pr / n_points << '\n';
            }
            std::cout << "Loss-term ablation CSV: "
                      << ablation_output << '\n';
        }
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "Eckardt diffuser sweep failed: " << error.what() << '\n';
        return 1;
    }
}
#endif
