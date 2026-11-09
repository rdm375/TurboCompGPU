#include <turbocompgpu/complete_stage_solver.hpp>
#include <turbocompgpu/station.hpp>

#include <thermomodel/caloric.hpp>
#include <thermomodel/peng_robinson.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace turbocompgpu;

namespace {

struct Point {
    double T0, p0, rpm, mass, measured_pr, measured_eta;
};

Point read_point(const std::filesystem::path& path)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open " + path.string());
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line.front() == '#') continue;
        std::stringstream row(line);
        Point p{};
        char c1, c2, c3, c4, c5;
        if (!(row >> p.T0 >> c1 >> p.p0 >> c2 >> p.rpm >> c3 >> p.mass >> c4
                  >> p.measured_pr >> c5 >> p.measured_eta) ||
            c1 != ',' || c2 != ',' || c3 != ',' || c4 != ',' || c5 != ',')
            throw std::runtime_error("Invalid design-point CSV row: " + line);
        return p;
    }
    throw std::runtime_error("No design-point row in " + path.string());
}

thermomodel::PengRobinson make_r134a()
{
    thermomodel::PRParameters p;
    p.molar_mass = 0.102032;
    p.Tc = 374.21;
    p.pc = 4.0593e6;
    p.acentric = 0.32684;
    p.Tref = 300.0;
    p.pref = 101325.0;
    p.Tmin = 170.0;
    p.Tmax = 700.0;
    p.pmin = 1.0e3;
    p.pmax = 20.0e6;
    p.root_policy = thermomodel::PRRootPolicy::Vapor;

    // Deliberately simple ideal-gas caloric term for a runnable screen.
    // This is not a reference-quality R134a property model.
    return thermomodel::PengRobinson(
        std::make_shared<thermomodel::ConstantCp>(850.0), p);
}

ImpellerSystemGeometry schiffmann_geometry()
{
    constexpr double hub_r = 0.0020;
    constexpr double shroud_r = 0.0056;
    constexpr double exit_r = 0.0100;
    constexpr double blade_t = 0.0001;
    constexpr double full_blades = 9.0;
    constexpr double splitters = 9.0;
    const double mean_inlet_r = std::sqrt(0.5 * (hub_r*hub_r + shroud_r*shroud_r));
    const double inlet_blockage = full_blades * blade_t /
        (2.0 * std::numbers::pi * mean_inlet_r);
    const double exit_blockage = (full_blades + splitters) * blade_t /
        (2.0 * std::numbers::pi * exit_r);
    return {
        {2.0*hub_r, 2.0*shroud_r, inlet_blockage},
        {2.0*exit_r, 0.0012, exit_blockage},
        {2.0*exit_r, 45.0, 2.0*mean_inlet_r}
    };
}

thermomodel::State inlet_static_from_total(
    const thermomodel::Model& eos, const Point& point, double area)
{
    const auto total_result = eos.evaluate_Tp(point.T0, point.p0);
    if (!total_result) throw std::runtime_error("R134a total-state EOS failed: " + total_result.message);
    const auto total = total_result.value;

    auto state_at_entropy = [&](double temperature) {
        double lo = 1.0e3, hi = point.p0;
        const auto low = eos.evaluate_Tp(temperature, lo);
        const auto high = eos.evaluate_Tp(temperature, hi);
        if (!low || !high || low.value.s < total.s || high.value.s > total.s)
            throw std::runtime_error("Could not bracket R134a static inlet pressure");
        for (int i = 0; i < 48; ++i) {
            const double mid = 0.5 * (lo + hi);
            const auto candidate = eos.evaluate_Tp(temperature, mid);
            if (!candidate) throw std::runtime_error("R134a EOS failed in pressure inversion");
            if (candidate.value.s > total.s) lo = mid;
            else hi = mid;
        }
        const auto result = eos.evaluate_Tp(temperature, 0.5 * (lo + hi));
        if (!result) throw std::runtime_error("R134a static inlet state failed");
        return result.value;
    };

    auto state_for_velocity = [&](double velocity) {
        const double target_h = total.h - 0.5 * velocity * velocity;
        double hi = point.T0;
        double interval = 1.0e-3;
        double lo = hi - interval;
        auto low = state_at_entropy(lo);
        while (low.h > target_h && interval < 40.0) {
            interval = std::min(40.0, 2.0*interval);
            lo = point.T0 - interval;
            low = state_at_entropy(lo);
        }
        if (low.h > target_h) throw std::runtime_error("Could not bracket R134a static inlet temperature");
        for (int i = 0; i < 44; ++i) {
            const double mid = 0.5 * (lo + hi);
            const auto candidate = state_at_entropy(mid);
            if (candidate.h < target_h) lo = mid;
            else hi = mid;
        }
        return state_at_entropy(0.5 * (lo + hi));
    };

    double lo = 0.0, hi = 150.0;
    for (int i = 0; i < 36; ++i) {
        const double velocity = 0.5 * (lo + hi);
        const auto state = state_for_velocity(velocity);
        if (state.rho * area * velocity < point.mass) lo = velocity;
        else hi = velocity;
    }
    const double velocity = 0.5 * (lo + hi);
    return state_for_velocity(velocity);
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const std::filesystem::path input = argc > 1 ? argv[1] :
            "data/validation/schiffmann_favrat_design_point.csv";
        const std::filesystem::path output = argc > 2 ? argv[2] :
            "schiffmann_favrat_design_point_screen.csv";
        const auto point = read_point(input);
        const auto eos = make_r134a();
        const auto geometry = schiffmann_geometry();
        const double inlet_area = effective_flow_area(
            inlet_annulus_area(geometry.inlet), geometry.inlet.blockage);
        const auto inlet_static = inlet_static_from_total(eos, point, inlet_area);

        // Effective slip is a provisional correlation input; no fit is made.
        constexpr double slip = 1.0 -
            std::sqrt(std::cos(45.0 * std::numbers::pi / 180.0)) /
            std::pow(9.0, 0.7);
        constexpr double pressure_scale = 500000.0;
        const OperatingPointRequest request{
            {inlet_static.T, inlet_static.p, point.mass, point.rpm, slip, 0.0, 0.0},
            {inlet_static.T + 45.0,
             std::log(2.0 * point.p0 / pressure_scale)}
        };
        NonlinearResidualScales scales{30000.0, 100.0, pressure_scale};
        NonlinearSolverOptions solver_options;
        solver_options.residual_tolerance = 1.0e-10;
        const auto stage = evaluate_complete_stage(
            eos, eos, geometry,
            {geometry.exit.diameter/2.0,
             1.65*geometry.exit.diameter/2.0,
             0.75*geometry.exit.width, 0.0},
            request, scales, solver_options);
        if (!stage.complete())
            throw std::runtime_error("R134a stage screen failed: " + stage.message);

        const auto& station3 = *stage.stage->station3;
        const FlowStation outlet_flow{
            station3.state,
            {station3.meridional_velocity, station3.tangential_velocity}
        };
        const auto outlet_total = evaluate_stagnation_state(eos, outlet_flow);
        const auto inlet_total = eos.evaluate_Tp(point.T0, point.p0);
        if (!inlet_total) throw std::runtime_error("R134a inlet total state failed");
        const double predicted_pr = outlet_total.stagnation_state.p / point.p0;
        const double predicted_eta =
            (invert_temperature(eos, outlet_total.stagnation_state.p,
                inlet_total.value.s, Property::Entropy, point.T0 + 30.0).state.h -
             inlet_total.value.h) /
            (outlet_total.stagnation_enthalpy - inlet_total.value.h);
        const double pr_error_pct = 100.0*(predicted_pr/point.measured_pr - 1.0);
        const double eta_error_pp = 100.0*(predicted_eta - point.measured_eta);

        if (!output.parent_path().empty())
            std::filesystem::create_directories(output.parent_path());
        std::ofstream csv(output);
        if (!csv) throw std::runtime_error("Cannot write " + output.string());
        csv << std::setprecision(12)
            << "rpm,mass_flow_kg_s,inlet_total_T_K,inlet_total_p_Pa,"
               "inlet_static_T_K,inlet_static_p_Pa,slip_factor,"
               "measured_PRtt,predicted_PRtt,PR_error_percent,"
               "measured_eta_tt,predicted_eta_tt,eta_error_percentage_points\n"
            << point.rpm << ',' << point.mass << ',' << point.T0 << ',' << point.p0 << ','
            << inlet_static.T << ',' << inlet_static.p << ',' << slip << ','
            << point.measured_pr << ',' << predicted_pr << ',' << pr_error_pct << ','
            << point.measured_eta << ',' << predicted_eta << ',' << eta_error_pp << '\n';
        std::cout << std::fixed << std::setprecision(4)
                  << "Schiffmann-Favrat R134a summary point: PR measured="
                  << point.measured_pr << ", model=" << predicted_pr
                  << " (error " << pr_error_pct << "%), eta measured="
                  << point.measured_eta << ", model=" << predicted_eta
                  << " (error " << eta_error_pp << " pp)\nCSV: " << output << '\n';
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "Schiffmann-Favrat screen failed: " << error.what() << '\n';
        return 1;
    }
}
