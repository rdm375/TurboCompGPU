#include <turbocompgpu/complete_stage_solver.hpp>
#include <turbocompgpu/station.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <string>

using namespace turbocompgpu;

namespace {

constexpr std::array<double, 11> flows{
    2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0
};

constexpr std::array<double, 21> speeds{
    9000.0, 9300.0, 9600.0, 9900.0, 10200.0,
    10500.0, 10800.0, 11100.0, 11400.0, 11700.0,
    12000.0, 12300.0, 12600.0, 12900.0, 13200.0,
    13500.0, 13800.0, 14100.0, 14400.0, 14700.0,
    15000.0
};

constexpr int boundary_grid_intervals = 20;

ImpellerSystemGeometry geometry()
{
    return {
        {0.10, 0.20, 0.05},
        {0.40, 0.025, 0.10},
        {0.40, 45.0, 0.15}
    };
}

NonlinearResidualScales scales()
{
    return {50000.0, 100.0, 500000.0};
}

OperatingPointRequest request(double flow, double rpm)
{
    return {
        {300.0, 500000.0, flow, rpm, 0.90, 20.0, 20.0},
        {340.0, std::log(700000.0 / 500000.0)}
    };
}

NonlinearSolverOptions options()
{
    NonlinearSolverOptions result;
    result.residual_tolerance = 1.0e-13;
    return result;
}

const char* stage_status(CompleteStageStatus status)
{
    switch (status) {
    case CompleteStageStatus::Complete:
        return "Complete";
    case CompleteStageStatus::ImpellerFailed:
        return "ImpellerFailed";
    case CompleteStageStatus::DiffuserFailed:
        return "DiffuserFailed";
    }
    return "Unknown";
}

const char* diffuser_status(DiffuserBranchStatus status)
{
    switch (status) {
    case DiffuserBranchStatus::Converged:
        return "Converged";
    case DiffuserBranchStatus::Critical:
        return "Critical";
    case DiffuserBranchStatus::AboveCriticalCapacity:
        return "AboveCriticalCapacity";
    case DiffuserBranchStatus::BranchUnavailable:
        return "BranchUnavailable";
    case DiffuserBranchStatus::ThermodynamicDomainLimit:
        return "ThermodynamicDomainLimit";
    }
    return "Unknown";
}

template<class Model>
VanelessDiffuserGeometry calibrate(
    const Model& model,
    const ImpellerSystemGeometry& g,
    const NonlinearResidualScales& s)
{
    const auto op = request(5.0, 12000.0);

    const auto impeller = evaluate_operating_point(
        model, g, op, s, options());

    if (!(impeller.converged &&
          impeller.independently_validated &&
          impeller.basic_physics_valid)) {
        throw std::runtime_error(
            "Reference impeller calibration failed");
    }

    const auto state2 = evaluate(
        model,
        impeller.outlet_temperature,
        impeller.outlet_pressure);

    const double mass = op.operating.mass_flow;
    const double r2 = g.exit.diameter / 2.0;
    const double r3 = 1.5 * r2;

    const double area2 = effective_flow_area(
        impeller_exit_area(g.exit), g.exit.blockage);

    const double vm2 = mass / (state2.rho * area2);

    const double omega =
        2.0 * std::numbers::pi *
        op.operating.rotational_speed_rpm / 60.0;

    const double beta =
        g.impeller.exit_blade_angle_deg *
        std::numbers::pi / 180.0;

    const double whirl2 =
        op.operating.slip_factor *
        (omega * r2 - vm2 / std::tan(beta));

    const double whirl3 = r2 * whirl2 / r3;

    const FlowStation station2{
        state2, {vm2, whirl2}
    };

    const auto stagnation2 =
        evaluate_stagnation_state(model, station2);

    const double p_reference =
        state2.p + 0.10 *
        (stagnation2.stagnation_state.p - state2.p);

    const auto reference = invert_temperature(
        model,
        p_reference,
        state2.s,
        Property::Entropy,
        state2.T);

    const double vm3_squared =
        2.0 * (
            stagnation2.stagnation_enthalpy -
            reference.state.h
        ) - whirl3 * whirl3;

    if (!(vm3_squared > 0.0)) {
        throw std::runtime_error(
            "Invalid diffuser calibration velocity");
    }

    const double area3 =
        mass / (reference.state.rho * std::sqrt(vm3_squared));

    const double width3 =
        area3 / (2.0 * std::numbers::pi * r3);

    return {r2, r3, width3, 0.0};
}


struct BoundaryPoint {
    double mass = 0.0;
    double residual = 0.0;
    bool valid = false;
};

struct BoundaryResult {
    bool bracketed = false;
    bool converged = false;
    double mass = std::numeric_limits<double>::quiet_NaN();
    double residual = std::numeric_limits<double>::quiet_NaN();
    double bracket_width = std::numeric_limits<double>::quiet_NaN();
    int iterations = 0;
};

enum class BoundaryKind {
    Recovery,
    Capacity
};

template<class Model>
BoundaryPoint evaluate_boundary_point(
    const Model& model,
    const ImpellerSystemGeometry& g,
    const VanelessDiffuserGeometry& diffuser,
    const NonlinearResidualScales& s,
    double mass,
    double rpm,
    BoundaryKind kind)
{
    const auto result = evaluate_complete_stage(
        model, model, g, diffuser,
        request(mass, rpm), s, options());

    if (!result.impeller_valid())
        return {mass, 0.0, false};

    if (kind == BoundaryKind::Capacity) {
        if (!result.stage)
            return {mass, 0.0, false};

        const double critical =
            result.stage->diffuser_critical_mass_flow;

        if (!(std::isfinite(critical) && critical > 0.0))
            return {mass, 0.0, false};

        return {
            mass,
            critical - mass,
            true
        };
    }

    if (!result.complete() || !result.station3_valid())
        return {mass, 0.0, false};

    return {
        mass,
        result.stage->station3->pressure -
            result.stage->station2_static.p,
        true
    };
}

template<class Model>
BoundaryResult find_boundary(
    const Model& model,
    const ImpellerSystemGeometry& g,
    const VanelessDiffuserGeometry& diffuser,
    const NonlinearResidualScales& s,
    double rpm,
    BoundaryKind kind)
{
    BoundaryResult result;

    // First locate a sign-changing interval.
    //
    // A failed intermediate evaluation breaks continuity
    // of the numerical evidence: we do not bracket across
    // an invalid point.

    BoundaryPoint previous;
    bool previous_valid = false;

    for (double mass : flows) {
        const auto current = evaluate_boundary_point(
            model, g, diffuser, s, mass, rpm, kind);

        if (!current.valid) {
            previous_valid = false;
            continue;
        }

        if (current.residual == 0.0) {
            result.bracketed = true;
            result.converged = true;
            result.mass = mass;
            result.residual = 0.0;
            result.bracket_width = 0.0;
            return result;
        }

        if (previous_valid &&
            std::signbit(previous.residual) !=
            std::signbit(current.residual)) {

            double lo = previous.mass;
            double hi = current.mass;
            double flo = previous.residual;

            result.bracketed = true;

            constexpr int max_iterations = 100;
            constexpr double relative_mass_tolerance = 1.0e-10;

            for (int iteration = 1;
                 iteration <= max_iterations;
                 ++iteration) {

                const double mid = lo + 0.5 * (hi - lo);

                if (mid == lo || mid == hi)
                    break;

                const auto point = evaluate_boundary_point(
                    model, g, diffuser, s, mid, rpm, kind);

                if (!point.valid) {
                    // Do not claim convergence if the solver
                    // cannot evaluate the bracket interior.
                    result.iterations = iteration;
                    result.bracket_width = hi - lo;
                    return result;
                }

                result.iterations = iteration;

                if (point.residual == 0.0) {
                    result.converged = true;
                    result.mass = mid;
                    result.residual = 0.0;
                    result.bracket_width = 0.0;
                    return result;
                }

                if (std::signbit(point.residual) ==
                    std::signbit(flo)) {
                    lo = mid;
                    flo = point.residual;
                } else {
                    hi = mid;
                }

                if ((hi - lo) <=
                    relative_mass_tolerance *
                    std::max(1.0, std::abs(mid))) {

                    const double final_mass =
                        lo + 0.5 * (hi - lo);

                    const auto final_point =
                        evaluate_boundary_point(
                            model, g, diffuser, s,
                            final_mass, rpm, kind);

                    if (!final_point.valid)
                        return result;

                    result.converged = true;
                    result.mass = final_mass;
                    result.residual = final_point.residual;
                    result.bracket_width = hi - lo;
                    return result;
                }
            }

            return result;
        }

        previous = current;
        previous_valid = true;
    }

    return result;
}

template<class Model>
void analyze(
    const Model& model,
    const std::string& name,
    std::ofstream& csv,
    std::ofstream& samples,
    const VanelessDiffuserGeometry& diffuser)
{
    const auto g = geometry();
    const auto s = scales();

    // Exactly the same reference calibration as 8B.5b.1.
    // Use the shared, externally calibrated geometry.

    std::cout << "\n" << name
              << " fixed diffuser width = "
              << diffuser.outlet_width << " m\n";

    for (double rpm : speeds) {
        const auto recovery = find_boundary(
            model, g, diffuser, s, rpm,
            BoundaryKind::Recovery);

        const auto capacity = find_boundary(
            model, g, diffuser, s, rpm,
            BoundaryKind::Capacity);

        const double nan =
            std::numeric_limits<double>::quiet_NaN();

        const double ratio =
            recovery.converged && capacity.converged
                ? recovery.mass / capacity.mass
                : nan;

        csv
            << name << ','
            << rpm << ','
            << diffuser.outlet_width << ','
            << recovery.bracketed << ','
            << recovery.converged << ','
            << recovery.mass << ','
            << recovery.residual << ','
            << recovery.bracket_width << ','
            << recovery.iterations << ','
            << capacity.bracketed << ','
            << capacity.converged << ','
            << capacity.mass << ','
            << capacity.residual << ','
            << capacity.bracket_width << ','
            << capacity.iterations << ','
            << ratio << '\n';

        if (recovery.converged && capacity.converged &&
            capacity.mass > recovery.mass) {
            for (int i = 0; i <= boundary_grid_intervals; ++i) {
                const double fraction =
                    static_cast<double>(i) /
                    static_cast<double>(boundary_grid_intervals);
                const double mass = recovery.mass + fraction *
                    (capacity.mass - recovery.mass);
                const auto op = request(mass, rpm);
                const auto stage = evaluate_complete_stage(
                    model, model, g, diffuser,
                    op, s, options(), DiffuserBranch::Subsonic);

                if (!stage.complete() || !stage.station3_valid()) {
                    throw std::runtime_error(
                        "Boundary-grid sample failed complete-stage validation");
                }

                const auto inlet_state = evaluate(
                    model,
                    op.operating.inlet_temperature,
                    op.operating.inlet_pressure);
                const double inlet_area = effective_flow_area(
                    inlet_annulus_area(g.inlet), g.inlet.blockage);
                const double inlet_vm =
                    mass / (inlet_state.rho * inlet_area);
                const FlowStation inlet_station{
                    inlet_state,
                    {inlet_vm, op.operating.inlet_prewhirl}
                };
                const auto inlet_total =
                    evaluate_stagnation_state(model, inlet_station);

                const auto& outlet = *stage.stage->station3;
                const FlowStation outlet_station{
                    outlet.state,
                    {outlet.meridional_velocity,
                     outlet.tangential_velocity}
                };
                const auto outlet_total =
                    evaluate_stagnation_state(model, outlet_station);

                const char* boundary_kind =
                    i == 0 ? "zero_recovery" :
                    i == boundary_grid_intervals
                        ? "critical_diffuser_capacity" : "interior";

                samples
                    << name << ','
                    << rpm << ','
                    << mass << ','
                    << fraction << ','
                    << recovery.mass << ','
                    << capacity.mass << ','
                    << boundary_kind << ','
                    << stage_status(stage.status) << ','
                    << diffuser_status(stage.stage->status) << ','
                    << stage.impeller.converged << ','
                    << stage.impeller.independently_validated << ','
                    << stage.impeller.basic_physics_valid << ','
                    << outlet_total.stagnation_state.T << ','
                    << outlet_total.stagnation_state.p << ','
                    << outlet_total.stagnation_state.p /
                       inlet_total.stagnation_state.p << ','
                    << stage.stage->impeller_specific_work << ','
                    << stage.stage->impeller_shaft_power << ','
                    << stage.stage->combined_energy_residual << ','
                    << stage.stage->diffuser.mass_residual << ','
                    << stage.stage->diffuser.entropy_residual << ','
                    << stage.stage->diffuser.angular_momentum_residual
                    << '\n';
            }
        }

        std::cout
            << "N=" << std::setw(6) << rpm
            << " recovery=";

        if (recovery.converged)
            std::cout << recovery.mass;
        else
            std::cout << "UNRESOLVED";

        std::cout << " capacity=";

        if (capacity.converged)
            std::cout << capacity.mass;
        else
            std::cout << "UNRESOLVED";

        std::cout << " ratio=" << ratio << '\n';
    }
}

} // namespace

int main()
{
    try {
        const std::filesystem::path output =
            "benchmarks/results-complete-stage-characteristic/"
            "boundaries-common-geometry.csv";
        const std::filesystem::path sample_output =
            "benchmarks/results-complete-stage-characteristic/"
            "complete-stage-boundary-grid.csv";

        std::filesystem::create_directories(
            output.parent_path());

        std::ofstream csv(output);
        std::ofstream samples(sample_output);

        if (!csv || !samples)
            throw std::runtime_error("Cannot open output CSV");

        csv << std::setprecision(17);
        samples << std::setprecision(17);

        csv
            << "model,rpm,diffuser_width_m,"
            << "recovery_bracketed,recovery_converged,"
            << "recovery_mass_kg_s,recovery_residual_pa,"
            << "recovery_bracket_width_kg_s,recovery_iterations,"
            << "capacity_bracketed,capacity_converged,"
            << "capacity_mass_kg_s,capacity_residual_kg_s,"
            << "capacity_bracket_width_kg_s,capacity_iterations,"
            << "recovery_to_capacity_ratio\n";

        samples
            << "model,rotational_speed_rpm,mass_flow,flow_fraction,"
            << "zero_recovery_flow_kg_s,critical_diffuser_capacity_kg_s,"
            << "boundary_kind,stage_status,diffuser_status,"
            << "converged,validated,basic_physics_valid,"
            << "outlet_total_temperature,outlet_total_pressure,"
            << "total_pressure_ratio,specific_work,shaft_power,"
            << "combined_energy_residual,diffuser_mass_residual,"
            << "diffuser_entropy_residual,diffuser_angular_momentum_residual\n";

        const thermomodel::IdealGas ideal;
        const auto common_diffuser =
            calibrate(ideal, geometry(), scales());

        analyze(ideal, "IdealGas", csv, samples, common_diffuser);

        const auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        thermomodel::PRParameters parameters;
        parameters.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        const thermomodel::PengRobinson pr(
            caloric, parameters);

        analyze(pr, "PengRobinson", csv, samples, common_diffuser);

        if (!csv || !samples)
            throw std::runtime_error("CSV write failed");

        std::cout << "\nCSV: " << output << '\n';
        std::cout << "Boundary-grid samples: "
                  << sample_output << '\n';

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
