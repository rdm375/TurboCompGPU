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

constexpr std::array<double, 7> flows{
    2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0
};

constexpr std::array<double, 5> speeds{
    9000.0, 10500.0, 12000.0, 13500.0, 15000.0
};

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
    std::ofstream& csv)
{
    const auto g = geometry();
    const auto s = scales();

    // Exactly the same reference calibration as 8B.5b.1.
    const auto diffuser = calibrate(model, g, s);

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
            "boundaries.csv";

        std::filesystem::create_directories(
            output.parent_path());

        std::ofstream csv(output);

        if (!csv)
            throw std::runtime_error("Cannot open output CSV");

        csv << std::setprecision(17);

        csv
            << "model,rpm,diffuser_width_m,"
            << "recovery_bracketed,recovery_converged,"
            << "recovery_mass_kg_s,recovery_residual_pa,"
            << "recovery_bracket_width_kg_s,recovery_iterations,"
            << "capacity_bracketed,capacity_converged,"
            << "capacity_mass_kg_s,capacity_residual_kg_s,"
            << "capacity_bracket_width_kg_s,capacity_iterations,"
            << "recovery_to_capacity_ratio\n";

        const thermomodel::IdealGas ideal;

        analyze(ideal, "IdealGas", csv);

        const auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        thermomodel::PRParameters parameters;
        parameters.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        const thermomodel::PengRobinson pr(
            caloric, parameters);

        analyze(pr, "PengRobinson", csv);

        if (!csv)
            throw std::runtime_error("CSV write failed");

        std::cout << "\nCSV: " << output << '\n';

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
