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

    // Extended experimental search interval.
    // Invalid evaluations break the bracket chain.
    constexpr double minimum_mass = 0.25;
    constexpr double maximum_mass = 25.0;
    constexpr double scan_increment = 0.05;

    for (int index = 0;
         minimum_mass + index * scan_increment
             <= maximum_mass + 1.0e-12;
         ++index) {

        const double mass =
            minimum_mass + index * scan_increment;

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


struct CapacityPoint {
    double mass = 0.0;
    double residual = 0.0;
    double critical_mass = 0.0;
    double critical_pressure = 0.0;
    double sonic_mach = 0.0;
    bool valid = false;
};

struct CapacityBoundary {
    bool bracketed = false;
    bool converged = false;
    bool search_limited = false;
    double mass = std::numeric_limits<double>::quiet_NaN();
    double residual = std::numeric_limits<double>::quiet_NaN();
    double critical_pressure =
        std::numeric_limits<double>::quiet_NaN();
    double sonic_mach = std::numeric_limits<double>::quiet_NaN();
    double bracket_width =
        std::numeric_limits<double>::quiet_NaN();
    int iterations = 0;
    int evaluations = 0;
    std::string reason = "NotBracketed";
};

template<class Model>
CapacityPoint evaluate_capacity(
    const Model& model,
    const ImpellerSystemGeometry& g,
    const VanelessDiffuserGeometry& diffuser,
    const NonlinearResidualScales& s,
    double mass,
    double rpm)
{
    CapacityPoint point;
    point.mass = mass;

    const auto result = evaluate_complete_stage(
        model, model, g, diffuser,
        request(mass, rpm), s, options());

    if (!result.impeller_valid() || !result.stage)
        return point;

    const auto& branch = result.stage->diffuser;

    if (!(std::isfinite(branch.critical_mass_flow) &&
          branch.critical_mass_flow > 0.0 &&
          std::isfinite(branch.critical.pressure) &&
          branch.critical.pressure > 0.0 &&
          std::isfinite(branch.critical.meridional_mach))) {
        return point;
    }

    point.critical_mass = branch.critical_mass_flow;
    point.critical_pressure = branch.critical.pressure;
    point.sonic_mach = branch.critical.meridional_mach;
    point.residual = point.critical_mass - mass;

    // The critical state must genuinely be sonic.
    if (std::abs(point.sonic_mach - 1.0) > 1.0e-8)
        return point;

    point.valid = true;
    return point;
}

template<class Model>
CapacityBoundary find_capacity(
    const Model& model,
    const ImpellerSystemGeometry& g,
    const VanelessDiffuserGeometry& diffuser,
    const NonlinearResidualScales& s,
    double rpm)
{
    CapacityBoundary result;

    auto evaluate_checked = [&](double mass) {
        ++result.evaluations;
        return evaluate_capacity(
            model, g, diffuser, s, mass, rpm);
    };

    // Begin at the established reference-grid lower bound.
    auto previous = evaluate_checked(0.25);

    if (!previous.valid) {
        result.reason = "InvalidInitialPoint";
        return result;
    }

    if (previous.residual <= 0.0) {
        result.reason = "AlreadyAboveCapacityAtLowerBound";
        return result;
    }

    // Continue outward in mass flow.
    //
    // The 25 kg/s ceiling is an experimental search limit,
    // not a physical compressor operating limit.
    constexpr double maximum_mass = 25.0;
    constexpr double increment = 0.25;

    CapacityPoint lower = previous;
    CapacityPoint upper;

    for (double mass = 0.25 + increment;
         mass <= maximum_mass + 1.0e-12;
         mass += increment) {

        const auto current = evaluate_checked(mass);

        if (!current.valid) {
            result.search_limited = true;
            result.reason = "ImpellerOrDiffuserEvaluationFailed";
            return result;
        }

        if (current.residual <= 0.0) {
            lower = previous;
            upper = current;
            result.bracketed = true;
            break;
        }

        previous = current;
    }

    if (!result.bracketed) {
        result.search_limited = true;
        result.reason = "SearchCeilingReached";
        return result;
    }

    double lo = lower.mass;
    double hi = upper.mass;
    double flo = lower.residual;

    constexpr int maximum_iterations = 100;
    constexpr double relative_mass_tolerance = 1.0e-10;
    constexpr double relative_capacity_tolerance = 1.0e-9;

    for (int iteration = 1;
         iteration <= maximum_iterations;
         ++iteration) {

        const double mid = lo + 0.5 * (hi - lo);

        if (mid == lo || mid == hi) {
            result.reason = "FloatingPointResolutionLimit";
            return result;
        }

        const auto point = evaluate_checked(mid);
        result.iterations = iteration;

        if (!point.valid) {
            result.reason = "InvalidBisectionPoint";
            return result;
        }

        if (point.residual > 0.0) {
            lo = mid;
            flo = point.residual;
        } else {
            hi = mid;
        }

        const double width = hi - lo;

        if (width <=
            relative_mass_tolerance *
            std::max(1.0, std::abs(mid))) {

            const double root = lo + 0.5 * width;
            const auto final_point = evaluate_checked(root);

            if (!final_point.valid) {
                result.reason = "InvalidFinalPoint";
                return result;
            }

            result.mass = root;
            result.residual = final_point.residual;
            result.critical_pressure =
                final_point.critical_pressure;
            result.sonic_mach = final_point.sonic_mach;
            result.bracket_width = width;

            const double residual_scale =
                std::max(1.0, root);

            result.converged =
                std::abs(result.residual) <=
                relative_capacity_tolerance * residual_scale;

            result.reason = result.converged
                ? "Converged"
                : "ResidualToleranceNotMet";

            return result;
        }
    }

    result.reason = "IterationLimit";
    return result;
}


constexpr std::array<double, 6> width_ratios{
    0.5, 0.75, 1.0, 1.25, 1.5, 2.0
};

template<class Model>
void analyze_widths(
    const Model& model,
    const std::string& name,
    const VanelessDiffuserGeometry& reference,
    std::ofstream& csv)
{
    const auto g = geometry();
    const auto s = scales();

    for (double ratio : width_ratios) {
        auto diffuser = reference;
        diffuser.outlet_width =
            reference.outlet_width * ratio;

        std::cout
            << "\n" << name
            << " width_ratio=" << ratio
            << " width_m=" << diffuser.outlet_width
            << '\n';

        for (double rpm : speeds) {
            const auto recovery = find_boundary(
                model, g, diffuser, s, rpm,
                BoundaryKind::Recovery);

            const auto capacity = find_capacity(
                model, g, diffuser, s, rpm);

            const double nan =
                std::numeric_limits<double>::quiet_NaN();

            const double chi =
                recovery.converged && capacity.converged
                    ? recovery.mass / capacity.mass
                    : nan;

            csv
                << name << ','
                << ratio << ','
                << diffuser.outlet_width << ','
                << rpm << ','
                << recovery.bracketed << ','
                << recovery.converged << ','
                << recovery.mass << ','
                << recovery.residual << ','
                << recovery.bracket_width << ','
                << capacity.bracketed << ','
                << capacity.converged << ','
                << capacity.search_limited << ','
                << capacity.mass << ','
                << capacity.residual << ','
                << capacity.critical_pressure << ','
                << capacity.sonic_mach << ','
                << capacity.bracket_width << ','
                << capacity.evaluations << ','
                << capacity.reason << ','
                << chi << '\n';

            std::cout
                << "  N=" << std::setw(6) << rpm
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

            std::cout
                << " chi=" << chi
                << " capacity_reason=" << capacity.reason
                << '\n';
        }
    }
}

} // namespace

int main()
{
    try {
        const std::filesystem::path output =
            "benchmarks/results-complete-stage-characteristic/"
            "width-sensitivity.csv";

        std::filesystem::create_directories(
            output.parent_path());

        std::ofstream csv(output);

        if (!csv)
            throw std::runtime_error(
                "Cannot open width sensitivity CSV");

        csv << std::setprecision(17);

        csv
            << "model,width_ratio,diffuser_width_m,rpm,"
            << "recovery_bracketed,recovery_converged,"
            << "recovery_mass_kg_s,recovery_residual_pa,"
            << "recovery_bracket_width_kg_s,"
            << "capacity_bracketed,capacity_converged,"
            << "capacity_search_limited,"
            << "capacity_mass_kg_s,capacity_residual_kg_s,"
            << "critical_pressure_pa,sonic_mach,"
            << "capacity_bracket_width_kg_s,"
            << "capacity_evaluations,capacity_reason,"
            << "recovery_to_capacity_ratio\n";

        const thermomodel::IdealGas ideal;

        // One physical diffuser geometry for both EOS models.
        const auto reference =
            calibrate(ideal, geometry(), scales());

        analyze_widths(
            ideal, "IdealGas", reference, csv);

        const auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        thermomodel::PRParameters parameters;
        parameters.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        const thermomodel::PengRobinson pr(
            caloric, parameters);

        analyze_widths(
            pr, "PengRobinson", reference, csv);

        if (!csv)
            throw std::runtime_error(
                "Width sensitivity CSV write failed");

        std::cout
            << "\nCSV: " << output << '\n';

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
