#include <turbocompgpu/vaneless_diffuser_branches.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>

using namespace turbocompgpu;

namespace {

void check(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void close(
    double actual,
    double expected,
    double tolerance,
    const char* message)
{
    const double scale =
        std::max({1.0, std::abs(actual), std::abs(expected)});

    if (!std::isfinite(actual) ||
        std::abs(actual - expected) > tolerance * scale) {
        throw std::runtime_error(message);
    }
}

template <typename ModelType>
void run_case(
    const ModelType& model,
    const char* name)
{
    const ImpellerSystemGeometry geometry{
        {0.10, 0.20, 0.05},
        {0.40, 0.025, 0.10},
        {0.40, 45.0, 0.15}
    };

    const NonlinearResidualScales scales{
        50000.0, 100.0, 500000.0
    };

    const OperatingPointRequest request{
        {
            300.0, 500000.0, 5.0, 12000.0,
            0.90, 20.0, 20.0
        },
        {
            340.0,
            std::log(700000.0 / scales.pressure)
        }
    };

    NonlinearSolverOptions options;
    options.residual_tolerance = 1.0e-13;

    const auto impeller = evaluate_operating_point(
        model, geometry, request, scales, options);

    check(
        impeller.converged &&
        impeller.independently_validated &&
        impeller.basic_physics_valid,
        "Impeller operating point failed");

    const auto state2 = evaluate(
        model,
        impeller.outlet_temperature,
        impeller.outlet_pressure);

    const double mass = request.operating.mass_flow;

    const double area2 = effective_flow_area(
        impeller_exit_area(geometry.exit),
        geometry.exit.blockage);

    const double vm2 =
        mass / (state2.rho * area2);

    const double omega =
        2.0 * std::numbers::pi *
        request.operating.rotational_speed_rpm / 60.0;

    const double r2 = geometry.exit.diameter / 2.0;
    const double r3 = 1.5 * r2;

    const double beta =
        geometry.impeller.exit_blade_angle_deg *
        std::numbers::pi / 180.0;

    const double whirl2 =
        request.operating.slip_factor *
        (omega * r2 - vm2 / std::tan(beta));

    const double whirl3 = r2 * whirl2 / r3;

    const auto stagnation = evaluate_stagnation_state(
        model, {state2, {vm2, whirl2}});

    const double p_reference =
        state2.p + 0.10 *
        (stagnation.stagnation_state.p - state2.p);

    const auto reference = invert_temperature(
        model, p_reference, state2.s,
        Property::Entropy, state2.T);

    const double vm_reference = std::sqrt(
        2.0 * (
            stagnation.stagnation_enthalpy -
            reference.state.h) -
        whirl3 * whirl3);

    const double area3 =
        mass / (reference.state.rho * vm_reference);

    const VanelessDiffuserGeometry diffuser_geometry{
        r2,
        r3,
        area3 / (2.0 * std::numbers::pi * r3),
        0.0
    };

    const VanelessDiffuserInlet inlet{
        state2, vm2, whirl2, mass
    };

    const auto subsonic = evaluate_vaneless_diffuser_branch(
        model, model, diffuser_geometry, inlet,
        DiffuserBranch::Subsonic);

    const auto supersonic = evaluate_vaneless_diffuser_branch(
        model, model, diffuser_geometry, inlet,
        DiffuserBranch::Supersonic);

    check(
        subsonic.status == DiffuserBranchStatus::Converged,
        "Subsonic branch did not converge");

    check(
        supersonic.status == DiffuserBranchStatus::Converged,
        "Supersonic branch did not converge");

    check(
        subsonic.outlet.meridional_mach < 1.0,
        "Subsonic Mach classification failed");

    check(
        supersonic.outlet.meridional_mach > 1.0,
        "Supersonic Mach classification failed");

    check(
        supersonic.outlet.pressure <
        subsonic.critical.pressure,
        "Supersonic pressure ordering failed");

    check(
        subsonic.outlet.pressure >
        subsonic.critical.pressure,
        "Subsonic pressure ordering failed");

    close(
        subsonic.critical.meridional_mach,
        1.0, 1.0e-8,
        "Critical Mach mismatch");

    close(
        subsonic.critical_mass_flow,
        supersonic.critical_mass_flow,
        1.0e-10,
        "Critical capacity mismatch");

    close(
        subsonic.mass_residual,
        0.0, 1.0e-8,
        "Subsonic mass residual");

    close(
        supersonic.mass_residual,
        0.0, 1.0e-8,
        "Supersonic mass residual");

    close(
        subsonic.energy_residual,
        0.0, 1.0e-7,
        "Subsonic energy residual");

    close(
        supersonic.energy_residual,
        0.0, 1.0e-7,
        "Supersonic energy residual");

    const auto legacy = evaluate_vaneless_diffuser(
        model, diffuser_geometry, inlet);

    close(
        subsonic.outlet.pressure,
        legacy.outlet_static.p,
        1.0e-8,
        "Legacy high-pressure solution mismatch");

    // Above-critical request with the same fixed inlet
    // thermodynamic and velocity state. This is a capacity
    // query, not a new coupled impeller operating point.

    auto excessive = inlet;
    excessive.mass_flow =
        1.05 * subsonic.critical_mass_flow;

    const auto above = evaluate_vaneless_diffuser_branch(
        model, model, diffuser_geometry, excessive,
        DiffuserBranch::Subsonic);

    check(
        above.status ==
            DiffuserBranchStatus::AboveCriticalCapacity,
        "Above-critical request not detected");


    // --------------------------------------------------------
    // Stage 8B.4d: exact critical capacity.
    // --------------------------------------------------------

    auto critical_inlet = inlet;
    critical_inlet.mass_flow =
        subsonic.critical_mass_flow;

    const auto exact_critical =
        evaluate_vaneless_diffuser_branch(
            model, model, diffuser_geometry,
            critical_inlet, DiffuserBranch::Subsonic);

    check(
        exact_critical.status ==
            DiffuserBranchStatus::Critical,
        "Exact-critical status failed");

    close(
        exact_critical.outlet.meridional_mach,
        1.0, 1.0e-8,
        "Exact-critical Mach failed");

    close(
        exact_critical.outlet.pressure,
        subsonic.critical.pressure,
        1.0e-9,
        "Exact-critical pressure failed");

    // --------------------------------------------------------
    // Near-critical roots must converge toward sonic pressure.
    //
    // The two roots become increasingly close as the mass
    // deficit approaches zero.
    // --------------------------------------------------------

    const double deficits[] = {
        1.0e-2,
        1.0e-3,
        1.0e-4,
        1.0e-5
    };

    double previous_sub_distance =
        std::numeric_limits<double>::infinity();

    double previous_super_distance =
        std::numeric_limits<double>::infinity();

    for (double deficit : deficits) {
        auto near_inlet = inlet;

        near_inlet.mass_flow =
            subsonic.critical_mass_flow *
            (1.0 - deficit);

        const auto near_sub =
            evaluate_vaneless_diffuser_branch(
                model, model, diffuser_geometry,
                near_inlet, DiffuserBranch::Subsonic);

        const auto near_super =
            evaluate_vaneless_diffuser_branch(
                model, model, diffuser_geometry,
                near_inlet, DiffuserBranch::Supersonic);

        check(
            near_sub.status ==
                DiffuserBranchStatus::Converged,
            "Near-critical subsonic solve failed");

        check(
            near_super.status ==
                DiffuserBranchStatus::Converged,
            "Near-critical supersonic solve failed");

        const double sub_distance =
            near_sub.outlet.pressure -
            near_sub.critical.pressure;

        const double super_distance =
            near_super.critical.pressure -
            near_super.outlet.pressure;

        check(
            sub_distance > 0.0 &&
            super_distance > 0.0,
            "Near-critical pressure ordering failed");

        check(
            sub_distance < previous_sub_distance &&
            super_distance < previous_super_distance,
            "Near-critical roots do not converge");

        check(
            near_sub.outlet.meridional_mach < 1.0 &&
            near_super.outlet.meridional_mach > 1.0,
            "Near-critical Mach classification failed");

        check(
            near_sub.relative_pressure_bracket < 1.0e-10,
            "Subsonic pressure bracket too wide");

        check(
            near_super.relative_pressure_bracket < 1.0e-10,
            "Supersonic pressure bracket too wide");

        previous_sub_distance = sub_distance;
        previous_super_distance = super_distance;
    }

    std::cout
        << name
        << ": critical and near-critical tests PASS\n";

    std::cout << std::setprecision(12)
              << name
              << ": p_sonic=" << subsonic.critical.pressure
              << " m_critical=" << subsonic.critical_mass_flow
              << " p_sub=" << subsonic.outlet.pressure
              << " M_sub=" << subsonic.outlet.meridional_mach
              << " p_super=" << supersonic.outlet.pressure
              << " M_super=" << supersonic.outlet.meridional_mach
              << '\n';
}

} // namespace

int main()
{
    try {
        thermomodel::IdealGas ideal;
        run_case(ideal, "IdealGas");

        auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        thermomodel::PRParameters parameters;
        parameters.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        thermomodel::PengRobinson pr(
            caloric, parameters);

        run_case(pr, "PengRobinson");

        std::cout
            << "PASS: Stage 8B.4c branch-aware diffuser\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
