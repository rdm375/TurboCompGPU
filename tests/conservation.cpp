#include <turbocompgpu/conservation.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {

using namespace turbocompgpu;

void check_close(
    const char* name,
    double actual,
    double expected,
    double tolerance = 1.0e-12)
{
    const double scale =
        std::max(1.0, std::abs(expected));

    const double error =
        std::abs(actual - expected) / scale;

    if (!std::isfinite(actual) ||
        error > tolerance) {

        std::cerr
            << "FAIL: " << name
            << " actual=" << actual
            << " expected=" << expected
            << " relative_error=" << error
            << '\n';

        throw std::runtime_error(name);
    }
}

template<class Function>
void expect_failure(
    const char* name,
    Function&& function)
{
    bool failed = false;

    try {
        function();
    }
    catch (const std::exception&) {
        failed = true;
    }

    if (!failed) {
        throw std::runtime_error(name);
    }

    std::cout << "PASS: " << name << '\n';
}

ConservationGeometry geometry()
{
    return {0.10, 0.20};
}

ConservationState consistent_state()
{
    return {
        5.0,       // inlet mass flow
        5.0,       // outlet mass flow

        20.0,      // inlet whirl
        200.0,     // outlet whirl

        300000.0,  // inlet stagnation enthalpy
        338000.0,  // outlet stagnation enthalpy

        190.0,     // shaft torque
        1000.0     // angular speed
    };
}

void test_consistent_state()
{
    const auto result =
        evaluate_conservation(
            geometry(),
            consistent_state());

    check_close("Mass conservation", result.mass, 0.0);

    check_close(
        "Angular momentum conservation",
        result.angular_momentum,
        0.0);

    check_close(
        "Stagnation energy conservation",
        result.energy,
        0.0);

    std::cout
        << "PASS: Manufactured consistent state\n";
}

void test_mass_perturbation()
{
    auto state = consistent_state();

    state.outlet_mass_flow += 0.5;

    const auto result =
        evaluate_conservation(geometry(), state);

    check_close("Mass perturbation", result.mass, 0.5);

    // The other balances must respond to the altered
    // outlet flux; they are not independently forced to zero.
    check_close(
        "Angular momentum response",
        result.angular_momentum,
        -20.0);

    check_close(
        "Energy response",
        result.energy,
        169000.0);

    std::cout
        << "PASS: Mass-flow perturbation\n";
}

void test_whirl_perturbation()
{
    auto state = consistent_state();

    state.outlet_whirl += 10.0;

    const auto result =
        evaluate_conservation(geometry(), state);

    check_close("Whirl mass residual", result.mass, 0.0);

    check_close(
        "Whirl angular momentum residual",
        result.angular_momentum,
        -10.0);

    check_close(
        "Whirl energy residual",
        result.energy,
        0.0);

    std::cout
        << "PASS: Outlet-whirl perturbation\n";
}

void test_enthalpy_perturbation()
{
    auto state = consistent_state();

    state.outlet_stagnation_h += 1000.0;

    const auto result =
        evaluate_conservation(geometry(), state);

    check_close(
        "Enthalpy mass residual",
        result.mass,
        0.0);

    check_close(
        "Enthalpy angular momentum residual",
        result.angular_momentum,
        0.0);

    check_close(
        "Enthalpy energy residual",
        result.energy,
        5000.0);

    std::cout
        << "PASS: Outlet-enthalpy perturbation\n";
}

void test_torque_perturbation()
{
    auto state = consistent_state();

    state.shaft_torque += 10.0;

    const auto result =
        evaluate_conservation(geometry(), state);

    check_close(
        "Torque angular momentum residual",
        result.angular_momentum,
        10.0);

    check_close(
        "Torque energy residual",
        result.energy,
        -10000.0);

    std::cout
        << "PASS: Shaft-torque perturbation\n";
}

void test_normalization()
{
    const auto scales =
        make_conservation_scales(
            5.0,
            38000.0,
            1000.0);

    check_close("Mass scale", scales.mass, 5.0);

    check_close(
        "Angular momentum scale",
        scales.angular_momentum,
        190.0);

    check_close(
        "Energy scale",
        scales.energy,
        190000.0);

    const ConservationResiduals residuals{
        0.5,
        19.0,
        19000.0
    };

    const auto normalized =
        normalize_conservation(residuals, scales);

    check_close("Normalized mass", normalized.mass, 0.1);

    check_close(
        "Normalized angular momentum",
        normalized.angular_momentum,
        0.1);

    check_close(
        "Normalized energy",
        normalized.energy,
        0.1);

    std::cout
        << "PASS: Residual normalization\n";
}


void test_euler_consistency()
{
    const auto g = geometry();
    auto state = consistent_state();

    const auto check_identity = [&]() {
        const auto residual =
            evaluate_conservation(g, state);

        const double U1 =
            state.angular_speed * g.inlet_radius;

        const double U2 =
            state.angular_speed * g.outlet_radius;

        // General identity: valid even when mass flows differ.
        const double general_expected =
            state.outlet_mass_flow *
                (state.outlet_stagnation_h -
                 U2 * state.outlet_whirl)
            -
            state.inlet_mass_flow *
                (state.inlet_stagnation_h -
                 U1 * state.inlet_whirl);

        const double combined_residual =
            residual.energy +
            state.angular_speed *
                residual.angular_momentum;

        check_close(
            "General Euler consistency identity",
            combined_residual,
            general_expected);

        // The specific-work form additionally requires
        // continuity.
        if (state.inlet_mass_flow ==
            state.outlet_mass_flow) {

            const double euler_work =
                U2 * state.outlet_whirl -
                U1 * state.inlet_whirl;

            const double specific_work_defect =
                state.outlet_stagnation_h -
                state.inlet_stagnation_h -
                euler_work;

            check_close(
                "Euler specific-work identity",
                combined_residual,
                state.inlet_mass_flow *
                    specific_work_defect);
        }
    };

    // Fully consistent state.
    check_identity();

    // Independent enthalpy perturbation.
    state.outlet_stagnation_h += 1000.0;
    check_identity();

    // Independent whirl perturbation.
    state = consistent_state();
    state.outlet_whirl += 10.0;
    check_identity();

    // Independent shaft-torque perturbation.
    // Torque must cancel from the combined residual.
    state = consistent_state();
    state.shaft_torque += 10.0;
    check_identity();

    // Deliberate continuity violation.
    // Only the general identity applies.
    state = consistent_state();
    state.outlet_mass_flow += 0.5;
    check_identity();

    std::cout
        << "PASS: Euler consistency identities\n";
}

void test_invalid_inputs()
{
    auto state = consistent_state();

    state.inlet_mass_flow = -1.0;

    expect_failure(
        "Negative mass flow",
        [&] {
            evaluate_conservation(geometry(), state);
        });

    state = consistent_state();

    state.angular_speed = 0.0;

    expect_failure(
        "Zero angular speed",
        [&] {
            evaluate_conservation(geometry(), state);
        });

    state = consistent_state();

    state.outlet_whirl =
        std::numeric_limits<double>::quiet_NaN();

    expect_failure(
        "NaN outlet whirl",
        [&] {
            evaluate_conservation(geometry(), state);
        });

    state = consistent_state();

    state.shaft_torque =
        std::numeric_limits<double>::infinity();

    expect_failure(
        "Infinite shaft torque",
        [&] {
            evaluate_conservation(geometry(), state);
        });

    expect_failure(
        "Invalid geometry",
        [&] {
            evaluate_conservation(
                {0.20, 0.10},
                consistent_state());
        });

    expect_failure(
        "Invalid normalization scale",
        [&] {
            normalize_conservation(
                {1.0, 1.0, 1.0},
                {1.0, 0.0, 1.0});
        });
}

} // namespace

int main()
{
    try {
        test_consistent_state();
        test_mass_perturbation();
        test_whirl_perturbation();
        test_enthalpy_perturbation();
        test_torque_perturbation();
        test_normalization();
        test_euler_consistency();
        test_invalid_inputs();

        std::cout
            << "\nPASS: All conservation tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
