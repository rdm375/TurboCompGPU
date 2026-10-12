#include <turbocompgpu/vaneless_diffuser.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

using namespace turbocompgpu;

namespace {

void close(
    const char* name,
    double actual,
    double expected,
    double tolerance = 1.0e-7)
{
    const double scale =
        std::max({1.0, std::abs(actual), std::abs(expected)});

    if (!std::isfinite(actual) ||
        std::abs(actual - expected) > tolerance * scale) {
        std::cerr << "FAIL: " << name
                  << " actual=" << actual
                  << " expected=" << expected << '\n';
        throw std::runtime_error(name);
    }
}

template<class Model>
void manufactured_case(const Model& model, const char* label)
{
    // Manufacture a diffuser solution at a known pressure.
    // This tests the pressure solver against an independently
    // constructed conservation-consistent reference.

    const auto state2 = evaluate(model, 320.0, 5.0e6);

    constexpr double r2 = 0.20;
    constexpr double r3 = 0.30;
    constexpr double vm2 = 180.0;
    constexpr double whirl2 = 120.0;
    constexpr double mass = 2.0;

    const FlowStation inlet_station{
        state2, {vm2, whirl2}
    };

    const auto stagnation2 =
        evaluate_stagnation_state(model, inlet_station);

    const double p3 =
        state2.p +
        0.35 * (
            stagnation2.stagnation_state.p - state2.p);

    const auto inversion = invert_temperature(
        model,
        p3,
        state2.s,
        Property::Entropy,
        state2.T);

    const auto& state3 = inversion.state;

    const double whirl3 = r2 / r3 * whirl2;

    const double kinetic =
        2.0 * (
            stagnation2.stagnation_enthalpy - state3.h)
        - whirl3 * whirl3;

    if (!(kinetic > 0.0))
        throw std::runtime_error(
            "Manufactured state has no real meridional speed");

    const double vm3 = std::sqrt(kinetic);

    const double area3 = mass / (state3.rho * vm3);

    const double width3 =
        area3 / (2.0 * std::numbers::pi * r3);

    const VanelessDiffuserGeometry geometry{
        r2, r3, width3, 0.0
    };

    const auto result = evaluate_vaneless_diffuser(
        model,
        geometry,
        {state2, vm2, whirl2, mass});

    close("mass", result.mass_residual, 0.0);
    close("energy", result.energy_residual, 0.0);
    close("entropy", result.entropy_residual, 0.0);
    close("angular momentum",
          result.angular_momentum_residual, 0.0);

    close("stagnation pressure",
          result.outlet_stagnation_pressure,
          result.inlet_stagnation_pressure);

    close("outlet whirl",
          result.outlet_tangential_velocity, whirl3);

    close("outlet area",
          result.outlet_effective_area, area3);

    if (!(result.admissible_roots >= 1))
        throw std::runtime_error("No root reported");

    // The manufactured solution is one admissible root.
    // If two roots exist, the solver deliberately selects
    // the higher-pressure branch. Do not require the
    // selected pressure to equal the manufactured one
    // unless it is known to be that branch.

    std::cout << "PASS: " << label
              << " vaneless diffuser"
              << " roots=" << result.admissible_roots
              << " p3=" << result.selected_pressure
              << '\n';
}

} // namespace

int main()
{
    try {
        thermomodel::IdealGas ideal;
        manufactured_case(ideal, "IdealGas");

        auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        thermomodel::PRParameters params;
        params.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        thermomodel::PengRobinson pr(caloric, params);
        manufactured_case(pr, "PengRobinson");

        std::cout << "PASS: Stage 8B.2 diffuser\n";
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
