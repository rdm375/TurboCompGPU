#include <turbocompgpu/impeller_system_loss.hpp>
#include <turbocompgpu/reference_shift_model.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace {

using namespace turbocompgpu;

constexpr double pi = 3.141592653589793238462643383279502884;

void check_scaled(
    const char* name,
    double residual,
    double scale,
    double tolerance)
{
    if (!(scale > 0.0) ||
        !std::isfinite(residual) ||
        std::abs(residual) > tolerance * scale) {
        std::cerr << std::setprecision(17)
                  << "FAIL: " << name
                  << " residual=" << residual
                  << " scale=" << scale
                  << " normalized=" << residual / scale
                  << '\n';
        throw std::runtime_error(name);
    }
}

ImpellerSystemGeometry geometry()
{
    return {
        {0.10, 0.20, 0.05},
        {0.40, 0.025, 0.10},
        {0.40, 45.0, 0.15}
    };
}

ImpellerTrialState trial()
{
    return {
        300.0, 500000.0,
        340.0, 700000.0,
        5.0, 5.0,
        12000.0, 0.90,
        20.0, 180.0,
        200.0
    };
}

thermomodel::State evaluate(
    const thermomodel::Model& model,
    double T,
    double pressure)
{
    const auto result = model.evaluate_Tp(T, pressure);

    if (!result) {
        throw std::runtime_error(
            "Thermodynamic evaluation failed: " +
            result.message);
    }

    if (!(result.value.rho > 0.0) ||
        !std::isfinite(result.value.h) ||
        !std::isfinite(result.value.s)) {
        throw std::runtime_error(
            "Invalid thermodynamic state");
    }

    return result.value;
}

struct ManufacturedContext {
    double inlet_h0;
    double inlet_s;
    double inlet_whirl;
    double mass_flow;
    double omega;
    double U1;
    double U2;
    double A2;
    double slip_factor;
    double cot_beta;
    double prescribed_entropy;
};

ManufacturedContext make_context(
    const thermomodel::Model& model,
    const ImpellerSystemGeometry& g,
    const ImpellerTrialState& t,
    double entropy_generation)
{
    const auto inlet = evaluate(
        model, t.inlet_temperature, t.inlet_pressure);

    const double A1 =
        0.25 * pi *
        (g.inlet.shroud_diameter * g.inlet.shroud_diameter -
         g.inlet.hub_diameter * g.inlet.hub_diameter) *
        (1.0 - g.inlet.blockage);

    const double A2 =
        pi * g.exit.diameter * g.exit.width *
        (1.0 - g.exit.blockage);

    if (!(A1 > 0.0) || !(A2 > 0.0)) {
        throw std::runtime_error("Invalid effective flow area");
    }

    const double Vm1 = t.inlet_mass_flow / (inlet.rho * A1);
    const double omega = 2.0 * pi * t.rotational_speed_rpm / 60.0;

    const double beta =
        g.impeller.exit_blade_angle_deg * pi / 180.0;

    return {
        inlet.h +
            0.5 * (Vm1 * Vm1 +
                   t.inlet_prewhirl * t.inlet_prewhirl),
        inlet.s,
        t.inlet_prewhirl,
        t.inlet_mass_flow,
        omega,
        omega * g.impeller.inlet_diameter / 2.0,
        omega * g.impeller.exit_diameter / 2.0,
        A2,
        t.slip_factor,
        1.0 / std::tan(beta),
        entropy_generation
    };
}

double energy_defect(
    const thermomodel::State& outlet,
    const ManufacturedContext& c)
{
    const double Vm2 =
        c.mass_flow / (outlet.rho * c.A2);

    const double whirl =
        c.slip_factor *
        (c.U2 - Vm2 * c.cot_beta);

    const double h02 =
        outlet.h + 0.5 * (Vm2 * Vm2 + whirl * whirl);

    const double euler_work =
        c.U2 * whirl - c.U1 * c.inlet_whirl;

    return h02 - c.inlet_h0 - euler_work;
}

// Solve the entropy equation for pressure at fixed T.
// Pressure is bracketed and bisected in logarithmic coordinates.
double pressure_from_entropy(
    const thermomodel::Model& model,
    const ManufacturedContext& c,
    double T)
{
    const double target_s =
        c.inlet_s + c.prescribed_entropy;

    // Bounds chosen for this regression scenario, not a
    // general-purpose thermodynamic inversion routine.
    double lo = std::log(50000.0);
    double hi = std::log(30000000.0);

    auto entropy_defect = [&](double log_pressure) {
        return evaluate(model, T, std::exp(log_pressure)).s
             - target_s;
    };

    double flo = entropy_defect(lo);
    double fhi = entropy_defect(hi);

    if (!(flo > 0.0 && fhi < 0.0)) {
        throw std::runtime_error(
            "Entropy pressure root is not bracketed");
    }

    for (int iteration = 0; iteration < 100; ++iteration) {
        const double mid = lo + 0.5 * (hi - lo);
        const double fm = entropy_defect(mid);

        if (fm > 0.0) {
            lo = mid;
            flo = fm;
        } else {
            hi = mid;
            fhi = fm;
        }

        if (hi - lo <= 2.0e-14) {
            return std::exp(lo + 0.5 * (hi - lo));
        }
    }

    throw std::runtime_error(
        "Entropy pressure solve did not converge");
}

struct OutletRoot {
    double T;
    double pressure;
    double energy_defect;
    double entropy_defect;
};

OutletRoot manufacture_outlet(
    const thermomodel::Model& model,
    const ManufacturedContext& c)
{
    auto at_temperature = [&](double T) {
        const double pressure =
            pressure_from_entropy(model, c, T);

        const auto outlet = evaluate(model, T, pressure);

        return OutletRoot {
            T,
            pressure,
            energy_defect(outlet, c),
            outlet.s - c.inlet_s - c.prescribed_entropy
        };
    };

    double lo = 280.0;
    double hi = 420.0;

    auto left = at_temperature(lo);
    auto right = at_temperature(hi);

    if (!(left.energy_defect < 0.0 &&
          right.energy_defect > 0.0)) {
        throw std::runtime_error(
            "Energy temperature root is not bracketed");
    }

    for (int iteration = 0; iteration < 100; ++iteration) {
        const double mid = lo + 0.5 * (hi - lo);
        const auto candidate = at_temperature(mid);

        if (candidate.energy_defect < 0.0) {
            lo = mid;
            left = candidate;
        } else {
            hi = mid;
            right = candidate;
        }

        if (hi - lo <= 2.0e-13 * std::max(1.0, mid)) {
            return at_temperature(lo + 0.5 * (hi - lo));
        }
    }

    throw std::runtime_error(
        "Energy temperature solve did not converge");
}

struct ManufacturedResult {
    ImpellerTrialState trial_state;
    ImpellerSystemWithLossResult result;
};

ManufacturedResult manufacture(
    const thermomodel::Model& model,
    double entropy_generation)
{
    const auto g = geometry();
    auto t = trial();

    t.outlet_mass_flow = t.inlet_mass_flow;

    const auto c =
        make_context(model, g, t, entropy_generation);

    const auto root = manufacture_outlet(model, c);

    t.outlet_temperature = root.T;
    t.outlet_pressure = root.pressure;

    const auto outlet = evaluate(
        model, t.outlet_temperature, t.outlet_pressure);

    const double Vm2 =
        c.mass_flow / (outlet.rho * c.A2);

    t.outlet_whirl =
        c.slip_factor * (c.U2 - Vm2 * c.cot_beta);

    const double r1 = g.impeller.inlet_diameter / 2.0;
    const double r2 = g.impeller.exit_diameter / 2.0;

    t.shaft_torque =
        c.mass_flow *
        (r2 * t.outlet_whirl - r1 * t.inlet_prewhirl);

    return {
        t,
        evaluate_impeller_system_with_loss(
            model, g, t, {entropy_generation})
    };
}

void verify(
    const ManufacturedResult& manufactured,
    double entropy_generation,
    const char* label)
{
    const auto& t = manufactured.trial_state;
    const auto& r = manufactured.result;

    const double omega =
        2.0 * pi * t.rotational_speed_rpm / 60.0;

    const double mass_scale =
        std::max(1.0, t.inlet_mass_flow);

    const double torque_scale =
        std::max(1.0, std::abs(t.shaft_torque));

    const double power_scale =
        std::max(1.0, std::abs(t.shaft_torque * omega));

    const double whirl_scale =
        std::max(1.0, std::abs(t.outlet_whirl));

    const double entropy_scale =
        std::max(1.0, std::abs(entropy_generation));

    check_scaled(
        "mass", r.system.conservation.mass,
        mass_scale, 1.0e-12);

    check_scaled(
        "angular momentum",
        r.system.conservation.angular_momentum,
        torque_scale, 1.0e-12);

    check_scaled(
        "energy", r.system.conservation.energy,
        power_scale, 1.0e-10);

    check_scaled(
        "slip", r.system.slip_residual,
        whirl_scale, 1.0e-12);

    check_scaled(
        "entropy", r.loss.entropy_residual,
        entropy_scale, 1.0e-10);

    // The reversible limit has exactly zero entropy generation.
    // Floating-point subtraction may produce a small negative
    // entropy difference even when the entropy residual converges.
    //
    // Keep the independent entropy-residual tolerance unchanged.

    constexpr double entropy_admissibility_tolerance = 1.0e-10;

    if (r.loss.actual_entropy_increase <
        -entropy_admissibility_tolerance) {
        throw std::runtime_error(
            "Manufactured state violates entropy inequality");
    }

    std::cout << std::setprecision(17)
              << label << '\n'
              << "  T2               = "
              << t.outlet_temperature << " K\n"
              << "  p2               = "
              << t.outlet_pressure << " Pa\n"
              << "  Vtheta2          = "
              << t.outlet_whirl << " m/s\n"
              << "  shaft torque     = "
              << t.shaft_torque << " Nm\n"
              << "  p02/p01          = "
              << r.system.outlet_stagnation.stagnation_state.p /
                 r.system.inlet_stagnation.stagnation_state.p
              << '\n'
              << "  mass residual    = "
              << r.system.conservation.mass << '\n'
              << "  torque residual  = "
              << r.system.conservation.angular_momentum << '\n'
              << "  energy residual  = "
              << r.system.conservation.energy << '\n'
              << "  slip residual    = "
              << r.system.slip_residual << '\n'
              << "  entropy residual = "
              << r.loss.entropy_residual << '\n';

    std::cout << "PASS: " << label << '\n';
}

void test_model(
    const thermomodel::Model& model,
    const char* name)
{
    const auto baseline = manufacture(model, 20.0);

    verify(baseline, 20.0, name);

    // Additional nonnegative loss targets test the
    // coupled thermodynamic/aerodynamic solution.
    const auto zero_loss = manufacture(model, 0.0);
    const auto higher_loss = manufacture(model, 40.0);

    verify(zero_loss, 0.0, "Zero-loss manufactured state");
    verify(higher_loss, 40.0, "Higher-loss manufactured state");

    // Do not impose the fixed-work pressure-penalty
    // formula here: Euler work changes with outlet density.
    const auto& a = baseline.result.system;
    const auto& b = higher_loss.result.system;

    if (!std::isfinite(a.euler_specific_work) ||
        !std::isfinite(b.euler_specific_work)) {
        throw std::runtime_error(
            "Non-finite coupled Euler work");
    }
}


void test_full_reference_invariance(
    const thermomodel::Model& model,
    const char* label)
{
    constexpr double Ch = 1.0e6;
    constexpr double Cs = 1.0e4;
    constexpr double entropy_generation = 20.0;

    const ReferenceShiftModel shifted_model(
        model, Ch, Cs);

    const auto baseline =
        manufacture(model, entropy_generation);

    const auto shifted =
        manufacture(shifted_model, entropy_generation);

    // Verify the shifted manufactured state independently.
    verify(
        shifted,
        entropy_generation,
        "Reference-shifted manufactured state");

    const auto& a = baseline.trial_state;
    const auto& b = shifted.trial_state;

    const auto& ra = baseline.result;
    const auto& rb = shifted.result;

    auto compare = [](
        const char* name,
        double actual,
        double expected,
        double scale,
        double tolerance)
    {
        check_scaled(
            name,
            actual - expected,
            scale,
            tolerance);
    };

    compare(
        "Reference-invariant outlet temperature",
        b.outlet_temperature,
        a.outlet_temperature,
        std::max(1.0, std::abs(a.outlet_temperature)),
        1.0e-10);

    compare(
        "Reference-invariant outlet pressure",
        b.outlet_pressure,
        a.outlet_pressure,
        std::max(1.0, std::abs(a.outlet_pressure)),
        1.0e-10);

    compare(
        "Reference-invariant outlet whirl",
        b.outlet_whirl,
        a.outlet_whirl,
        std::max(1.0, std::abs(a.outlet_whirl)),
        1.0e-10);

    compare(
        "Reference-invariant shaft torque",
        b.shaft_torque,
        a.shaft_torque,
        std::max(1.0, std::abs(a.shaft_torque)),
        1.0e-10);

    compare(
        "Reference-invariant Euler work",
        rb.system.euler_specific_work,
        ra.system.euler_specific_work,
        std::max(
            1.0,
            std::abs(ra.system.euler_specific_work)),
        1.0e-10);

    const double pressure_ratio_a =
        ra.system.outlet_stagnation.stagnation_state.p /
        ra.system.inlet_stagnation.stagnation_state.p;

    const double pressure_ratio_b =
        rb.system.outlet_stagnation.stagnation_state.p /
        rb.system.inlet_stagnation.stagnation_state.p;

    compare(
        "Reference-invariant stagnation pressure ratio",
        pressure_ratio_b,
        pressure_ratio_a,
        std::max(1.0, std::abs(pressure_ratio_a)),
        1.0e-10);

    // Absolute thermodynamic properties should shift
    // by the prescribed reference constants.

    compare(
        "Inlet stagnation enthalpy reference shift",
        rb.system.inlet_stagnation.stagnation_enthalpy -
            ra.system.inlet_stagnation.stagnation_enthalpy,
        Ch,
        std::abs(Ch),
        1.0e-10);

    compare(
        "Outlet stagnation enthalpy reference shift",
        rb.system.outlet_stagnation.stagnation_enthalpy -
            ra.system.outlet_stagnation.stagnation_enthalpy,
        Ch,
        std::abs(Ch),
        1.0e-10);

    compare(
        "Inlet stagnation entropy reference shift",
        rb.system.inlet_stagnation.stagnation_state.s -
            ra.system.inlet_stagnation.stagnation_state.s,
        Cs,
        std::abs(Cs),
        1.0e-10);

    compare(
        "Outlet stagnation entropy reference shift",
        rb.system.outlet_stagnation.stagnation_state.s -
            ra.system.outlet_stagnation.stagnation_state.s,
        Cs,
        std::abs(Cs),
        1.0e-10);

    // Conservation and closure residuals must remain
    // invariant within the established numerical accuracy.

    compare(
        "Reference-invariant mass residual",
        rb.system.conservation.mass,
        ra.system.conservation.mass,
        std::max(1.0, a.inlet_mass_flow),
        1.0e-12);

    compare(
        "Reference-invariant torque residual",
        rb.system.conservation.angular_momentum,
        ra.system.conservation.angular_momentum,
        std::max(1.0, std::abs(a.shaft_torque)),
        1.0e-12);

    const double omega =
        2.0 * pi * a.rotational_speed_rpm / 60.0;

    compare(
        "Reference-invariant energy residual",
        rb.system.conservation.energy,
        ra.system.conservation.energy,
        std::max(1.0, std::abs(a.shaft_torque * omega)),
        1.0e-10);

    compare(
        "Reference-invariant slip residual",
        rb.system.slip_residual,
        ra.system.slip_residual,
        std::max(1.0, std::abs(a.outlet_whirl)),
        1.0e-12);

    compare(
        "Reference-invariant entropy residual",
        rb.loss.entropy_residual,
        ra.loss.entropy_residual,
        std::max(1.0, entropy_generation),
        1.0e-10);

    std::cout
        << "PASS: Full reference invariance: "
        << label << '\n';
}

} // namespace

int main()
{
    try {
        const thermomodel::IdealGas ideal_gas;

        const auto caloric =
            std::make_shared<thermomodel::ConstantCp>(
                2200.0);

        const thermomodel::PengRobinson methane(caloric);

        test_model(ideal_gas, "Ideal-gas five-residual state");
        test_model(methane, "Peng-Robinson five-residual state");

        test_full_reference_invariance(
            ideal_gas,
            "Ideal gas");

        test_full_reference_invariance(
            methane,
            "Peng-Robinson methane");

        std::cout
            << "\nPASS: All five-residual manufactured states\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
