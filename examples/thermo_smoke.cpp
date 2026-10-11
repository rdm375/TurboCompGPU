#include <thermogpu/thermo.hpp>
#include <thermomodel/ideal_gas.hpp>

#include <cmath>
#include <iostream>

int main()
{
    constexpr double T = 300.0;
    constexpr double P = 5.0e6;

    thermomodel::IdealGas model;

    const auto result = model.evaluate_Tp(T, P);

    if (!result) {
        std::cerr << "Thermodynamic evaluation failed: "
                  << result.message << '\n';
        return 1;
    }

    const auto& state = result.value;

    std::cout.precision(17);

    std::cout
        << "T   = " << state.T << " K\n"
        << "P   = " << state.p << " Pa\n"
        << "rho = " << state.rho << " kg/m3\n"
        << "h   = " << state.h << " J/kg\n"
        << "s   = " << state.s << " J/(kg K)\n"
        << "cp  = " << state.cp << " J/(kg K)\n"
        << "cv  = " << state.cv << " J/(kg K)\n"
        << "Z   = " << state.Z << '\n';

    if (!std::isfinite(state.rho) ||
        !std::isfinite(state.h) ||
        !std::isfinite(state.s) ||
        !std::isfinite(state.cp) ||
        !std::isfinite(state.cv)) {
        std::cerr << "Non-finite thermodynamic property\n";
        return 1;
    }

    if (std::abs(state.Z - 1.0) > 1.0e-12) {
        std::cerr << "Ideal-gas compressibility invariant failed\n";
        return 1;
    }


    // Independent ideal-gas consistency checks.
    constexpr double R_reference = 287.05;

    const double R_calculated = state.cp - state.cv;
    const double rho_reference = P / (R_reference * T);

    const double density_relative_error =
        std::abs(state.rho - rho_reference) / rho_reference;

    if (std::abs(R_calculated - R_reference) > 1.0e-10) {
        std::cerr << "Ideal-gas heat-capacity identity failed\n";
        return 1;
    }

    if (density_relative_error > 1.0e-12) {
        std::cerr << "Ideal-gas density identity failed: relative error = "
                  << density_relative_error << '\n';
        return 1;
    }

    std::cout << "Density relative error = "
              << density_relative_error << '\n';

    std::cout << "PASS: ThermoGPU integration\n";

    return 0;
}
