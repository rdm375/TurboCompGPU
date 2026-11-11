#include "lin2026_natural_gas.hpp"

#include <array>
#include <iomanip>
#include <iostream>

int main()
{
    const auto gas = turbocompgpu::benchmark_data::lin2026_natural_gas();
    constexpr double reference_temperature_k = 298.15;
    constexpr double reference_pressure_pa = 101325.0;
    const auto reference = gas.evaluate_Tp(
        reference_temperature_k, reference_pressure_pa);
    if (!reference) {
        std::cerr << "PR natural-gas enthalpy reference evaluation failed: "
                  << reference.message << '\n';
        return 1;
    }
    constexpr std::array pressures_mpa{6.5, 6.95, 7.1, 7.5, 8.0};
    constexpr std::array temperatures_k{293.15, 298.15, 303.15, 308.15};
    std::cout << std::setprecision(12)
              << "temperature_K,pressure_MPa,Z,density_kg_m3,enthalpy_J_kg,"
                 "enthalpy_relative_to_reference_J_kg,cp_J_kg_K,cv_J_kg_K,"
                 "cp_over_cv\n";
    for (double T : temperatures_k) {
        for (double p_mpa : pressures_mpa) {
            const auto state = gas.evaluate_Tp(T, p_mpa*1.0e6);
            if (!state) {
                std::cerr << "PR natural-gas property evaluation failed at T="
                          << T << " K, p=" << p_mpa << " MPa: "
                          << state.message << '\n';
                return 1;
            }
            std::cout << T << ',' << p_mpa << ',' << state.value.Z << ','
                      << state.value.rho << ',' << state.value.h << ','
                      << state.value.h - reference.value.h << ','
                      << state.value.cp << ','
                      << state.value.cv << ','
                      << state.value.cp/state.value.cv << '\n';
        }
    }
}
