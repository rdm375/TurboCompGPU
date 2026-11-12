#include <turbocompgpu/polytropic_efficiency.hpp>

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {

void require_near(double actual, double expected, double tolerance,
                  const char* label)
{
    if (!(std::isfinite(actual) &&
          std::abs(actual - expected) <= tolerance)) {
        std::cerr << label << ": actual=" << actual
                  << " expected=" << expected
                  << " tolerance=" << tolerance << '\n';
        std::exit(1);
    }
}

struct GasCase {
    double R;
    double cp;
    double compressibility_slope_per_pa;

    turbocompgpu::PolytropicPathProperties operator()(
        double temperature, double pressure) const
    {
        const double Z = 1.0 + compressibility_slope_per_pa * pressure;
        return {pressure / (R * temperature * Z), 0.0, cp};
    }

    double temperature_after(double T1, double p1, double p2,
                             double eta) const
    {
        const double a = R / (cp * eta);
        const double log_ratio = std::log(p2 / p1) +
            compressibility_slope_per_pa * (p2 - p1);
        return T1 * std::exp(a * log_ratio);
    }
};

void check_case(const GasCase& gas, double label)
{
    constexpr double T1 = 300.0;
    constexpr double p1 = 2.0e5;
    constexpr double p2 = 1.4e6;
    constexpr double eta = 0.82;
    const double T2 = gas.temperature_after(T1, p1, p2, eta);
    const auto result = turbocompgpu::calculate_polytropic_efficiency(
        T1, p1, gas.cp * T1, T2, p2, gas.cp * T2,
        [&gas](double T, double p) { return gas(T, p); }, 256, 56);

    require_near(result.polytropic_efficiency, eta, 2e-9,
                 label == 0.0 ? "ideal-gas eta_p" : "real-gas eta_p");
    require_near(result.endpoint_temperature_residual_k, 0.0, 2e-9,
                 "endpoint shooting residual");
    const auto coarse = turbocompgpu::calculate_polytropic_efficiency(
        T1, p1, gas.cp * T1 + 1.0e7, T2, p2, gas.cp * T2 + 1.0e7,
        [&gas](double T, double p) { return gas(T, p); }, 64, 48);
    require_near(coarse.polytropic_efficiency,
                 result.polytropic_efficiency, 2e-9,
                 "pressure-step refinement and enthalpy reference shift");
    if (!(std::isfinite(result.polytropic_head_j_kg) &&
          result.polytropic_head_j_kg > 0.0 &&
          result.pressure_steps == 256 &&
          result.efficiency_iterations == 56)) {
        std::cerr << "Invalid polytropic diagnostic result\n";
        std::exit(1);
    }
}

} // namespace

int main()
{
    constexpr double R = 287.05;
    constexpr double cp = 1004.5;
    check_case({R, cp, 0.0}, 0.0);
    // A simple compressible nonideal EOS, Z(p)=1+k p, gives a closed-form
    // constant-eta path and confirms the pressure-dependent EOS is used.
    check_case({R, cp, 1.0e-7}, 1.0);
    std::cout << "Ideal- and nonideal-gas polytropic checks passed\n";
    return 0;
}
