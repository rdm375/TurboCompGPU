#include <turbocompgpu/compressor.hpp>
#include <thermomodel/ideal_gas.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void check_close(
    const char* name,
    double actual,
    double expected,
    double relative_tolerance = 1.0e-10)
{
    const double scale = std::max(1.0, std::abs(expected));
    const double error = std::abs(actual - expected) / scale;

    if (!std::isfinite(actual) ||
        error > relative_tolerance) {
        std::cerr << "FAIL: " << name
                  << " actual=" << actual
                  << " expected=" << expected
                  << " relative_error=" << error << '\n';

        throw std::runtime_error("Regression check failed");
    }
}

} // namespace

int main()
{
    try {
        thermomodel::IdealGas model;

        constexpr double T01 = 300.0;
        constexpr double p01 = 5.0e6;
        constexpr double pressure_ratio = 2.0;
        constexpr double efficiency = 0.8;
        constexpr double mass_flow = 1.0;

        const turbocompgpu::OperatingPoint input{
            T01,
            p01,
            p01 * pressure_ratio,
            mass_flow,
            efficiency
        };

        const auto result =
            turbocompgpu::solve_compressor(model, input);

        const auto& gas = model.parameters();

        const double gamma = gas.cp / (gas.cp - gas.R);

        const double T02s_reference =
            T01 * std::pow(
                pressure_ratio,
                (gamma - 1.0) / gamma);

        const double T02_reference =
            T01 + (T02s_reference - T01) / efficiency;

        const double work_reference =
            gas.cp * (T02_reference - T01);

        check_close(
            "Isentropic outlet temperature",
            result.outlet_isentropic.T,
            T02s_reference);

        check_close(
            "Actual outlet temperature",
            result.outlet.T,
            T02_reference);

        check_close(
            "Specific work",
            result.specific_work,
            work_reference);

        check_close(
            "Power",
            result.power,
            mass_flow * work_reference);

        check_close(
            "Isentropic entropy",
            result.outlet_isentropic.s,
            result.inlet.s);

        check_close(
            "Energy balance",
            result.outlet.h - result.inlet.h,
            result.specific_work);

        const double efficiency_recovered =
            (result.outlet_isentropic.h - result.inlet.h) /
            (result.outlet.h - result.inlet.h);

        check_close(
            "Recovered isentropic efficiency",
            efficiency_recovered,
            efficiency);


        std::cout.precision(17);

        std::cout
            << "Entropy inversion:\n"
            << "  iterations    = "
            << result.entropy_inversion.iterations << '\n'
            << "  evaluations   = "
            << result.entropy_inversion.evaluations << '\n'
            << "  residual      = "
            << result.entropy_inversion.residual << '\n'
            << "  bracket width = "
            << result.entropy_inversion.bracket_width << " K\n"            << "  resolution limited = "
            << result.entropy_inversion.resolution_limited << '\n'
            << "Enthalpy inversion:\n"
            << "  iterations    = "
            << result.enthalpy_inversion.iterations << '\n'
            << "  evaluations   = "
            << result.enthalpy_inversion.evaluations << '\n'
            << "  residual      = "
            << result.enthalpy_inversion.residual << '\n'
            << "  bracket width = "
            << result.enthalpy_inversion.bracket_width << " K\n"            << "  resolution limited = "
            << result.enthalpy_inversion.resolution_limited << '\n';


        std::cout
            << "T01       = " << result.inlet.T << " K\n"
            << "T02s      = " << result.outlet_isentropic.T << " K\n"
            << "T02       = " << result.outlet.T << " K\n"
            << "Delta h0  = " << result.specific_work << " J/kg\n"
            << "Power     = " << result.power << " W\n"
            << "Efficiency= " << efficiency_recovered << '\n'
            << "PASS: Ideal-gas compressor validation\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
