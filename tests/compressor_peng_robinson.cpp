#include <turbocompgpu/compressor.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {

void check_close(
    const char* name,
    double actual,
    double expected,
    double relative_tolerance = 1.0e-9)
{
    const double scale = std::max(1.0, std::abs(expected));
    const double error = std::abs(actual - expected) / scale;

    if (!std::isfinite(actual) || error > relative_tolerance) {
        std::cerr << "FAIL: " << name
                  << " actual=" << actual
                  << " expected=" << expected
                  << " relative_error=" << error << '\n';
        throw std::runtime_error("Regression check failed");
    }
}

void check_vapor(
    const char* name,
    const thermomodel::State& state)
{
    if (state.phase != thermomodel::Phase::Vapor &&
        state.phase != thermomodel::Phase::Supercritical) {
        throw std::runtime_error(
            std::string(name) + ": unexpected phase");
    }
}

} // namespace

int main()
{
    try {
        auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        thermomodel::PRParameters params;
        params.root_policy = thermomodel::PRRootPolicy::Vapor;

        thermomodel::PengRobinson model(caloric, params);

        const turbocompgpu::OperatingPoint input{
            300.0,
            5.0e6,
            10.0e6,
            1.0,
            0.8
        };

        const auto result =
            turbocompgpu::solve_compressor(model, input);

        check_vapor("Inlet", result.inlet);
        check_vapor("Isentropic outlet", result.outlet_isentropic);
        check_vapor("Actual outlet", result.outlet);

        check_close(
            "Isentropic entropy",
            result.outlet_isentropic.s,
            result.inlet.s);

        check_close(
            "Actual enthalpy balance",
            result.outlet.h - result.inlet.h,
            result.specific_work);

        check_close(
            "Power",
            result.power,
            input.mass_flow * result.specific_work);

        const double efficiency_recovered =
            (result.outlet_isentropic.h - result.inlet.h) /
            (result.outlet.h - result.inlet.h);

        check_close(
            "Recovered efficiency",
            efficiency_recovered,
            input.isentropic_efficiency);

        if (!(result.outlet.T > result.outlet_isentropic.T) ||
            !(result.outlet_isentropic.T > result.inlet.T) ||
            !(result.specific_work > 0.0)) {
            throw std::runtime_error(
                "Unexpected compressor temperature or work ordering");
        }


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
            << "Z01       = " << result.inlet.Z << '\n'
            << "Z02s      = " << result.outlet_isentropic.Z << '\n'
            << "Z02       = " << result.outlet.Z << '\n'
            << "Delta h0  = " << result.specific_work << " J/kg\n"
            << "Power     = " << result.power << " W\n"
            << "Efficiency= " << efficiency_recovered << '\n'
            << "PASS: Peng-Robinson compressor validation\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
