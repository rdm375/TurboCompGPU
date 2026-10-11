#include <turbocompgpu/compressor.hpp>
#include <thermomodel/ideal_gas.hpp>

#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

template<class Function>
void expect_failure(const char* name, Function&& function)
{
    try {
        function();
    }
    catch (const std::exception&) {
        std::cout << "PASS: " << name << '\n';
        return;
    }

    throw std::runtime_error(
        std::string("Expected failure: ") + name);
}

} // namespace

int main()
{
    try {
        thermomodel::IdealGas model;

        const turbocompgpu::OperatingPoint valid{
            300.0, 5.0e6, 10.0e6, 1.0, 0.8
        };

        auto invalid = valid;

        invalid.isentropic_efficiency = 0.0;
        expect_failure("Zero efficiency", [&] {
            turbocompgpu::solve_compressor(model, invalid);
        });

        invalid = valid;
        invalid.isentropic_efficiency = 1.1;
        expect_failure("Efficiency above unity", [&] {
            turbocompgpu::solve_compressor(model, invalid);
        });

        invalid = valid;
        invalid.outlet_p0 = invalid.inlet_p0;
        expect_failure("No pressure rise", [&] {
            turbocompgpu::solve_compressor(model, invalid);
        });

        invalid = valid;
        invalid.mass_flow = -1.0;
        expect_failure("Negative mass flow", [&] {
            turbocompgpu::solve_compressor(model, invalid);
        });

        invalid = valid;
        invalid.inlet_T0 =
            std::numeric_limits<double>::quiet_NaN();

        expect_failure("NaN inlet temperature", [&] {
            turbocompgpu::solve_compressor(model, invalid);
        });

        expect_failure("Unreachable entropy", [&] {
            turbocompgpu::invert_temperature(
                model, 5.0e6, 1.0e9,
                turbocompgpu::Property::Entropy, 300.0);
        });

        expect_failure("Unreachable enthalpy", [&] {
            turbocompgpu::invert_temperature(
                model, 5.0e6, -1.0e9,
                turbocompgpu::Property::Enthalpy, 300.0);
        });

        expect_failure("Invalid initial temperature", [&] {
            turbocompgpu::invert_temperature(
                model, 5.0e6, 1000.0,
                turbocompgpu::Property::Entropy, -10.0);
        });

        std::cout << "PASS: Compressor failure tests\n";
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
