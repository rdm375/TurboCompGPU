#include <turbocompgpu/compressor.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace {

class ShiftedReferenceModel final : public thermomodel::Model {
public:
    ShiftedReferenceModel(
        const thermomodel::Model& base,
        double enthalpy_shift,
        double entropy_shift)
        : base_(base),
          enthalpy_shift_(enthalpy_shift),
          entropy_shift_(entropy_shift)
    {
    }

    thermomodel::Result<thermomodel::State> evaluate_Tp(
        double T,
        double pressure,
        const thermomodel::Composition& composition = {}) const override
    {
        auto result = base_.evaluate_Tp(T, pressure, composition);

        if (result) {
            result.value.h += enthalpy_shift_;
            result.value.s += entropy_shift_;
        }

        return result;
    }

    thermomodel::Result<thermomodel::State> evaluate_Trho(
        double T,
        double density,
        const thermomodel::Composition& composition = {}) const override
    {
        auto result = base_.evaluate_Trho(T, density, composition);

        if (result) {
            result.value.h += enthalpy_shift_;
            result.value.s += entropy_shift_;
        }

        return result;
    }

    thermomodel::Capabilities capabilities() const noexcept override
    {
        return base_.capabilities();
    }

    thermomodel::Metadata metadata() const override
    {
        return base_.metadata();
    }

private:
    const thermomodel::Model& base_;
    double enthalpy_shift_;
    double entropy_shift_;
};

void check_close(
    const char* name,
    double actual,
    double expected,
    double relative_tolerance = 1.0e-9)
{
    const double scale =
        std::max(1.0, std::abs(expected));

    const double error =
        std::abs(actual - expected) / scale;

    if (!std::isfinite(actual) ||
        error > relative_tolerance) {
        std::cerr
            << "FAIL: " << name
            << " actual=" << actual
            << " expected=" << expected
            << " relative_error=" << error
            << '\n';

        throw std::runtime_error(
            "Reference-shift invariance failure");
    }
}

void test_model(
    const char* name,
    const thermomodel::Model& base)
{
    const turbocompgpu::OperatingPoint input{
        300.0,
        5.0e6,
        10.0e6,
        1.0,
        0.8
    };

    const auto reference =
        turbocompgpu::solve_compressor(base, input);

    struct Shift {
        double h;
        double s;
    };

    constexpr std::array<Shift, 5> shifts{{
        {0.0, 0.0},
        {1.0e3, 1.0e2},
        {-1.0e3, -1.0e2},
        {1.0e6, 1.0e4},
        {-1.0e6, -1.0e4}
    }};

    std::cout << "\nModel: " << name << '\n';

    std::cout
        << std::setw(15) << "h shift"
        << std::setw(15) << "s shift"
        << std::setw(20) << "dT02s (K)"
        << std::setw(20) << "dT02 (K)"
        << std::setw(20) << "dWork (J/kg)"
        << '\n';

    for (const auto& shift : shifts) {
        ShiftedReferenceModel model(
            base, shift.h, shift.s);

        const auto result =
            turbocompgpu::solve_compressor(model, input);

        const double delta_T02s =
            result.outlet_isentropic.T -
            reference.outlet_isentropic.T;

        const double delta_T02 =
            result.outlet.T -
            reference.outlet.T;

        const double delta_work =
            result.specific_work -
            reference.specific_work;

        std::cout
            << std::setw(15) << shift.h
            << std::setw(15) << shift.s
            << std::setw(20) << delta_T02s
            << std::setw(20) << delta_T02
            << std::setw(20) << delta_work
            << '\n';

        check_close(
            "Isentropic outlet temperature",
            result.outlet_isentropic.T,
            reference.outlet_isentropic.T);

        check_close(
            "Actual outlet temperature",
            result.outlet.T,
            reference.outlet.T);

        check_close(
            "Specific work",
            result.specific_work,
            reference.specific_work);

        check_close(
            "Power",
            result.power,
            reference.power);
    }

    std::cout
        << "PASS: " << name
        << " reference-shift invariance\n";
}

} // namespace

int main()
{
    try {
        std::cout << std::setprecision(12);

        thermomodel::IdealGas ideal_gas;

        test_model("Ideal gas", ideal_gas);

        auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        thermomodel::PRParameters parameters;
        parameters.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        thermomodel::PengRobinson peng_robinson(
            caloric, parameters);

        test_model("Peng-Robinson", peng_robinson);

        std::cout
            << "\nPASS: All reference-shift tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
