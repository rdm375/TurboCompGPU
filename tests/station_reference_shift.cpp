#include <turbocompgpu/station.hpp>

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
        auto result =
            base_.evaluate_Tp(T, pressure, composition);

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
        auto result =
            base_.evaluate_Trho(T, density, composition);

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
    double relative_tolerance = 1.0e-10)
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

        throw std::runtime_error(name);
    }
}

void test_model(
    const char* name,
    const thermomodel::Model& base)
{
    using namespace turbocompgpu;

    const auto static_state =
        evaluate(base, 300.0, 5.0e6);

    const AbsoluteVelocity velocity{
        120.0,
        160.0
    };

    const auto reference =
        evaluate_stagnation_state(
            base,
            {static_state, velocity});

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
        << std::setw(20) << "dT0 (K)"
        << std::setw(20) << "dp0 (Pa)"
        << '\n';

    for (const auto& shift : shifts) {
        ShiftedReferenceModel model(
            base, shift.h, shift.s);

        const auto shifted_static =
            evaluate(model, 300.0, 5.0e6);

        const auto result =
            evaluate_stagnation_state(
                model,
                {shifted_static, velocity});

        const double delta_T =
            result.stagnation_state.T -
            reference.stagnation_state.T;

        const double delta_p =
            result.stagnation_state.p -
            reference.stagnation_state.p;

        std::cout
            << std::setw(15) << shift.h
            << std::setw(15) << shift.s
            << std::setw(20) << delta_T
            << std::setw(20) << delta_p
            << '\n';

        check_close(
            "Stagnation temperature invariance",
            result.stagnation_state.T,
            reference.stagnation_state.T);

        check_close(
            "Stagnation pressure invariance",
            result.stagnation_state.p,
            reference.stagnation_state.p);

        check_close(
            "Stagnation enthalpy shift",
            result.stagnation_enthalpy -
                reference.stagnation_enthalpy,
            shift.h,
            1.0e-9);

        check_close(
            "Stagnation entropy shift",
            result.stagnation_entropy -
                reference.stagnation_entropy,
            shift.s,
            1.0e-9);
    }

    std::cout
        << "PASS: " << name
        << " station reference-shift invariance\n";
}

} // namespace

int main()
{
    try {
        std::cout << std::setprecision(17);

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
            << "\nPASS: All station reference-shift tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
