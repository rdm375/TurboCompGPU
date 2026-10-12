#include <turbocompgpu/passage_turning_sensitivity.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {

using namespace turbocompgpu;

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void close(
    double actual,
    double expected,
    double atol,
    double rtol,
    const char* message)
{
    require(
        std::isfinite(actual) &&
        std::isfinite(expected) &&
        std::abs(actual-expected) <=
            atol + rtol*std::max(std::abs(actual),std::abs(expected)),
        message);
}

const ImpellerSystemGeometry geometry{
    {0.10, 0.20, 0.05},
    {0.40, 0.025, 0.10},
    {0.40, 45.0, 0.15}
};

const ThroatGeometry throat{
    12, 0.15, 0.0035, 45.0, 0.002, 0.10
};

template<class Model>
void run_case(
    const Model& model,
    const char* name,
    double rpm)
{
    const PassageTurningChokingInputs base{
        300.0, 500000.0, rpm, 20.0, throat, 0.0
    };

    const auto curvature =
        evaluate_passage_turning_curvature(
            model, model, geometry, base);

    close(
        curvature.implicit_second_derivative_lambda,
        curvature.direct_second_derivative_lambda,
        2e-4,
        2e-2,
        "Implicit/direct curvature disagreement");

    require(
        curvature.implicit_second_derivative_lambda < 0.0,
        "Expected negative curvature for test operating point");

    std::cout
        << std::setprecision(12)
        << name
        << " rpm=" << rpm
        << " curvature_implicit="
        << curvature.implicit_second_derivative_lambda
        << " curvature_direct="
        << curvature.direct_second_derivative_lambda
        << " kg/s\n";

    for (double lambda : {0.0, 0.25, 0.5, 0.75, 1.0}) {
        auto input = base;
        input.turning_fraction = lambda;

        const auto sensitivity =
            evaluate_passage_turning_sensitivity(
                model, model, geometry, input);

        const double implicit =
            sensitivity.implicit_derivative_lambda;

        const double direct =
            sensitivity.finite_difference_derivative_lambda;

        close(
            implicit,
            direct,
            2e-5,
            2e-2,
            "Implicit/direct first derivative disagreement");

        if (lambda == 0.0) {
            require(
                std::abs(implicit) < 1e-8,
                "Zero turning-fraction derivative not zero");

            require(
                std::abs(direct) < 2e-5,
                "Direct endpoint derivative too large");
        }

        std::cout
            << std::setprecision(12)
            << "  lambda=" << lambda
            << " dmdot_dlambda_implicit=" << implicit
            << " dmdot_dlambda_direct=" << direct
            << " denominator="
            << sensitivity.implicit_denominator
            << '\n';
    }
}

} // namespace

int main()
{
    try {
        thermomodel::IdealGas air;

        thermomodel::PRParameters pr_params;
        pr_params.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        thermomodel::IdealGasParameters methane_params;
        methane_params.R =
            8.31446261815324 / pr_params.molar_mass;
        methane_params.cp = 2200.0;
        methane_params.Tref = pr_params.Tref;
        methane_params.pref = pr_params.pref;

        thermomodel::IdealGas methane_ideal(methane_params);

        auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        thermomodel::PengRobinson methane_pr(
            caloric, pr_params);

        for (double rpm : {0.0, 9000.0, 15000.0}) {
            run_case(air, "AirLikeIdeal", rpm);
            run_case(methane_ideal, "MethaneIdeal", rpm);
            run_case(methane_pr, "MethanePR", rpm);
        }

        std::cout
            << "PASS: passage-turning sensitivity tests\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
