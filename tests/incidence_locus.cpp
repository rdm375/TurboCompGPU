#include <turbocompgpu/incidence_locus.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <string>

using namespace turbocompgpu;

namespace {

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void check_close(
    double actual,
    double expected,
    double tolerance = 1e-11)
{
    const double scale =
        std::max({1.0, std::abs(actual), std::abs(expected)});

    require(
        std::isfinite(actual) &&
        std::abs(actual - expected) <= tolerance * scale,
        "Numerical mismatch"
    );
}

const ImpellerSystemGeometry geometry {
    {0.10, 0.20, 0.05},
    {0.40, 0.025, 0.10},
    {0.40, 45.0, 0.15}
};

constexpr std::array<double, 5> speeds {
    9000.0,
    10500.0,
    12000.0,
    13500.0,
    15000.0
};

void test_model(
    const thermomodel::Model& model,
    const char* name)
{
    const auto inlet = model.evaluate_Tp(
        300.0, 500000.0
    );

    require(
        static_cast<bool>(inlet),
        "Inlet EOS evaluation failed"
    );

    const double density = inlet.value.rho;

    const OperatingIncidenceConditions design {
        5.0,
        12000.0,
        density,
        20.0
    };

    const double blade_angle =
        calibrate_inlet_blade_angle_deg(
            geometry, design
        );

    const double D1 =
        geometry.impeller.inlet_diameter;

    const double U_design =
        std::numbers::pi * D1 *
        design.rotational_speed_rpm / 60.0;

    const double prewhirl =
        design.inlet_prewhirl;

    std::cout
        << "\n"
        << name
        << " density=" << density
        << " blade_angle=" << blade_angle
        << '\n';

    for (double speed : speeds) {
        const double U =
            std::numbers::pi * D1 *
            speed / 60.0;

        // Independent analytical prediction,
        // eliminating blade angle and density.
        const double predicted =
            design.mass_flow *
            (U - prewhirl) /
            (U_design - prewhirl);

        const double calculated =
            zero_incidence_mass_flow(
                geometry,
                density,
                speed,
                prewhirl,
                blade_angle
            );

        check_close(calculated, predicted);

        // Substitute the analytical mass flow
        // into the actual incidence evaluator.
        const OperatingIncidenceConditions operating {
            calculated,
            speed,
            density,
            prewhirl
        };

        const auto result =
            evaluate_operating_incidence_loss(
                geometry,
                operating,
                {blade_angle, 1.0}
            );

        check_close(
            result.loss.incidence_angle_deg,
            0.0
        );

        check_close(
            result.loss.specific_loss,
            0.0
        );

        // Verify that perturbations on either side
        // produce strictly positive incidence loss.
        for (double factor : {0.9, 1.1}) {
            auto perturbed = operating;
            perturbed.mass_flow *= factor;

            const auto off_design =
                evaluate_operating_incidence_loss(
                    geometry,
                    perturbed,
                    {blade_angle, 1.0}
                );

            require(
                off_design.loss.specific_loss > 0.0,
                "Perturbed incidence loss not positive"
            );
        }

        std::cout
            << "  speed=" << speed
            << " predicted_flow=" << predicted
            << " calculated_flow=" << calculated
            << " incidence_deg="
            << result.loss.incidence_angle_deg
            << '\n';
    }
}

} // namespace

int main()
{
    try {
        const thermomodel::IdealGas air_like;

        thermomodel::PRParameters pr_params;

        thermomodel::IdealGasParameters ideal_params;
        ideal_params.R =
            8.31446261815324 / pr_params.molar_mass;
        ideal_params.cp = 2200.0;
        ideal_params.Tref = pr_params.Tref;
        ideal_params.pref = pr_params.pref;

        const thermomodel::IdealGas methane_ideal(
            ideal_params
        );

        const auto caloric =
            std::make_shared<thermomodel::ConstantCp>(
                2200.0
            );

        const thermomodel::PengRobinson methane_pr(
            caloric, pr_params
        );

        test_model(air_like, "AirLikeIdeal");
        test_model(methane_ideal, "MethaneIdeal");
        test_model(methane_pr, "MethanePR");

        std::cout
            << "\nPASS: analytical zero-incidence locus\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
