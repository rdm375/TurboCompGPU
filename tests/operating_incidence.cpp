#include <turbocompgpu/operating_incidence.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

using namespace turbocompgpu;

void check_close(double actual, double expected)
{
    const double scale =
        std::max({1.0, std::abs(actual), std::abs(expected)});

    if (!std::isfinite(actual) ||
        std::abs(actual - expected) > 1e-11 * scale) {
        throw std::runtime_error(
            "Operating-incidence mismatch"
        );
    }
}

void check_model(
    const thermomodel::Model& model,
    const char* name)
{
    const ImpellerSystemGeometry geometry {
        {0.10, 0.20, 0.05},
        {0.40, 0.025, 0.10},
        {0.40, 45.0, 0.15}
    };

    const auto state = model.evaluate_Tp(
        300.0, 500000.0
    );

    if (!state)
        throw std::runtime_error(
            "Inlet EOS evaluation failed"
        );

    const OperatingIncidenceConditions design {
        5.0,
        12000.0,
        state.value.rho,
        20.0
    };

    const double blade_angle =
        calibrate_inlet_blade_angle_deg(
            geometry, design
        );

    const IncidenceLossParameters params {
        blade_angle,
        1.0
    };

    const auto reference =
        evaluate_operating_incidence_loss(
            geometry, design, params
        );

    check_close(
        reference.loss.incidence_angle_deg,
        0.0
    );

    check_close(
        reference.loss.specific_loss,
        0.0
    );

    // Verify agreement with the existing impeller model.
    const auto impeller = evaluate_impeller(
        geometry.impeller,
        {
            design.rotational_speed_rpm,
            0.0,
            0.90,
            reference.inlet_meridional_velocity,
            design.inlet_prewhirl
        }
    );

    check_close(
        reference.inlet_blade_speed,
        impeller.inlet_blade_speed
    );

    check_close(
        reference.inlet_relative_velocity,
        impeller.inlet_relative_speed
    );

    check_close(
        reference.inlet_relative_angle_deg,
        impeller.inlet_relative_angle_deg
    );

    // Off-design losses should be positive.
    for (double flow : {2.0, 3.0, 4.0, 6.0, 7.0, 8.0}) {
        auto operating = design;
        operating.mass_flow = flow;

        const auto result =
            evaluate_operating_incidence_loss(
                geometry, operating, params
            );

        if (!(result.loss.specific_loss > 0.0))
            throw std::runtime_error(
                "Off-design incidence loss not positive"
            );
    }

    // Verify the analytic loss expression using the
    // velocity component perpendicular to the design
    // relative-flow direction.
    const double beta_design =
        blade_angle * std::numbers::pi / 180.0;

    const double Wtheta_design =
        design.inlet_prewhirl -
        reference.inlet_blade_speed;

    for (double flow : {2.0, 5.0, 8.0}) {
        auto operating = design;
        operating.mass_flow = flow;

        const auto result =
            evaluate_operating_incidence_loss(
                geometry, operating, params
            );

        const double transverse =
            result.inlet_meridional_velocity *
                std::cos(beta_design)
            +
            Wtheta_design *
                std::sin(beta_design);

        check_close(
            result.loss.specific_loss,
            0.5 * transverse * transverse
        );
    }

    // Check quadratic flow departure at fixed inlet
    // density, speed and prewhirl.
    auto low = design;
    low.mass_flow = 4.0;

    auto high = design;
    high.mass_flow = 6.0;

    const auto low_result =
        evaluate_operating_incidence_loss(
            geometry, low, params
        );

    const auto high_result =
        evaluate_operating_incidence_loss(
            geometry, high, params
        );

    check_close(
        low_result.loss.specific_loss,
        high_result.loss.specific_loss
    );

    std::cout
        << name
        << " density=" << state.value.rho
        << " blade_angle_deg=" << blade_angle
        << " design_relative_velocity="
        << reference.inlet_relative_velocity
        << " loss_at_4kg_s="
        << low_result.loss.specific_loss
        << '\n';
}

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

        check_model(air_like, "AirLikeIdeal");
        check_model(methane_ideal, "MethaneIdeal");
        check_model(methane_pr, "MethanePR");

        std::cout
            << "PASS: operating-incidence calibration\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';

        return 1;
    }
}
