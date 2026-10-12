#include <turbocompgpu/complete_stage_solver.hpp>
#include <turbocompgpu/station.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <string>

using namespace turbocompgpu;

namespace {

constexpr std::array<double, 7> flows{
    2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0
};

constexpr std::array<double, 5> speeds{
    9000.0, 10500.0, 12000.0, 13500.0, 15000.0
};

ImpellerSystemGeometry geometry()
{
    return {
        {0.10, 0.20, 0.05},
        {0.40, 0.025, 0.10},
        {0.40, 45.0, 0.15}
    };
}

NonlinearResidualScales scales()
{
    return {50000.0, 100.0, 500000.0};
}

OperatingPointRequest request(double flow, double rpm)
{
    return {
        {300.0, 500000.0, flow, rpm, 0.90, 20.0, 20.0},
        {340.0, std::log(700000.0 / 500000.0)}
    };
}

NonlinearSolverOptions options()
{
    NonlinearSolverOptions result;
    result.residual_tolerance = 1.0e-13;
    return result;
}

const char* stage_status(CompleteStageStatus status)
{
    switch (status) {
    case CompleteStageStatus::Complete:
        return "Complete";
    case CompleteStageStatus::ImpellerFailed:
        return "ImpellerFailed";
    case CompleteStageStatus::DiffuserFailed:
        return "DiffuserFailed";
    }
    return "Unknown";
}

const char* diffuser_status(DiffuserBranchStatus status)
{
    switch (status) {
    case DiffuserBranchStatus::Converged:
        return "Converged";
    case DiffuserBranchStatus::Critical:
        return "Critical";
    case DiffuserBranchStatus::AboveCriticalCapacity:
        return "AboveCriticalCapacity";
    case DiffuserBranchStatus::BranchUnavailable:
        return "BranchUnavailable";
    case DiffuserBranchStatus::ThermodynamicDomainLimit:
        return "ThermodynamicDomainLimit";
    }
    return "Unknown";
}

template<class Model>
VanelessDiffuserGeometry calibrate(
    const Model& model,
    const ImpellerSystemGeometry& g,
    const NonlinearResidualScales& s)
{
    const auto op = request(5.0, 12000.0);

    const auto impeller = evaluate_operating_point(
        model, g, op, s, options());

    if (!(impeller.converged &&
          impeller.independently_validated &&
          impeller.basic_physics_valid)) {
        throw std::runtime_error(
            "Reference impeller calibration failed");
    }

    const auto state2 = evaluate(
        model,
        impeller.outlet_temperature,
        impeller.outlet_pressure);

    const double mass = op.operating.mass_flow;
    const double r2 = g.exit.diameter / 2.0;
    const double r3 = 1.5 * r2;

    const double area2 = effective_flow_area(
        impeller_exit_area(g.exit), g.exit.blockage);

    const double vm2 = mass / (state2.rho * area2);

    const double omega =
        2.0 * std::numbers::pi *
        op.operating.rotational_speed_rpm / 60.0;

    const double beta =
        g.impeller.exit_blade_angle_deg *
        std::numbers::pi / 180.0;

    const double whirl2 =
        op.operating.slip_factor *
        (omega * r2 - vm2 / std::tan(beta));

    const double whirl3 = r2 * whirl2 / r3;

    const FlowStation station2{
        state2, {vm2, whirl2}
    };

    const auto stagnation2 =
        evaluate_stagnation_state(model, station2);

    const double p_reference =
        state2.p + 0.10 *
        (stagnation2.stagnation_state.p - state2.p);

    const auto reference = invert_temperature(
        model,
        p_reference,
        state2.s,
        Property::Entropy,
        state2.T);

    const double vm3_squared =
        2.0 * (
            stagnation2.stagnation_enthalpy -
            reference.state.h
        ) - whirl3 * whirl3;

    if (!(vm3_squared > 0.0)) {
        throw std::runtime_error(
            "Invalid diffuser calibration velocity");
    }

    const double area3 =
        mass / (reference.state.rho * std::sqrt(vm3_squared));

    const double width3 =
        area3 / (2.0 * std::numbers::pi * r3);

    return {r2, r3, width3, 0.0};
}

struct Counts {
    int complete = 0;
    int impeller_failed = 0;
    int diffuser_failed = 0;
    int capacity_limited = 0;
    int domain_limited = 0;
    double maximum_energy_residual = 0.0;
};

template<class Model>
void sweep(
    const Model& model,
    const std::string& name,
    const std::filesystem::path& output)
{
    const auto g = geometry();
    const auto s = scales();
    const auto diffuser_geometry = calibrate(model, g, s);

    std::filesystem::create_directories(output.parent_path());

    std::ofstream csv(output);

    if (!csv)
        throw std::runtime_error("Cannot open " + output.string());

    csv << std::setprecision(17);

    csv
        << "model,rpm,mass_flow,stage_status,diffuser_status,"
        << "impeller_converged,impeller_validated,physics_valid,"
        << "impeller_iterations,impeller_residual,"
        << "p2_pa,p3_pa,p01_pa,p03_pa,"
        << "total_pressure_ratio,static_pressure_recovery,"
        << "specific_work_j_kg,shaft_power_w,"
        << "diffuser_critical_mass_flow_kg_s,"
        << "diffuser_capacity_margin,meridional_mach3,"
        << "combined_energy_residual_j_kg,"
        << "mass_residual_kg_s,"
        << "entropy_residual_j_kg_k,"
        << "angular_momentum_residual_m2_s\n";

    Counts counts;
    constexpr double nan =
        std::numeric_limits<double>::quiet_NaN();

    for (double rpm : speeds) {
        for (double mass : flows) {
            const auto op = request(mass, rpm);

            const auto result = evaluate_complete_stage(
                model, model, g, diffuser_geometry,
                op, s, options());

            if (result.complete())
                ++counts.complete;
            else if (result.status ==
                     CompleteStageStatus::ImpellerFailed)
                ++counts.impeller_failed;
            else
                ++counts.diffuser_failed;

            double p2 = nan;
            double p3 = nan;
            double p01 = nan;
            double p03 = nan;
            double ratio = nan;
            double recovery = nan;
            double critical_mass = nan;
            double margin = nan;
            double mach3 = nan;
            double energy = nan;
            double mass_residual = nan;
            double entropy = nan;
            double angular_momentum = nan;

            std::string diffuser_class = "NotAttempted";

            if (result.impeller_valid()) {
                p2 = result.impeller.outlet_pressure;

                if (result.stage) {
                    const auto& stage = *result.stage;
                    diffuser_class =
                        diffuser_status(stage.status);

                    if (std::isfinite(
                            stage.diffuser_critical_mass_flow) &&
                        stage.diffuser_critical_mass_flow > 0.0) {
                        critical_mass =
                            stage.diffuser_critical_mass_flow;

                        margin = 1.0 - mass / critical_mass;
                    }

                    if (stage.status ==
                        DiffuserBranchStatus::AboveCriticalCapacity)
                        ++counts.capacity_limited;

                    if (stage.status ==
                        DiffuserBranchStatus::ThermodynamicDomainLimit)
                        ++counts.domain_limited;
                }
            }

            if (result.complete()) {
                const auto& stage = *result.stage;
                const auto& station3 = *stage.station3;

                p3 = station3.pressure;
                recovery = p3 / p2;
                mach3 = station3.meridional_mach;

                energy = stage.combined_energy_residual;
                mass_residual = stage.mass_residual;
                entropy = stage.diffuser_entropy_residual;
                angular_momentum =
                    stage.diffuser_angular_momentum_residual;

                counts.maximum_energy_residual =
                    std::max(
                        counts.maximum_energy_residual,
                        std::abs(energy));

                // Stagnation-pressure inversion is deliberately
                // performed only for completed stage solutions.
                try {
                    const auto state1 = evaluate(
                        model,
                        op.operating.inlet_temperature,
                        op.operating.inlet_pressure);

                    const double area1 = effective_flow_area(
                        inlet_annulus_area(g.inlet),
                        g.inlet.blockage);

                    const double vm1 =
                        mass / (state1.rho * area1);

                    const FlowStation station1{
                        state1,
                        {vm1, op.operating.inlet_prewhirl}
                    };

                    const FlowStation station3_flow{
                        station3.state,
                        {
                            station3.meridional_velocity,
                            station3.tangential_velocity
                        }
                    };

                    p01 = evaluate_stagnation_state(
                        model, station1).stagnation_state.p;

                    p03 = evaluate_stagnation_state(
                        model, station3_flow).stagnation_state.p;

                    ratio = p03 / p01;
                }
                catch (const std::exception&) {
                    // A stagnation-property failure is distinct
                    // from the completed stage solution.
                    // Leave derived pressures as NaN.
                }
            }

            csv
                << name << ','
                << rpm << ','
                << mass << ','
                << stage_status(result.status) << ','
                << diffuser_class << ','
                << result.impeller.converged << ','
                << result.impeller.independently_validated << ','
                << result.impeller.basic_physics_valid << ','
                << result.impeller.iterations << ','
                << result.impeller.residual_norm << ','
                << p2 << ','
                << p3 << ','
                << p01 << ','
                << p03 << ','
                << ratio << ','
                << recovery << ','
                << result.impeller.specific_work << ','
                << result.impeller.shaft_power << ','
                << critical_mass << ','
                << margin << ','
                << mach3 << ','
                << energy << ','
                << mass_residual << ','
                << entropy << ','
                << angular_momentum << '\n';
        }
    }

    if (!csv)
        throw std::runtime_error("CSV write failure");

    std::cout
        << name
        << " width3=" << diffuser_geometry.outlet_width
        << " complete=" << counts.complete
        << " impeller_failed=" << counts.impeller_failed
        << " diffuser_failed=" << counts.diffuser_failed
        << " capacity_limited=" << counts.capacity_limited
        << " domain_limited=" << counts.domain_limited
        << " max_energy_residual="
        << counts.maximum_energy_residual
        << " output=" << output.string()
        << '\n';
}

} // namespace

int main()
{
    try {
        const std::filesystem::path directory =
            "benchmarks/results-complete-stage-characteristic";

        const thermomodel::IdealGas ideal;

        sweep(
            ideal,
            "IdealGas",
            directory / "ideal-gas.csv");

        const auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        thermomodel::PRParameters parameters;
        parameters.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        const thermomodel::PengRobinson pr(
            caloric, parameters);

        sweep(
            pr,
            "PengRobinson",
            directory / "peng-robinson.csv");

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
