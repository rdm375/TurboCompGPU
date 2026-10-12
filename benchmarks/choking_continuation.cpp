#include <turbocompgpu/choking_continuation.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace turbocompgpu;

int main(int argc, char** argv)
{
    try {
        const std::string filename =
            argc > 1 ? argv[1] : "choking_continuation.csv";

        const ImpellerSystemGeometry geometry{
            {0.10, 0.20, 0.05},
            {0.40, 0.025, 0.10},
            {0.40, 45.0, 0.15}
        };

        const ThroatGeometry throat{
            12, 0.15, 0.0035, 45.0, 0.002, 0.10
        };

        thermomodel::PRParameters pr_params;
        pr_params.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        thermomodel::PengRobinson methane_pr(
            caloric, pr_params);

        const PassageTurningChokingInputs baseline{
            300.0, 500000.0, 12000.0, 20.0, throat, 0.5
        };

        std::ofstream out(filename);

        if (!out)
            throw std::runtime_error("Cannot open CSV output");

        out << std::setprecision(17);

        out
            << "sweep,parameter,mdot_predicted,mdot_critical,"
            << "predictor_error,Wtau,Wtau_transport,"
            << "implicit_denominator,Cprime,gprime,"
            << "mass_residual,rothalpy_residual,"
            << "continuity_residual,sonic_residual,"
            << "throat_meridional_velocity,iterations,"
            << "near_singular\n";

        struct Sweep {
            const char* name;
            ChokingSweepParameter parameter;
            std::vector<double> values;
        };

        const std::vector<Sweep> sweeps{
            {
                "rpm",
                ChokingSweepParameter::rotational_speed_rpm,
                {0, 3000, 6000, 9000, 12000, 15000, 18000}
            },
            {
                "prewhirl",
                ChokingSweepParameter::inlet_prewhirl,
                {-40, -20, 0, 20, 40, 60}
            },
            {
                "blade_angle",
                ChokingSweepParameter::throat_blade_angle_deg,
                {25, 30, 35, 40, 45, 50, 55, 60, 65}
            },
            {
                "lambda",
                ChokingSweepParameter::turning_fraction,
                {0, 0.125, 0.25, 0.375, 0.5,
                 0.625, 0.75, 0.875, 1.0}
            }
        };

        for (const auto& sweep : sweeps) {
            const auto points =
                evaluate_choking_continuation(
                    methane_pr,
                    geometry,
                    baseline,
                    sweep.parameter,
                    sweep.values);

            for (const auto& point : points) {
                out
                    << sweep.name << ','
                    << point.parameter_value << ','
                    << point.predicted_mass_flow << ','
                    << point.critical_mass_flow << ','
                    << point.predictor_error << ','
                    << point.tangential_velocity << ','
                    << point.transported_tangential_velocity << ','
                    << point.implicit_denominator << ','
                    << point.capacity_derivative_w << ','
                    << point.transport_derivative_m << ','
                    << point.mass_flow_residual << ','
                    << point.rothalpy_residual << ','
                    << point.continuity_residual << ','
                    << point.sonic_residual << ','
                    << point.throat_meridional_velocity << ','
                    << point.iterations << ','
                    << (point.near_singular ? 1 : 0)
                    << '\n';

                std::cout
                    << sweep.name
                    << " x=" << point.parameter_value
                    << " mdot=" << point.critical_mass_flow
                    << " D=" << point.implicit_denominator
                    << " predictor_error="
                    << point.predictor_error
                    << '\n';
            }
        }

        std::cout << "Wrote " << filename << '\n';

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
