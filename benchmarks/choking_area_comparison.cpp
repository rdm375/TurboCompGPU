#include <turbocompgpu/choking_area_comparison.hpp>

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
            argc > 1
                ? argv[1]
                : "choking_area_comparison.csv";

        const ImpellerSystemGeometry geometry{
            {0.10, 0.20, 0.05},
            {0.40, 0.025, 0.10},
            {0.40, 45.0, 0.15}
        };

        const ThroatGeometry throat{
            12, 0.15, 0.0035, 45.0, 0.002, 0.10
        };

        const PassageTurningChokingInputs baseline{
            300.0, 500000.0, 12000.0, 20.0, throat, 0.5
        };

        thermomodel::PRParameters parameters;
        parameters.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        auto caloric =
            std::make_shared<thermomodel::ConstantCp>(
                2200.0);

        thermomodel::PengRobinson methane(
            caloric, parameters);

        const std::vector<double> angles{
            25, 30, 35, 40, 45, 50, 55, 60, 65
        };

        const auto points =
            evaluate_choking_area_comparison(
                methane, geometry, baseline, angles);

        std::ofstream out(filename);

        if (!out)
            throw std::runtime_error("Cannot open CSV output");

        out << std::setprecision(17);

        out
            << "blade_angle_deg,"
            << "physical_height,constant_area_height,"
            << "physical_area,constant_area,"
            << "physical_mdot,constant_area_mdot,"
            << "physical_mass_flux,constant_area_mass_flux,"
            << "physical_Wtau,constant_area_Wtau,"
            << "physical_mass_residual,constant_area_mass_residual,"
            << "physical_rothalpy_residual,"
            << "constant_area_rothalpy_residual,"
            << "mdot_difference,mdot_ratio\n";

        for (const auto& p : points) {
            out
                << p.blade_angle_deg << ','
                << p.physical_height << ','
                << p.constant_area_height << ','
                << p.physical_area << ','
                << p.constant_area << ','
                << p.physical_mass_flow << ','
                << p.constant_area_mass_flow << ','
                << p.physical_mass_flux << ','
                << p.constant_area_mass_flux << ','
                << p.physical_tangential_velocity << ','
                << p.constant_area_tangential_velocity << ','
                << p.physical_mass_residual << ','
                << p.constant_area_mass_residual << ','
                << p.physical_rothalpy_residual << ','
                << p.constant_area_rothalpy_residual << ','
                << p.mass_flow_difference << ','
                << p.mass_flow_ratio
                << '\n';

            std::cout
                << std::setprecision(12)
                << "beta=" << p.blade_angle_deg
                << " Aphysical=" << p.physical_area
                << " Aconstant=" << p.constant_area
                << " bconstant=" << p.constant_area_height
                << " mdot_physical=" << p.physical_mass_flow
                << " mdot_constant=" << p.constant_area_mass_flow
                << " ratio=" << p.mass_flow_ratio
                << '\n';
        }

        std::cout << "Wrote " << filename << '\n';
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
