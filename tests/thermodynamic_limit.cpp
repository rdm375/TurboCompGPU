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

constexpr double Ru = 8.31446261815324;
constexpr double cp_ideal = 2200.0;
constexpr double temperature = 300.0;

struct Difference {
    double density_relative;
    double compressibility;
    double cp_relative;
};

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

Difference compare(
    const thermomodel::Model& ideal,
    const thermomodel::Model& pr,
    double pressure)
{
    const auto a = ideal.evaluate_Tp(temperature, pressure);
    const auto b = pr.evaluate_Tp(temperature, pressure);

    require(
        static_cast<bool>(a) && static_cast<bool>(b),
        "Thermodynamic evaluation failed at p=" +
            std::to_string(pressure)
    );

    const auto& ig = a.value;
    const auto& real = b.value;

    require(
        std::isfinite(ig.rho) &&
        std::isfinite(real.rho) &&
        ig.rho > 0.0 &&
        real.rho > 0.0,
        "Invalid density"
    );

    require(
        std::isfinite(real.Z) &&
        std::isfinite(real.cp),
        "Invalid real-gas properties"
    );

    return {
        std::abs(real.rho / ig.rho - 1.0),
        std::abs(real.Z - 1.0),
        std::abs(real.cp / ig.cp - 1.0)
    };
}

} // namespace

int main()
{
    try {
        thermomodel::PRParameters pr_params;

        thermomodel::IdealGasParameters ideal_params;
        ideal_params.R = Ru / pr_params.molar_mass;
        ideal_params.cp = cp_ideal;
        ideal_params.Tref = pr_params.Tref;
        ideal_params.pref = pr_params.pref;

        const thermomodel::IdealGas ideal(ideal_params);

        const auto caloric =
            std::make_shared<thermomodel::ConstantCp>(
                cp_ideal
            );

        const thermomodel::PengRobinson pr(
            caloric, pr_params
        );

        // Stay within the models' documented pressure ranges.
        constexpr std::array<double, 6> pressures {
            1000.0,
            10000.0,
            100000.0,
            500000.0,
            1000000.0,
            5000000.0
        };

        std::cout << std::setprecision(12);

        std::cout
            << "pressure_Pa,"
            << "density_relative_error,"
            << "Z_error,"
            << "cp_relative_error\n";

        std::array<Difference, pressures.size()> errors {};

        for (std::size_t i = 0; i < pressures.size(); ++i) {
            errors[i] = compare(
                ideal, pr, pressures[i]
            );

            std::cout
                << pressures[i] << ','
                << errors[i].density_relative << ','
                << errors[i].compressibility << ','
                << errors[i].cp_relative
                << '\n';
        }

        // Check approach toward the ideal-gas limit.
        // Do not require exact equality at finite pressure.
        const auto& low = errors.front();
        const auto& high = errors.back();

        require(
            low.density_relative <
                high.density_relative,
            "Density does not approach ideal behavior"
        );

        require(
            low.compressibility <
                high.compressibility,
            "Compressibility does not approach unity"
        );

        require(
            low.cp_relative <
                high.cp_relative,
            "Heat capacity does not approach ideal behavior"
        );

        // Low-pressure deviations should be small.
        require(
            low.density_relative < 1.0e-3,
            "Low-pressure density deviation too large"
        );

        require(
            low.compressibility < 1.0e-3,
            "Low-pressure Z deviation too large"
        );

        require(
            low.cp_relative < 1.0e-3,
            "Low-pressure cp deviation too large"
        );

        // Compare reference-invariant property changes.
        // These are differences along an isothermal path,
        // so constant enthalpy/entropy offsets cancel.
        constexpr double p0 = 1000.0;
        constexpr double p1 = 10000.0;

        const auto ig0 = ideal.evaluate_Tp(temperature, p0);
        const auto ig1 = ideal.evaluate_Tp(temperature, p1);
        const auto pr0 = pr.evaluate_Tp(temperature, p0);
        const auto pr1 = pr.evaluate_Tp(temperature, p1);

        require(
            ig0 && ig1 && pr0 && pr1,
            "Reference-invariant comparison failed"
        );

        const double delta_h_ideal =
            ig1.value.h - ig0.value.h;

        const double delta_h_pr =
            pr1.value.h - pr0.value.h;

        const double delta_s_ideal =
            ig1.value.s - ig0.value.s;

        const double delta_s_pr =
            pr1.value.s - pr0.value.s;

        const double h_difference =
            delta_h_pr - delta_h_ideal;

        const double s_difference =
            delta_s_pr - delta_s_ideal;

        require(
            std::isfinite(h_difference) &&
            std::isfinite(s_difference),
            "Nonfinite reference-invariant differences"
        );

        std::cout
            << "\nReference-invariant differences:\n"
            << "delta_h_PR_minus_ideal = "
            << h_difference << " J/kg\n"
            << "delta_s_PR_minus_ideal = "
            << s_difference << " J/(kg K)\n";

        std::cout
            << "\nPASS: methane ideal-gas limiting behavior\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
