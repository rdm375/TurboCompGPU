#include <turbocompgpu/compressor_map.hpp>
#include <turbocompgpu/diffusion_operating_point.hpp>

#include <thermomodel/caloric.hpp>
#include <thermomodel/peng_robinson.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

using namespace turbocompgpu;

namespace {

constexpr std::array<double, 5> speeds{
    9000.0, 10500.0, 12000.0, 13500.0, 15000.0
};

constexpr std::array<double, 5> initial_peaks{
    0.5171554,
    0.8171665,
    1.1404950,
    1.4940273,
    1.8842576
};

constexpr std::array<double, 5> tolerances{
    1e-4, 1e-5, 1e-6, 1e-7, 1e-8
};

constexpr std::array<double, 9> steps{
    1e-2, 5e-3, 2e-3,
    1e-3, 5e-4, 2e-4,
    1e-4, 5e-5, 1e-5
};

const ImpellerSystemGeometry geometry{
    {0.10, 0.20, 0.05},
    {0.40, 0.025, 0.10},
    {0.40, 45.0, 0.15}
};

const NonlinearResidualScales scales{
    50000.0, 100.0, 500000.0
};

const MapReferenceConditions reference{
    288.15, 101325.0
};

struct Value {
    double flow;
    double pressure_ratio;
    double residual;
};

class Characteristic {
public:
    Characteristic(
        const thermomodel::PengRobinson& model,
        double speed,
        double blade_angle)
        : model_(model),
          speed_(speed),
          blade_angle_(blade_angle)
    {}

    Value operator()(double flow) const
    {
        const NonlinearOperatingPoint operating{
            300.0, 500000.0,
            flow, speed_,
            0.90, 20.0, 20.0
        };

        const NonlinearUnknowns guess{
            340.0,
            std::log(700000.0 / scales.pressure)
        };

        const IncidenceLossParameters incidence{
            blade_angle_, 1.0
        };

        const DiffusionLossParameters diffusion{
            1.0, 1.5
        };

        NonlinearSolverOptions options;
        options.residual_tolerance = 1e-13;

        const auto result =
            evaluate_diffusion_operating_point(
                model_,
                geometry,
                {operating, guess},
                scales,
                incidence,
                diffusion,
                options
            );

        const auto& point = result.validation;

        if (!(result.independently_validated &&
              point.converged &&
              point.independently_validated &&
              point.basic_physics_valid)) {
            throw std::runtime_error(
                "Invalid operating point at flow=" +
                std::to_string(flow));
        }

        const auto map =
            evaluate_compressor_map_point(
                model_, geometry, point,
                scales, reference
            );

        if (!std::isfinite(map.total_pressure_ratio))
            throw std::runtime_error(
                "Nonfinite pressure ratio");

        return {
            flow,
            map.total_pressure_ratio,
            point.validation_residual_norm
        };
    }

private:
    const thermomodel::PengRobinson& model_;
    double speed_;
    double blade_angle_;
};

struct Maximum {
    Value value;
    double width;
    int iterations;
};

Maximum maximize(
    const Characteristic& characteristic,
    double left,
    double right,
    double tolerance)
{
    constexpr double golden =
        0.6180339887498948482;

    double a = left;
    double b = right;

    double c = b - golden * (b - a);
    double d = a + golden * (b - a);

    auto fc = characteristic(c);
    auto fd = characteristic(d);

    int iterations = 0;

    while (b - a > tolerance && iterations < 150) {
        if (fc.pressure_ratio < fd.pressure_ratio) {
            a = c;
            c = d;
            fc = fd;
            d = a + golden * (b - a);
            fd = characteristic(d);
        }
        else {
            b = d;
            d = c;
            fd = fc;
            c = b - golden * (b - a);
            fc = characteristic(c);
        }

        ++iterations;
    }

    if (b - a > tolerance)
        throw std::runtime_error(
            "Optimizer did not converge");

    return {
        characteristic(0.5 * (a + b)),
        b - a,
        iterations
    };
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const std::filesystem::path output =
            argc > 1
                ? argv[1]
                : "benchmarks/results-characteristic-convergence";

        std::filesystem::create_directories(output);

        std::ofstream optimizer(
            output / "optimizer_convergence.csv");

        std::ofstream derivatives(
            output / "derivative_convergence.csv");

        if (!optimizer || !derivatives)
            throw std::runtime_error(
                "Cannot open output files");

        optimizer << std::setprecision(17);
        derivatives << std::setprecision(17);

        optimizer
            << "rpm,tolerance,peak_flow,peak_ratio,"
            << "bracket_width,iterations,residual\n";

        derivatives
            << "rpm,peak_flow,h,first_derivative,"
            << "second_derivative,left_ratio,"
            << "center_ratio,right_ratio,"
            << "left_residual,center_residual,"
            << "right_residual\n";

        thermomodel::PRParameters parameters;
        parameters.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        const thermomodel::PengRobinson methane(
            std::make_shared<thermomodel::ConstantCp>(
                2200.0),
            parameters
        );

        const auto inlet =
            methane.evaluate_Tp(300.0, 500000.0);

        if (!inlet)
            throw std::runtime_error(
                "Inlet EOS evaluation failed");

        const OperatingIncidenceConditions design{
            5.0, 12000.0, inlet.value.rho, 20.0
        };

        const double blade_angle =
            calibrate_inlet_blade_angle_deg(
                geometry, design);

        std::cout << std::setprecision(12);

        for (std::size_t i = 0; i < speeds.size(); ++i) {
            const double speed = speeds[i];

            Characteristic characteristic(
                methane, speed, blade_angle);

            // Use the Stage 7G result only to establish
            // an initial bracket, not as a fixed solution.
            const double center = initial_peaks[i];

            const double left = center - 0.025;
            const double right = center + 0.025;

            const auto fleft = characteristic(left);
            const auto fcenter = characteristic(center);
            const auto fright = characteristic(right);

            if (!(fcenter.pressure_ratio >
                      fleft.pressure_ratio &&
                  fcenter.pressure_ratio >
                      fright.pressure_ratio)) {
                throw std::runtime_error(
                    "Initial maximum not bracketed at RPM=" +
                    std::to_string(speed));
            }

            Maximum refined{};

            for (double tolerance : tolerances) {
                const auto maximum = maximize(
                    characteristic,
                    left,
                    right,
                    tolerance
                );

                optimizer
                    << speed << ','
                    << tolerance << ','
                    << maximum.value.flow << ','
                    << maximum.value.pressure_ratio << ','
                    << maximum.width << ','
                    << maximum.iterations << ','
                    << maximum.value.residual << '\n';

                refined = maximum;
            }

            const double peak = refined.value.flow;
            const auto center_value = characteristic(peak);

            for (double h : steps) {
                const auto minus =
                    characteristic(peak - h);

                const auto plus =
                    characteristic(peak + h);

                const double first =
                    (plus.pressure_ratio -
                     minus.pressure_ratio) / (2.0 * h);

                const double second =
                    (plus.pressure_ratio -
                     2.0 * center_value.pressure_ratio +
                     minus.pressure_ratio) / (h * h);

                derivatives
                    << speed << ','
                    << peak << ','
                    << h << ','
                    << first << ','
                    << second << ','
                    << minus.pressure_ratio << ','
                    << center_value.pressure_ratio << ','
                    << plus.pressure_ratio << ','
                    << minus.residual << ','
                    << center_value.residual << ','
                    << plus.residual << '\n';
            }

            std::cout
                << "RPM=" << speed
                << " refined_peak=" << peak
                << " peak_ratio=" << refined.value.pressure_ratio
                << " bracket=" << refined.width
                << " residual=" << refined.value.residual
                << '\n';
        }

        std::cout
            << "\nSTAGE 7H SUMMARY\n"
            << "speeds=5\n"
            << "optimizer_tolerances=" << tolerances.size()
            << '\n'
            << "derivative_steps=" << steps.size()
            << '\n'
            << "output=" << output << '\n';

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
