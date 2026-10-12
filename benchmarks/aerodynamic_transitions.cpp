#include <turbocompgpu/compressor_map.hpp>
#include <turbocompgpu/diffusion_operating_point.hpp>
#include <turbocompgpu/passage_turning_choking.hpp>

#include <thermomodel/caloric.hpp>
#include <thermomodel/peng_robinson.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

using namespace turbocompgpu;

namespace {

constexpr std::array<double, 5> speeds{
    9000, 10500, 12000, 13500, 15000
};

constexpr std::array<double, 3> incidence_coefficients{
    0.5, 1.0, 1.5
};

constexpr std::array<double, 4> critical_ratios{
    1.3, 1.5, 1.7, 2.0
};

constexpr std::array<double, 5> diffusion_coefficients{
    0.25, 0.5, 1.0, 1.5, 2.0
};

constexpr std::array<double, 5> baseline_peaks{
    0.517155458,
    0.817166005,
    1.140494942,
    1.494026868,
    1.884257360
};

constexpr double minimum_flow = 0.1;
constexpr double derivative_step = 0.001;
constexpr double coefficient_step = 0.005;
constexpr double transition_tolerance = 1e-7;

const ImpellerSystemGeometry geometry{
    {0.10, 0.20, 0.05},
    {0.40, 0.025, 0.10},
    {0.40, 45.0, 0.15}
};

const ThroatGeometry throat{
    12, 0.15, 0.0035, 45.0, 0.002, 0.10
};

const NonlinearResidualScales scales{
    50000.0, 100.0, 500000.0
};

const MapReferenceConditions reference{
    288.15, 101325.0
};

struct Evaluation {
    double ratio;
    double residual;
};

class Characteristic {
public:
    Characteristic(
        const thermomodel::PengRobinson& model,
        double speed,
        double blade_angle,
        double ki,
        double dcrit)
        : model_(model),
          speed_(speed),
          blade_angle_(blade_angle),
          ki_(ki),
          dcrit_(dcrit)
    {}

    Evaluation evaluate(double flow, double kd) const
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
            blade_angle_, ki_
        };

        const DiffusionLossParameters diffusion{
            kd, dcrit_
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
                "Invalid operating point");
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
            map.total_pressure_ratio,
            point.validation_residual_norm
        };
    }

private:
    const thermomodel::PengRobinson& model_;
    double speed_;
    double blade_angle_;
    double ki_;
    double dcrit_;
};

struct Derivative {
    double value;
    double max_residual;
};

Derivative boundary_derivative(
    const Characteristic& characteristic,
    double boundary,
    double kd,
    bool upper,
    double h = derivative_step)
{
    const double direction = upper ? -1.0 : 1.0;

    const auto f0 =
        characteristic.evaluate(boundary, kd);

    const auto f1 =
        characteristic.evaluate(
            boundary + direction * h, kd);

    const auto f2 =
        characteristic.evaluate(
            boundary + direction * 2.0 * h, kd);

    const double derivative =
        upper
            ? (3.0 * f0.ratio -
               4.0 * f1.ratio +
               f2.ratio) / (2.0 * h)
            : (-3.0 * f0.ratio +
               4.0 * f1.ratio -
               f2.ratio) / (2.0 * h);

    return {
        derivative,
        std::max({
            f0.residual,
            f1.residual,
            f2.residual
        })
    };
}

struct Transition {
    double coefficient;
    double derivative;
    double bracket_width;
    int iterations;
};

Transition refine_transition(
    const Characteristic& characteristic,
    double boundary,
    bool upper,
    double left,
    double right)
{
    double fa = boundary_derivative(
        characteristic, boundary, left, upper).value;

    double fb = boundary_derivative(
        characteristic, boundary, right, upper).value;

    if (!std::isfinite(fa) ||
        !std::isfinite(fb) ||
        fa * fb > 0.0) {
        throw std::runtime_error(
            "Transition is not bracketed");
    }

    int iterations = 0;

    while (right - left > transition_tolerance &&
           iterations < 100) {
        const double mid = 0.5 * (left + right);

        const double fm = boundary_derivative(
            characteristic, boundary, mid, upper).value;

        if (!std::isfinite(fm))
            throw std::runtime_error(
                "Nonfinite boundary derivative");

        if (fm == 0.0) {
            left = mid;
            right = mid;
            break;
        }

        if (std::signbit(fa) == std::signbit(fm)) {
            left = mid;
            fa = fm;
        }
        else {
            right = mid;
            fb = fm;
        }

        ++iterations;
    }

    const double kd = 0.5 * (left + right);

    return {
        kd,
        boundary_derivative(
            characteristic, boundary, kd, upper).value,
        right - left,
        iterations
    };
}

struct Sensitivity {
    double first;
    double second;
    double mixed;
    double peak_sensitivity;
    double max_residual;
};

Sensitivity peak_sensitivity(
    const Characteristic& characteristic,
    double flow,
    double kd)
{
    const double h = derivative_step;
    const double dk = coefficient_step;

    const auto f0 = characteristic.evaluate(flow, kd);
    const auto fm = characteristic.evaluate(flow - h, kd);
    const auto fp = characteristic.evaluate(flow + h, kd);

    const double first =
        (fp.ratio - fm.ratio) / (2.0 * h);

    const double second =
        (fp.ratio - 2.0 * f0.ratio + fm.ratio) /
        (h * h);

    const auto pp =
        characteristic.evaluate(flow + h, kd + dk);

    const auto pm =
        characteristic.evaluate(flow + h, kd - dk);

    const auto mp =
        characteristic.evaluate(flow - h, kd + dk);

    const auto mm =
        characteristic.evaluate(flow - h, kd - dk);

    const double mixed =
        (pp.ratio - pm.ratio -
         mp.ratio + mm.ratio) /
        (4.0 * h * dk);

    if (!(std::isfinite(second) && second < 0.0))
        throw std::runtime_error(
            "Invalid peak curvature");

    const double residual = std::max({
        f0.residual,
        fm.residual,
        fp.residual,
        pp.residual,
        pm.residual,
        mp.residual,
        mm.residual
    });

    return {
        first,
        second,
        mixed,
        -mixed / second,
        residual
    };
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const std::filesystem::path output =
            argc > 1
                ? argv[1]
                : "benchmarks/results-aerodynamic-transitions";

        std::filesystem::create_directories(output);

        std::ofstream boundaries(
            output / "boundary_derivatives.csv");

        std::ofstream transitions(
            output / "transitions.csv");

        std::ofstream sensitivities(
            output / "peak_sensitivities.csv");

        if (!boundaries || !transitions || !sensitivities)
            throw std::runtime_error(
                "Cannot open benchmark outputs");

        boundaries << std::setprecision(17);
        transitions << std::setprecision(17);
        sensitivities << std::setprecision(17);

        boundaries
            << "rpm,k_inc,d_crit,k_diff,boundary,"
            << "mass_flow,derivative,max_residual\n";

        transitions
            << "rpm,k_inc,d_crit,boundary,"
            << "k_diff_left,k_diff_right,"
            << "transition_coefficient,"
            << "derivative,bracket_width,iterations\n";

        sensitivities
            << "rpm,k_inc,d_crit,k_diff,peak_flow,"
            << "first_derivative,curvature,mixed_derivative,"
            << "dpeak_dkdiff,max_residual\n";

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

        int families = 0;
        int crossings = 0;

        for (std::size_t s = 0; s < speeds.size(); ++s) {
            const double speed = speeds[s];

            const PassageTurningChokingInputs input{
                300.0, 500000.0,
                speed, 20.0,
                throat, 0.5
            };

            const auto choking =
                evaluate_passage_turning_choking(
                    methane, methane, geometry, input);

            const double critical =
                choking.choking.coupled.critical_mass_flow;

            if (!(critical > minimum_flow +
                    2.0 * derivative_step)) {
                throw std::runtime_error(
                    "Insufficient flow interval");
            }

            for (double ki : incidence_coefficients) {
                for (double dcrit : critical_ratios) {
                    Characteristic characteristic(
                        methane, speed, blade_angle,
                        ki, dcrit
                    );

                    ++families;

                    for (int side = 0; side < 2; ++side) {
                        const bool upper = side == 1;

                        const double boundary =
                            upper ? critical : minimum_flow;

                        const char* label =
                            upper ? "capacity" : "low";

                        std::array<double, 5> slopes{};

                        for (std::size_t j = 0;
                             j < diffusion_coefficients.size();
                             ++j) {
                            const double kd =
                                diffusion_coefficients[j];

                            const auto derivative =
                                boundary_derivative(
                                    characteristic,
                                    boundary,
                                    kd,
                                    upper
                                );

                            slopes[j] = derivative.value;

                            boundaries
                                << speed << ','
                                << ki << ','
                                << dcrit << ','
                                << kd << ','
                                << label << ','
                                << boundary << ','
                                << derivative.value << ','
                                << derivative.max_residual
                                << '\n';
                        }

                        for (std::size_t j = 0;
                             j + 1 < slopes.size();
                             ++j) {
                            if (slopes[j] * slopes[j + 1] >= 0.0)
                                continue;

                            const double left =
                                diffusion_coefficients[j];

                            const double right =
                                diffusion_coefficients[j + 1];

                            const auto root =
                                refine_transition(
                                    characteristic,
                                    boundary,
                                    upper,
                                    left,
                                    right
                                );

                            transitions
                                << speed << ','
                                << ki << ','
                                << dcrit << ','
                                << label << ','
                                << left << ','
                                << right << ','
                                << root.coefficient << ','
                                << root.derivative << ','
                                << root.bracket_width << ','
                                << root.iterations
                                << '\n';

                            ++crossings;
                        }
                    }

                    // Baseline coefficient sensitivity:
                    // Kinc=1, Dcrit=1.5, Kdiff=1.
                    if (ki == 1.0 && dcrit == 1.5) {
                        const double flow = baseline_peaks[s];

                        const auto sensitivity =
                            peak_sensitivity(
                                characteristic, flow, 1.0
                            );

                        sensitivities
                            << speed << ','
                            << ki << ','
                            << dcrit << ','
                            << 1.0 << ','
                            << flow << ','
                            << sensitivity.first << ','
                            << sensitivity.second << ','
                            << sensitivity.mixed << ','
                            << sensitivity.peak_sensitivity << ','
                            << sensitivity.max_residual
                            << '\n';
                    }
                }
            }

            boundaries.flush();
            transitions.flush();
            sensitivities.flush();

            std::cout
                << "RPM=" << speed
                << " families=" << families
                << " crossings=" << crossings
                << '\n';
        }

        std::cout
            << "\nSTAGE 7I.2 SUMMARY\n"
            << "families=" << families << '\n'
            << "crossings=" << crossings << '\n'
            << "output=" << output << '\n';

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr
            << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
