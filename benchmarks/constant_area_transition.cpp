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
        double dcrit,
        const ImpellerSystemGeometry& geometry)
        : model_(model),
          geometry_(geometry),
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
                geometry_,
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
                model_, geometry_, point,
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
    const ImpellerSystemGeometry& geometry_;
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



namespace {

constexpr double speed_j = 15000.0;
constexpr double scale_left = 0.95;
constexpr double scale_right = 1.00;
constexpr double scale_tolerance = 1e-9;

constexpr std::array<double, 5> steps{
    0.004, 0.002, 0.001, 0.0005, 0.00025
};

struct ScaleEvaluation {
    double slope;
    double max_residual;
};

ScaleEvaluation evaluate_scale(
    const thermomodel::PengRobinson& methane,
    double blade_angle,
    double critical,
    double scale,
    double h)
{
    ImpellerSystemGeometry g = geometry;

    g.exit.diameter *= scale;
    g.impeller.exit_diameter *= scale;
    g.exit.width /= scale;

    const double area0 =
        effective_flow_area(
            impeller_exit_area(geometry.exit),
            geometry.exit.blockage);

    const double area =
        effective_flow_area(
            impeller_exit_area(g.exit),
            g.exit.blockage);

    if (std::abs(area / area0 - 1.0) > 1e-12)
        throw std::runtime_error(
            "Constant-area invariant violated");

    Characteristic characteristic(
        methane,
        speed_j,
        blade_angle,
        1.0,
        1.5,
        g
    );

    const auto result = boundary_derivative(
        characteristic,
        critical,
        1.0,
        true,
        h
    );

    return {
        result.value,
        result.max_residual
    };
}

struct ScaleRoot {
    double scale;
    double slope;
    double max_residual;
    double bracket_width;
    int iterations;
    double left_slope;
    double right_slope;
};

ScaleRoot solve_scale(
    const thermomodel::PengRobinson& methane,
    double blade_angle,
    double critical,
    double h)
{
    double left = scale_left;
    double right = scale_right;

    const auto initial_left =
        evaluate_scale(
            methane, blade_angle, critical, left, h);

    const auto initial_right =
        evaluate_scale(
            methane, blade_angle, critical, right, h);

    double fa = initial_left.slope;
    double fb = initial_right.slope;

    if (!std::isfinite(fa) ||
        !std::isfinite(fb) ||
        std::signbit(fa) == std::signbit(fb)) {
        throw std::runtime_error(
            "Constant-area transition is not bracketed");
    }

    int iterations = 0;

    while (right - left > scale_tolerance &&
           iterations < 100) {
        const double mid = 0.5 * (left + right);

        const double fm =
            evaluate_scale(
                methane, blade_angle, critical, mid, h
            ).slope;

        if (!std::isfinite(fm))
            throw std::runtime_error(
                "Nonfinite transition derivative");

        if (fm == 0.0) {
            left = mid;
            right = mid;
            break;
        }

        if (std::signbit(fm) == std::signbit(fa)) {
            left = mid;
            fa = fm;
        }
        else {
            right = mid;
            fb = fm;
        }

        ++iterations;
    }

    const double scale = 0.5 * (left + right);

    const auto final =
        evaluate_scale(
            methane, blade_angle, critical, scale, h);

    return {
        scale,
        final.slope,
        final.max_residual,
        right - left,
        iterations,
        initial_left.slope,
        initial_right.slope
    };
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const std::filesystem::path output =
            argc > 1
                ? argv[1]
                : "benchmarks/results-constant-area-transition";

        std::filesystem::create_directories(output);

        std::ofstream csv(
            output / "constant_area_transition.csv");

        if (!csv)
            throw std::runtime_error(
                "Cannot open transition CSV");

        csv << std::setprecision(17);

        csv
            << "rpm,h,critical_flow,"
            << "scale_left,scale_right,"
            << "scale_root,radius_root,width_root,"
            << "derivative,max_residual,"
            << "bracket_width,iterations,"
            << "left_derivative,right_derivative\n";

        thermomodel::PRParameters parameters;
        parameters.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        const thermomodel::PengRobinson methane(
            std::make_shared<thermomodel::ConstantCp>(2200.0),
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

        const PassageTurningChokingInputs input{
            300.0,
            500000.0,
            speed_j,
            20.0,
            throat,
            0.5
        };

        const auto choking =
            evaluate_passage_turning_choking(
                methane, methane, geometry, input);

        const double critical =
            choking.choking.coupled.critical_mass_flow;

        std::cout << std::setprecision(12);

        std::cout
            << "STAGE 7J.3 CONSTANT-AREA TRANSITION\n"
            << "RPM=" << speed_j << '\n'
            << "critical_flow=" << critical << '\n'
            << "baseline_radius="
            << 0.5 * geometry.exit.diameter << '\n'
            << "baseline_width="
            << geometry.exit.width << '\n';

        for (double h : steps) {
            const auto root =
                solve_scale(
                    methane, blade_angle, critical, h);

            const double radius =
                0.5 * geometry.exit.diameter * root.scale;

            const double width =
                geometry.exit.width / root.scale;

            csv
                << speed_j << ','
                << h << ','
                << critical << ','
                << scale_left << ','
                << scale_right << ','
                << root.scale << ','
                << radius << ','
                << width << ','
                << root.slope << ','
                << root.max_residual << ','
                << root.bracket_width << ','
                << root.iterations << ','
                << root.left_slope << ','
                << root.right_slope
                << '\n';

            std::cout
                << "h=" << h
                << " scale=" << root.scale
                << " radius_mm=" << 1000.0 * radius
                << " width_mm=" << 1000.0 * width
                << " slope=" << root.slope
                << " bracket=" << root.bracket_width
                << " residual=" << root.max_residual
                << '\n';
        }

        std::cout
            << "output=" << output << '\n';

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr
            << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
