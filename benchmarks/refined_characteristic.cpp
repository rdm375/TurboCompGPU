#include <turbocompgpu/compressor_map.hpp>
#include <turbocompgpu/diffusion_operating_point.hpp>
#include <turbocompgpu/passage_turning_choking.hpp>

#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

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
#include <vector>

namespace {

using namespace turbocompgpu;

constexpr std::array<double, 5> speeds{
    9000.0, 10500.0, 12000.0, 13500.0, 15000.0
};

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
    double flow;
    double ratio;
    double power;
    double incidence_loss;
    double diffusion_loss;
    double diffusion_ratio;
    double entropy;
    double residual;
};

struct Maximum {
    Evaluation point;
    double bracket_width;
    int iterations;
};

class Characteristic {
public:
    Characteristic(
        const thermomodel::PengRobinson& model,
        double speed,
        double blade_angle,
        std::ofstream& trace)
        : model_(model),
          speed_(speed),
          blade_angle_(blade_angle),
          trace_(trace)
    {}

    Evaluation evaluate(double flow)
    {
        if (!(std::isfinite(flow) && flow > 0.0))
            throw std::runtime_error("Invalid mass flow");

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
        options.residual_tolerance = 1.0e-13;

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
                "Invalid combined operating point at RPM=" +
                std::to_string(speed_) +
                " flow=" + std::to_string(flow));
        }

        const auto map =
            evaluate_compressor_map_point(
                model_, geometry, point,
                scales, reference
            );

        const Evaluation value{
            flow,
            map.total_pressure_ratio,
            map.shaft_power,
            result.incidence.loss.specific_loss,
            result.diffusion.loss.specific_loss,
            result.diffusion.loss.diffusion_ratio,
            result.total_entropy_generation,
            point.validation_residual_norm
        };

        trace_
            << speed_ << ','
            << value.flow << ','
            << value.ratio << ','
            << value.power << ','
            << value.incidence_loss << ','
            << value.diffusion_loss << ','
            << value.diffusion_ratio << ','
            << value.entropy << ','
            << value.residual << '\n';

        return value;
    }

private:
    const thermomodel::PengRobinson& model_;
    double speed_;
    double blade_angle_;
    std::ofstream& trace_;
};

// Golden-section maximization on a verified unimodal bracket.
// The calling code must first identify an interior sampled maximum.
Maximum maximize(
    Characteristic& characteristic,
    double left,
    double right,
    double flow_tolerance)
{
    constexpr double ratio =
        0.6180339887498948482;

    double a = left;
    double b = right;

    double c = b - ratio * (b - a);
    double d = a + ratio * (b - a);

    auto fc = characteristic.evaluate(c);
    auto fd = characteristic.evaluate(d);

    int iterations = 0;

    while (b - a > flow_tolerance &&
           iterations < 150) {
        if (fc.ratio < fd.ratio) {
            a = c;
            c = d;
            fc = fd;

            d = a + ratio * (b - a);
            fd = characteristic.evaluate(d);
        }
        else {
            b = d;
            d = c;
            fd = fc;

            c = b - ratio * (b - a);
            fc = characteristic.evaluate(c);
        }

        ++iterations;
    }

    if (b - a > flow_tolerance)
        throw std::runtime_error(
            "Maximum refinement did not converge");

    const auto peak =
        characteristic.evaluate(0.5 * (a + b));

    return {peak, b - a, iterations};
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const std::filesystem::path output =
            argc > 1
                ? argv[1]
                : "benchmarks/results-refined-characteristic";

        std::filesystem::create_directories(output);

        std::ofstream trace(output / "evaluations.csv");
        std::ofstream summary(output / "maxima.csv");

        if (!trace || !summary)
            throw std::runtime_error(
                "Cannot open benchmark output");

        trace << std::setprecision(17);
        summary << std::setprecision(17);

        trace
            << "rotational_speed_rpm,mass_flow,"
            << "total_pressure_ratio,shaft_power,"
            << "incidence_specific_loss,"
            << "diffusion_specific_loss,"
            << "diffusion_ratio,"
            << "total_entropy_generation,"
            << "validation_residual_norm\n";

        summary
            << "rotational_speed_rpm,"
            << "peak_found,"
            << "peak_mass_flow,"
            << "peak_total_pressure_ratio,"
            << "peak_shaft_power,"
            << "peak_diffusion_specific_loss,"
            << "peak_diffusion_ratio,"
            << "critical_mass_flow,"
            << "critical_total_pressure_ratio,"
            << "critical_shaft_power,"
            << "flow_interval,"
            << "relative_flow_interval,"
            << "peak_bracket_width,"
            << "peak_iterations,"
            << "curvature_estimate\n";

        thermomodel::PRParameters parameters;
        parameters.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        auto caloric =
            std::make_shared<thermomodel::ConstantCp>(
                2200.0);

        const thermomodel::PengRobinson methane(
            caloric, parameters);

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

        std::cout << std::setprecision(12)
                  << "blade_angle_deg="
                  << blade_angle << '\n';

        int peaks_found = 0;

        for (double speed : speeds) {
            Characteristic characteristic(
                methane, speed, blade_angle, trace);

            const PassageTurningChokingInputs input{
                300.0,
                500000.0,
                speed,
                20.0,
                throat,
                0.5
            };

            const auto choking =
                evaluate_passage_turning_choking(
                    methane, methane, geometry, input);

            const double critical =
                choking.choking.coupled.critical_mass_flow;

            const auto endpoint =
                characteristic.evaluate(critical);

            // Include flows below the Stage 7E lower bound.
            // Keep the scan within the provisional capacity.
            constexpr double minimum_flow = 0.1;
            constexpr double scan_step = 0.025;

            std::vector<Evaluation> samples;

            for (double flow = minimum_flow;
                 flow < critical;
                 flow += scan_step) {
                samples.push_back(
                    characteristic.evaluate(flow));
            }

            samples.push_back(endpoint);

            // Identify the highest sampled interior point.
            auto best = std::max_element(
                samples.begin(),
                samples.end(),
                [](const Evaluation& a,
                   const Evaluation& b) {
                    return a.ratio < b.ratio;
                });

            const auto index =
                static_cast<std::size_t>(
                    best - samples.begin());

            const bool interior =
                index > 0 &&
                index + 1 < samples.size();

            double peak_flow =
                std::numeric_limits<double>::quiet_NaN();

            double peak_ratio =
                std::numeric_limits<double>::quiet_NaN();

            double peak_power =
                std::numeric_limits<double>::quiet_NaN();

            double peak_loss =
                std::numeric_limits<double>::quiet_NaN();

            double peak_diffusion_ratio =
                std::numeric_limits<double>::quiet_NaN();

            double width =
                std::numeric_limits<double>::quiet_NaN();

            double curvature =
                std::numeric_limits<double>::quiet_NaN();

            int iterations = 0;

            if (interior) {
                const auto refined = maximize(
                    characteristic,
                    samples[index - 1].flow,
                    samples[index + 1].flow,
                    1.0e-7
                );

                const auto& peak = refined.point;

                peak_flow = peak.flow;
                peak_ratio = peak.ratio;
                peak_power = peak.power;
                peak_loss = peak.diffusion_loss;
                peak_diffusion_ratio =
                    peak.diffusion_ratio;

                width = refined.bracket_width;
                iterations = refined.iterations;

                // Symmetric curvature estimate. This is
                // a diagnostic, not part of the optimizer.
                constexpr double h = 0.002;

                const auto minus =
                    characteristic.evaluate(peak_flow - h);

                const auto plus =
                    characteristic.evaluate(peak_flow + h);

                curvature =
                    (plus.ratio -
                     2.0 * peak_ratio +
                     minus.ratio) / (h * h);

                if (!(curvature < 0.0))
                    throw std::runtime_error(
                        "Refined point is not a local maximum");

                ++peaks_found;
            }

            const double interval =
                critical - peak_flow;

            summary
                << speed << ','
                << interior << ','
                << peak_flow << ','
                << peak_ratio << ','
                << peak_power << ','
                << peak_loss << ','
                << peak_diffusion_ratio << ','
                << critical << ','
                << endpoint.ratio << ','
                << endpoint.power << ','
                << interval << ','
                << interval / critical << ','
                << width << ','
                << iterations << ','
                << curvature << '\n';

            std::cout
                << "RPM=" << speed
                << " peak_found=" << interior
                << " peak_flow=" << peak_flow
                << " critical_flow=" << critical
                << " interval=" << interval
                << " endpoint_Pi_tt=" << endpoint.ratio
                << " curvature=" << curvature
                << '\n';
        }

        std::cout
            << "\nSTAGE 7G SUMMARY\n"
            << "speeds=5\n"
            << "interior_peaks=" << peaks_found << '\n'
            << "output=" << output << '\n';

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
