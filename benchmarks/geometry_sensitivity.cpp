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
#include <vector>

using namespace turbocompgpu;

namespace {

constexpr std::array<double, 5> speeds{
    9000, 10500, 12000, 13500, 15000
};

constexpr std::array<double, 3> incidence_coefficients{
    0.5, 1.0, 1.5
};

constexpr std::array<double, 5> diffusion_coefficients{
    0.25, 0.5, 1.0, 1.5, 2.0
};

constexpr std::array<double, 4> critical_ratios{
    1.3, 1.5, 1.7, 2.0
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

constexpr double minimum_flow = 0.1;
constexpr double scan_step = 0.025;
constexpr double peak_tolerance = 1e-6;

struct Evaluation {
    double flow;
    double pressure_ratio;
    double incidence_loss;
    double diffusion_loss;
    double diffusion_ratio;
    double residual;
};

class Characteristic {
public:
    Characteristic(
        const thermomodel::PengRobinson& model,
        double speed,
        double blade_angle,
        double ki,
        double kd,
        double dcrit,
        const ImpellerSystemGeometry& geometry)
        : model_(model),
          geometry_(geometry),
          speed_(speed),
          blade_angle_(blade_angle),
          ki_(ki),
          kd_(kd),
          dcrit_(dcrit)
    {}

    Evaluation evaluate(double flow) const
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
            kd_, dcrit_
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
                "Operating-point validation failed");
        }

        const auto map =
            evaluate_compressor_map_point(
                model_, geometry_, point,
                scales, reference
            );

        return {
            flow,
            map.total_pressure_ratio,
            result.incidence.loss.specific_loss,
            result.diffusion.loss.specific_loss,
            result.diffusion.loss.diffusion_ratio,
            point.validation_residual_norm
        };
    }

private:
    const thermomodel::PengRobinson& model_;
    const ImpellerSystemGeometry& geometry_;
    double speed_;
    double blade_angle_;
    double ki_;
    double kd_;
    double dcrit_;
};

Evaluation maximize(
    const Characteristic& characteristic,
    double left,
    double right)
{
    constexpr double golden =
        0.6180339887498948482;

    double a = left;
    double b = right;

    double c = b - golden * (b - a);
    double d = a + golden * (b - a);

    auto fc = characteristic.evaluate(c);
    auto fd = characteristic.evaluate(d);

    for (int iteration = 0; iteration < 150; ++iteration) {
        if (b - a <= peak_tolerance)
            return characteristic.evaluate(
                0.5 * (a + b));

        if (fc.pressure_ratio < fd.pressure_ratio) {
            a = c;
            c = d;
            fc = fd;

            d = a + golden * (b - a);
            fd = characteristic.evaluate(d);
        }
        else {
            b = d;
            d = c;
            fd = fc;

            c = b - golden * (b - a);
            fc = characteristic.evaluate(c);
        }
    }

    throw std::runtime_error(
        "Peak refinement failed");
}

struct Result {
    bool peak_found = false;
    int valid_points = 0;
    int failed_points = 0;
    int candidate_peaks = 0;

    double peak_flow =
        std::numeric_limits<double>::quiet_NaN();

    double peak_ratio =
        std::numeric_limits<double>::quiet_NaN();

    double peak_diffusion_loss =
        std::numeric_limits<double>::quiet_NaN();

    double peak_diffusion_ratio =
        std::numeric_limits<double>::quiet_NaN();

    double endpoint_ratio =
        std::numeric_limits<double>::quiet_NaN();

    std::string classification = "unclassified";
};

Result analyze(
    const Characteristic& characteristic,
    double critical_flow)
{
    Result result;

    std::vector<Evaluation> samples;

    for (double flow = minimum_flow;
         flow < critical_flow;
         flow += scan_step) {
        try {
            samples.push_back(
                characteristic.evaluate(flow));
            ++result.valid_points;
        }
        catch (const std::exception&) {
            ++result.failed_points;
        }
    }

    try {
        const auto endpoint =
            characteristic.evaluate(critical_flow);

        result.endpoint_ratio =
            endpoint.pressure_ratio;

        samples.push_back(endpoint);
        ++result.valid_points;
    }
    catch (const std::exception&) {
        ++result.failed_points;
    }

    if (result.failed_points > 0) {
        result.classification = "incomplete";
        return result;
    }

    if (samples.size() < 3) {
        result.classification = "insufficient";
        return result;
    }

    // Identify all sampled interior local maxima.
    std::vector<std::size_t> candidates;

    for (std::size_t i = 1;
         i + 1 < samples.size();
         ++i) {
        if (samples[i].pressure_ratio >
                samples[i - 1].pressure_ratio &&
            samples[i].pressure_ratio >
                samples[i + 1].pressure_ratio) {
            candidates.push_back(i);
        }
    }

    result.candidate_peaks =
        static_cast<int>(candidates.size());

    if (candidates.empty()) {
        const auto best = std::max_element(
            samples.begin(),
            samples.end(),
            [](const auto& a, const auto& b) {
                return a.pressure_ratio <
                       b.pressure_ratio;
            });

        result.classification =
            best == samples.begin()
                ? "low_boundary"
                : best == samples.end() - 1
                    ? "capacity_boundary"
                    : "unresolved";

        return result;
    }

    // Refine every candidate, retaining the highest.
    for (auto index : candidates) {
        const auto refined = maximize(
            characteristic,
            samples[index - 1].flow,
            samples[index + 1].flow
        );

        if (!result.peak_found ||
            refined.pressure_ratio >
                result.peak_ratio) {
            result.peak_found = true;
            result.peak_flow = refined.flow;
            result.peak_ratio =
                refined.pressure_ratio;
            result.peak_diffusion_loss =
                refined.diffusion_loss;
            result.peak_diffusion_ratio =
                refined.diffusion_ratio;
        }
    }

    result.classification =
        candidates.size() == 1
            ? "interior_peak"
            : "multiple_peaks";

    return result;
}

} // namespace


int main(int argc, char** argv)
{
    try {
        const std::filesystem::path output =
            argc > 1
                ? argv[1]
                : "benchmarks/results-geometry-sensitivity";

        std::filesystem::create_directories(output);

        std::ofstream csv(output / "geometry_sensitivity.csv");

        if (!csv)
            throw std::runtime_error("Cannot open output CSV");

        csv << std::setprecision(17);

        csv
            << "rpm,variable,scale,"
            << "outlet_radius,outlet_width,"
            << "outlet_effective_area,"
            << "classification,candidate_peaks,"
            << "valid_points,failed_points,"
            << "peak_flow,peak_ratio,"
            << "critical_flow,endpoint_ratio,"
            << "relative_interval,"
            << "peak_diffusion_loss,peak_diffusion_ratio,"
            << "choking_mass_residual\n";

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
            throw std::runtime_error("Inlet EOS evaluation failed");

        const OperatingIncidenceConditions design{
            5.0, 12000.0, inlet.value.rho, 20.0
        };

        const double blade_angle =
            calibrate_inlet_blade_angle_deg(
                geometry, design);

        constexpr std::array<double, 5> scales_to_test{
            0.90, 0.95, 1.00, 1.05, 1.10
        };

        int combinations = 0;
        int incomplete = 0;
        int interior = 0;

        for (double speed : speeds) {
            for (int variable = 0; variable < 2; ++variable) {
                for (double scale : scales_to_test) {
                    // Do not duplicate the baseline geometry.
                    if (variable == 1 && scale == 1.0)
                        continue;

                    ImpellerSystemGeometry g = geometry;

                    const char* name =
                        variable == 0
                            ? "outlet_radius"
                            : "outlet_width";

                    if (variable == 0) {
                        g.exit.diameter *= scale;
                        g.impeller.exit_diameter *= scale;
                    }
                    else {
                        g.exit.width *= scale;
                    }

                    const double outlet_radius =
                        0.5 * g.exit.diameter;

                    const double outlet_width =
                        g.exit.width;

                    const double effective_area =
                        effective_flow_area(
                            impeller_exit_area(g.exit),
                            g.exit.blockage);

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
                            methane, methane, g, input);

                    const double critical =
                        choking.choking.coupled.critical_mass_flow;

                    Characteristic characteristic(
                        methane,
                        speed,
                        blade_angle,
                        1.0,
                        1.0,
                        1.5,
                        g
                    );

                    Result result;

                    try {
                        result = analyze(
                            characteristic, critical);
                    }
                    catch (const std::exception& e) {
                        result.classification =
                            "refinement_failure";

                        std::cerr
                            << "FAIL RPM=" << speed
                            << " variable=" << name
                            << " scale=" << scale
                            << ": " << e.what()
                            << '\n';
                    }

                    const double interval =
                        result.peak_found
                            ? 1.0 - result.peak_flow / critical
                            : std::numeric_limits<double>::quiet_NaN();

                    csv
                        << speed << ','
                        << name << ','
                        << scale << ','
                        << outlet_radius << ','
                        << outlet_width << ','
                        << effective_area << ','
                        << result.classification << ','
                        << result.candidate_peaks << ','
                        << result.valid_points << ','
                        << result.failed_points << ','
                        << result.peak_flow << ','
                        << result.peak_ratio << ','
                        << critical << ','
                        << result.endpoint_ratio << ','
                        << interval << ','
                        << result.peak_diffusion_loss << ','
                        << result.peak_diffusion_ratio << ','
                        << choking.mass_flow_residual
                        << '\n';

                    ++combinations;

                    if (result.peak_found)
                        ++interior;

                    if (result.classification == "incomplete" ||
                        result.classification == "refinement_failure" ||
                        result.classification == "unresolved")
                        ++incomplete;
                }
            }

            csv.flush();

            std::cout
                << "RPM=" << speed
                << " combinations=" << combinations
                << " interior=" << interior
                << " incomplete=" << incomplete
                << '\n';
        }

        std::cout
            << "\nSTAGE 7J SUMMARY\n"
            << "combinations=" << combinations << '\n'
            << "interior_peaks=" << interior << '\n'
            << "incomplete=" << incomplete << '\n'
            << "output=" << output << '\n';

        return incomplete == 0 ? 0 : 2;
    }
    catch (const std::exception& e) {
        std::cerr << "FATAL: " << e.what() << '\n';
        return 1;
    }
}
