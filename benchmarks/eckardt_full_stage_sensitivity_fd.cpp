#define TURBOCOMP_ECKARDT_NO_MAIN
#include "eckardt_diffuser_friction_sweep.cpp"

#include <array>


namespace {

using SensitivityVector = std::array<double, 11>;
using ScalarKernel = double (*)(const double*, double, double);

enum Input : std::size_t {
    ExitAngle = 0,
    ExitWidth,
    ExitBlockage,
    MixingChi,
    SkinScale,
    BladeScale,
    IncidenceScale,
    ClearanceScale,
    MixingScale,
    DiskScale,
    DiffuserFrictionScale
};

constexpr std::array<const char*, 11> parameter_names{
    "exit_blade_angle_deg", "exit_width_m", "exit_blockage",
    "mixing_chi", "skin_friction_scale", "blade_loading_scale",
    "incidence_scale", "clearance_scale", "mixing_scale",
    "disk_friction_scale", "diffuser_friction_scale"
};

constexpr std::array<Input, 11> differentiated_inputs{
    ExitAngle, ExitWidth, ExitBlockage, MixingChi, SkinScale, BladeScale,
    IncidenceScale, ClearanceScale, MixingScale, DiskScale,
    DiffuserFrictionScale
};

constexpr std::array<double, 11> finite_difference_steps{
    1.0e-4, 1.0e-7, 1.0e-6, 1.0e-5, 1.0e-5, 1.0e-5,
    1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5
};

template <bool WantEfficiency>
double full_stage(const SensitivityVector& z, double flow, double rpm)
{
    const MapPoint point{rpm, flow, 1.0};
    const thermomodel::IdealGas air;
    auto g = geometry();
    g.impeller.exit_blade_angle_deg = z[ExitAngle];
    g.exit.width = z[ExitWidth];
    g.exit.blockage = z[ExitBlockage];

    LossModelOptions losses;
    losses.mixing_chi = z[MixingChi];
    losses.skin_friction_multiplier = z[SkinScale];
    losses.blade_loading_multiplier = z[BladeScale];
    losses.incidence_multiplier = z[IncidenceScale];
    losses.clearance_multiplier = z[ClearanceScale];
    losses.mixing_multiplier = z[MixingScale];
    losses.disk_friction_multiplier = z[DiskScale];

    ThroatCapacityCache capacity_cache;
    const auto result = predict(
        point, air, g, -1.0, true, capacity_cache,
        true, true, losses, z[DiffuserFrictionScale]);
    if constexpr (WantEfficiency)
        return result.total_to_total_efficiency;
    else
        return result.stage_pr;
}

double stage_efficiency(const double* z, double flow, double rpm)
{
    SensitivityVector values{};
    std::copy_n(z, differentiated_inputs.size(), values.begin());
    return full_stage<true>(values, flow, rpm);
}

double stage_pressure_ratio(const double* z, double flow, double rpm)
{
    SensitivityVector values{};
    std::copy_n(z, differentiated_inputs.size(), values.begin());
    return full_stage<false>(values, flow, rpm);
}

SensitivityVector baseline()
{
    const auto g = geometry();
    SensitivityVector z{};
    z[ExitAngle] = g.impeller.exit_blade_angle_deg;
    z[ExitWidth] = g.exit.width;
    z[ExitBlockage] = g.exit.blockage;
    z[MixingChi] = 0.15;
    z[SkinScale] = 1.0;
    z[BladeScale] = 1.0;
    z[IncidenceScale] = 1.0;
    z[ClearanceScale] = 1.0;
    z[MixingScale] = 1.0;
    z[DiskScale] = 1.0;
    z[DiffuserFrictionScale] = 1.0;
    return z;
}

double centered_difference(ScalarKernel objective,
                           const SensitivityVector& z,
                           double flow,
                           double rpm,
                           Input parameter,
                           double step)
{
    auto plus = z;
    auto minus = z;
    plus[parameter] += step;
    minus[parameter] -= step;
    return (objective(plus.data(), flow, rpm) -
            objective(minus.data(), flow, rpm)) /
        (2.0 * step);
}

struct PointReport {
    double value;
    SensitivityVector gradient;
};

template <ScalarKernel Objective>
PointReport report_point(double flow, const char* objective_name)
{
    constexpr double rpm = 14000.0;
    const auto z = baseline();
    const double value = Objective(z.data(), flow, rpm);
    SensitivityVector gradient{};
    std::cout << std::setprecision(11)
              << "# flow=" << flow << " kg/s, objective=" << objective_name
              << ", base=" << value << '\n';
    for (std::size_t i = 0; i < differentiated_inputs.size(); ++i) {
        const double derivative = centered_difference(
            Objective, z, flow, rpm, differentiated_inputs[i],
            finite_difference_steps[i]);
        if (!std::isfinite(derivative))
            throw std::runtime_error(
                std::string("Non-finite full-stage finite difference for ") +
                objective_name + " / " + parameter_names[i]);
        gradient[differentiated_inputs[i]] = derivative;
        std::cout << flow << ',' << objective_name << ',' << parameter_names[i]
                  << ',' << derivative << '\n';
    }
    return {value, gradient};
}

template <std::size_t N>
double least_squares_slope(const std::array<double, N>& x,
                           const std::array<double, N>& y)
{
    double x_mean = 0.0;
    double y_mean = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        x_mean += x[i] / static_cast<double>(N);
        y_mean += y[i] / static_cast<double>(N);
    }
    double numerator = 0.0;
    double denominator = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        numerator += (x[i] - x_mean) * (y[i] - y_mean);
        denominator += (x[i] - x_mean) * (x[i] - x_mean);
    }
    return numerator / denominator;
}

double measured_pressure_ratio(double flow)
{
    // 14,000 rpm experimental speed-line points from the published map.
    // Interpolate only within the measured range.
    constexpr std::array<double, 4> measured_flow{3.50, 4.50, 5.30, 6.10};
    constexpr std::array<double, 4> measured_pr{1.9400, 1.9250, 1.8800, 1.8050};
    for (std::size_t i = 1; i < measured_flow.size(); ++i) {
        if (flow <= measured_flow[i]) {
            const double fraction = (flow - measured_flow[i - 1]) /
                (measured_flow[i] - measured_flow[i - 1]);
            return measured_pr[i - 1] + fraction *
                (measured_pr[i] - measured_pr[i - 1]);
        }
    }
    throw std::out_of_range("Flow is outside the measured 14,000 rpm speed line");
}

struct SweepMetrics {
    double eta_bias_pp;
    double eta_mae_pp;
    double eta_slope_pp_per_kg_s;
    double eta_slope_error_pp_per_kg_s;
    double pr_bias_percent;
    double pr_mape_percent;
};

SweepMetrics evaluate_scenario(double skin_scale, double mixing_chi)
{
    constexpr double rpm = 14000.0;
    constexpr std::array<double, 3> flows{3.99, 4.56, 5.11};
    // Corrected approximate digitization of the measured Figure 6(d) symbols.
    constexpr std::array<double, 3> measured_eta{0.865, 0.877, 0.890};
    std::array<double, 3> predicted_eta{};
    std::array<double, 3> predicted_pr{};
    double eta_bias = 0.0;
    double eta_absolute_error = 0.0;
    double pr_bias = 0.0;
    double pr_absolute_error = 0.0;
    auto z = baseline();
    z[SkinScale] = skin_scale;
    z[MixingChi] = mixing_chi;
    for (std::size_t i = 0; i < flows.size(); ++i) {
        predicted_eta[i] = stage_efficiency(z.data(), flows[i], rpm);
        predicted_pr[i] = stage_pressure_ratio(z.data(), flows[i], rpm);
        const double eta_error_pp = 100.0 * (predicted_eta[i] - measured_eta[i]);
        const double pr_error_percent = 100.0 *
            (predicted_pr[i] / measured_pressure_ratio(flows[i]) - 1.0);
        eta_bias += eta_error_pp / static_cast<double>(flows.size());
        eta_absolute_error += std::abs(eta_error_pp) /
            static_cast<double>(flows.size());
        pr_bias += pr_error_percent / static_cast<double>(flows.size());
        pr_absolute_error += std::abs(pr_error_percent) /
            static_cast<double>(flows.size());
    }
    const double measured_eta_slope = least_squares_slope(flows, measured_eta);
    const double model_eta_slope = least_squares_slope(flows, predicted_eta);
    return {
        eta_bias,
        eta_absolute_error,
        100.0 * model_eta_slope,
        100.0 * (model_eta_slope - measured_eta_slope),
        pr_bias,
        pr_absolute_error
    };
}

void print_sweep_row(const char* sweep, double skin_scale, double mixing_chi)
{
    const auto m = evaluate_scenario(skin_scale, mixing_chi);
    std::cout << sweep << ',' << skin_scale << ',' << mixing_chi << ','
              << m.eta_bias_pp << ',' << m.eta_mae_pp << ','
              << m.eta_slope_pp_per_kg_s << ','
              << m.eta_slope_error_pp_per_kg_s << ','
              << m.pr_bias_percent << ',' << m.pr_mape_percent << '\n';
}

} // namespace

int main()
{
    try {
        std::cout << std::setprecision(11);
        std::cout << "flow_kg_s,objective,parameter,centered_finite_difference\n";

        // Four retained 14,000 rpm points; omit the 5.69 kg/s symbol.
        constexpr std::array<double, 4> flows{3.30, 3.99, 4.56, 5.11};
        std::array<PointReport, 4> eta_reports{};
        std::array<PointReport, 4> pr_reports{};

        for (std::size_t i = 0; i < flows.size(); ++i) {
            eta_reports[i] = report_point<&stage_efficiency>(
                flows[i], "stage_eta_tt");
            pr_reports[i] = report_point<&stage_pressure_ratio>(
                flows[i], "stage_PR_tt");
        }

        std::array<double, 4> eta_values{};
        std::array<double, 4> pr_values{};
        for (std::size_t i = 0; i < flows.size(); ++i) {
            eta_values[i] = eta_reports[i].value;
            pr_values[i] = pr_reports[i].value;
        }

        const double eta_slope = least_squares_slope(flows, eta_values);
        const double pr_slope = least_squares_slope(flows, pr_values);

        std::cout << "\n# full-stage slope sensitivities over four retained 14,000 rpm points\n"
                     "# eta: percentage points/(kg/s); PR: ratio/(kg/s)\n"
                     "objective,base_slope,parameter,slope_derivative\n"
                  << "stage_eta_tt," << 100.0 * eta_slope << ",base,0\n"
                  << "stage_PR_tt," << pr_slope << ",base,0\n";

        for (std::size_t parameter = 0;
             parameter < differentiated_inputs.size(); ++parameter) {
            std::array<double, 4> eta_gradient_vs_flow{};
            std::array<double, 4> pr_gradient_vs_flow{};

            for (std::size_t i = 0; i < flows.size(); ++i) {
                const Input input = differentiated_inputs[parameter];
                eta_gradient_vs_flow[i] = eta_reports[i].gradient[input];
                pr_gradient_vs_flow[i] = pr_reports[i].gradient[input];
            }

            std::cout << "stage_eta_tt," << 100.0 * eta_slope << ','
                      << parameter_names[parameter] << ','
                      << 100.0 * least_squares_slope(
                             flows, eta_gradient_vs_flow) << '\n'
                      << "stage_PR_tt," << pr_slope << ','
                      << parameter_names[parameter] << ','
                      << least_squares_slope(
                             flows, pr_gradient_vs_flow) << '\n';
        }

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "Full-stage finite-difference sensitivity failed: "
                  << error.what() << '\n';
        return 1;
    }
}
