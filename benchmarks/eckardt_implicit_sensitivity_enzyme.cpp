#include <turbocompgpu/operating_point.hpp>
#include <turbocompgpu/station.hpp>

#include <thermomodel/ideal_gas.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>

// Enzyme replaces these calls during the Clang Enzyme pass.
extern int enzyme_dup;
extern int enzyme_const;

extern "C" double __enzyme_fwddiff(void*, ...);
extern "C" void __enzyme_autodiff(void*, ...);

namespace {

using Vector = std::array<double, 24>;
using ScalarKernel = double (*)(const double*);

enum Index : std::size_t {
    Temperature = 0,
    LogPressure,
    Density,
    Enthalpy,
    Entropy,
    MassFlow,
    SpeedRpm,
    Slip,
    Prewhirl,
    EntropyGeneration,
    ExitDiameter,
    ExitWidth,
    ExitBlockage,
    InletDiameter,
    ExitAngleDegrees,
    InletStagnationEnthalpy,
    InletEntropy,
    PressureScale,
    WorkScale,
    EntropyScale,
    InletStagnationPressure,
    InletStagnationTemperature,
    Cp,
    GasConstant
};

constexpr std::array<Index, 5> active_parameters{
    Slip, ExitAngleDegrees, ExitWidth, ExitBlockage,
    EntropyGeneration
};

constexpr std::array<const char*, 5> active_parameter_names{
    "slip_factor", "exit_blade_angle_deg", "exit_width_m",
    "exit_blockage", "entropy_generation_j_kg_k"
};

struct Point {
    double mass_flow;
    double rpm;
};

struct Parameters {
    double mass_flow;
    double rpm;
    double slip;
    double prewhirl;
    double entropy_generation;
    double exit_diameter;
    double exit_width;
    double exit_blockage;
    double inlet_diameter;
    double exit_angle_deg;
    double inlet_h01;
    double inlet_entropy;
    double pressure_scale;
    double work_scale;
    double entropy_scale;
    double inlet_p01;
    double inlet_T01;
    double cp;
    double gas_constant;
};

struct OperatingCase {
    turbocompgpu::ImpellerSystemGeometry geometry;
    turbocompgpu::OperatingPointRequest request;
    turbocompgpu::NonlinearResidualScales scales;
    turbocompgpu::OperatingPointResult solution;
    Vector packed;
    thermomodel::IdealGas air;
};

double normalized_energy(const double* z)
{
    constexpr double pi = std::numbers::pi_v<double>;
    const double area = pi * z[ExitDiameter] * z[ExitWidth] *
        (1.0 - z[ExitBlockage]);
    const double vm2 = z[MassFlow] / (z[Density] * area);
    const double omega = 2.0 * pi * z[SpeedRpm] / 60.0;
    const double U2 = omega * z[ExitDiameter] / 2.0;
    const double U1 = omega * z[InletDiameter] / 2.0;
    const double beta = z[ExitAngleDegrees] * pi / 180.0;
    const double whirl = z[Slip] *
        (U2 - vm2 * std::cos(beta) / std::sin(beta));
    const double euler_work = U2 * whirl - U1 * z[Prewhirl];
    const double h02 = z[Enthalpy] +
        0.5 * (vm2 * vm2 + whirl * whirl);
    return (h02 - z[InletStagnationEnthalpy] - euler_work) /
        z[WorkScale];
}

double normalized_entropy(const double* z)
{
    return (z[Entropy] - z[InletEntropy] - z[EntropyGeneration]) /
        z[EntropyScale];
}

double total_to_total_pressure_ratio(const double* z)
{
    constexpr double pi = std::numbers::pi_v<double>;
    const double p2 = z[PressureScale] * std::exp(z[LogPressure]);
    const double area = pi * z[ExitDiameter] * z[ExitWidth] *
        (1.0 - z[ExitBlockage]);
    const double vm2 = z[MassFlow] / (z[Density] * area);
    const double omega = 2.0 * pi * z[SpeedRpm] / 60.0;
    const double U2 = omega * z[ExitDiameter] / 2.0;
    const double beta = z[ExitAngleDegrees] * pi / 180.0;
    const double whirl = z[Slip] *
        (U2 - vm2 * std::cos(beta) / std::sin(beta));
    const double T02 = z[Temperature] +
        (vm2 * vm2 + whirl * whirl) / (2.0 * z[Cp]);
    const double gamma = z[Cp] / (z[Cp] - z[GasConstant]);
    const double p02 = p2 * std::pow(T02 / z[Temperature],
                                     gamma / (gamma - 1.0));
    return p02 / z[InletStagnationPressure];
}

double total_to_total_efficiency(const double* z)
{
    const double pr = total_to_total_pressure_ratio(z);
    const double exponent = z[GasConstant] / z[Cp];
    const double T02s = z[InletStagnationTemperature] *
        std::pow(pr, exponent);
    const double T02 = z[Temperature] +
        (z[MassFlow] / (z[Density] *
            (std::numbers::pi_v<double> * z[ExitDiameter] * z[ExitWidth] *
             (1.0 - z[ExitBlockage])))) *
        (z[MassFlow] / (z[Density] *
            (std::numbers::pi_v<double> * z[ExitDiameter] * z[ExitWidth] *
             (1.0 - z[ExitBlockage])))) /
        (2.0 * z[Cp]);
    // Add the whirl kinetic energy omitted from the radial term above.
    const double omega = 2.0 * std::numbers::pi_v<double> *
        z[SpeedRpm] / 60.0;
    const double U2 = omega * z[ExitDiameter] / 2.0;
    const double beta = z[ExitAngleDegrees] * std::numbers::pi_v<double> / 180.0;
    const double vm2 = z[MassFlow] /
        (z[Density] * std::numbers::pi_v<double> * z[ExitDiameter] *
         z[ExitWidth] * (1.0 - z[ExitBlockage]));
    const double whirl = z[Slip] *
        (U2 - vm2 * std::cos(beta) / std::sin(beta));
    const double T02_full = T02 + whirl * whirl / (2.0 * z[Cp]);
    return (T02s - z[InletStagnationTemperature]) /
        (T02_full - z[InletStagnationTemperature]);
}

double lagrangian_pressure(const double* z, double lambda_energy,
                           double lambda_entropy)
{
    return total_to_total_pressure_ratio(z) -
        lambda_energy * normalized_energy(z) -
        lambda_entropy * normalized_entropy(z);
}

double lagrangian_efficiency(const double* z, double lambda_energy,
                             double lambda_entropy)
{
    return total_to_total_efficiency(z) -
        lambda_energy * normalized_energy(z) -
        lambda_entropy * normalized_entropy(z);
}

template <ScalarKernel Function>
double forward(const Vector& z, const Vector& dz)
{
    return __enzyme_fwddiff(
        reinterpret_cast<void*>(Function), enzyme_dup,
        z.data(), dz.data());
}

using LagrangianKernel = double (*)(const double*, double, double);

template <LagrangianKernel Function>
Vector reverse(const Vector& z, double lambda_energy,
               double lambda_entropy)
{
    Vector gradient{};
    __enzyme_autodiff(
        reinterpret_cast<void*>(Function),
        enzyme_dup, z.data(), gradient.data(),
        enzyme_const, lambda_energy,
        enzyme_const, lambda_entropy);
    return gradient;
}

bool solve_2x2(const std::array<std::array<double, 2>, 2>& A,
               const std::array<double, 2>& rhs,
               std::array<double, 2>& answer)
{
    const double determinant = A[0][0] * A[1][1] -
        A[0][1] * A[1][0];
    const double scale = std::max({1.0, std::abs(A[0][0] * A[1][1]),
                                   std::abs(A[0][1] * A[1][0])});
    if (!std::isfinite(determinant) ||
        std::abs(determinant) <= 1.0e-13 * scale)
        return false;
    answer = {
        (rhs[0] * A[1][1] - A[0][1] * rhs[1]) / determinant,
        (A[0][0] * rhs[1] - rhs[0] * A[1][0]) / determinant
    };
    return std::isfinite(answer[0]) && std::isfinite(answer[1]);
}

bool solve_transpose_2x2(const std::array<std::array<double, 2>, 2>& A,
                         const std::array<double, 2>& rhs,
                         std::array<double, 2>& answer)
{
    const std::array<std::array<double, 2>, 2> AT{{
        {A[0][0], A[1][0]}, {A[0][1], A[1][1]}
    }};
    return solve_2x2(AT, rhs, answer);
}

turbocompgpu::ImpellerSystemGeometry make_geometry()
{
    constexpr double d2 = 0.400;
    constexpr double b2 = 0.026;
    constexpr double blade_count = 20.0;
    constexpr double blade_thickness_m = 0.00108;
    const double blockage = blade_count * blade_thickness_m /
        (2.0 * std::numbers::pi_v<double> * (d2 / 2.0));
    return {
        {0.120, 0.280, 0.0},
        {d2, b2, blockage},
        {d2, 60.0,
         2.0 * std::sqrt(0.5 * (0.140 * 0.140 + 0.060 * 0.060))}
    };
}

turbocompgpu::OperatingPointRequest make_request(
    Point point, const thermomodel::IdealGas& air,
    const turbocompgpu::ImpellerSystemGeometry& geometry,
    double slip, double entropy_generation)
{
    constexpr double T01 = 288.15;
    constexpr double p01 = 101000.0;
    const double gamma = air.parameters().cp /
        (air.parameters().cp - air.parameters().R);
    const double inlet_area = turbocompgpu::effective_flow_area(
        turbocompgpu::inlet_annulus_area(geometry.inlet),
        geometry.inlet.blockage);
    double lo = 0.0;
    double hi = 1.0;
    for (int i = 0; i < 80; ++i) {
        const double mach = 0.5 * (lo + hi);
        const double T = T01 /
            (1.0 + 0.5 * (gamma - 1.0) * mach * mach);
        const double p = p01 * std::pow(T / T01,
                                         gamma / (gamma - 1.0));
        const auto state = turbocompgpu::evaluate(air, T, p);
        const double flow = state.rho * inlet_area * mach *
            std::sqrt(gamma * air.parameters().R * T);
        if (flow < point.mass_flow)
            lo = mach;
        else
            hi = mach;
    }
    const double mach = 0.5 * (lo + hi);
    const double T1 = T01 /
        (1.0 + 0.5 * (gamma - 1.0) * mach * mach);
    const double p1 = p01 * std::pow(T1 / T01,
                                     gamma / (gamma - 1.0));
    return {
        {T1, p1, point.mass_flow, point.rpm, slip, 0.0,
         entropy_generation},
        {T1 * 1.15, std::log(1.5 * p1 / 50000.0)}
    };
}

Vector pack(const OperatingCase& op_case)
{
    const auto& request = op_case.request;
    const auto& result = op_case.solution;
    const auto state = turbocompgpu::evaluate(
        op_case.air, result.outlet_temperature, result.outlet_pressure);
    const auto inlet = turbocompgpu::evaluate(
        op_case.air, request.operating.inlet_temperature,
        request.operating.inlet_pressure);
    const double area1 = turbocompgpu::effective_flow_area(
        turbocompgpu::inlet_annulus_area(op_case.geometry.inlet),
        op_case.geometry.inlet.blockage);
    const double vm1 = request.operating.mass_flow /
        (inlet.rho * area1);
    const turbocompgpu::FlowStation inlet_station{
        inlet, {vm1, request.operating.inlet_prewhirl}};
    const auto inlet_total = turbocompgpu::evaluate_stagnation_state(
        op_case.air, inlet_station).stagnation_state;

    Vector z{};
    z[Temperature] = result.outlet_temperature;
    z[LogPressure] = std::log(result.outlet_pressure /
                              op_case.scales.pressure);
    z[Density] = state.rho;
    z[Enthalpy] = state.h;
    z[Entropy] = state.s;
    z[MassFlow] = request.operating.mass_flow;
    z[SpeedRpm] = request.operating.rotational_speed_rpm;
    z[Slip] = request.operating.slip_factor;
    z[Prewhirl] = request.operating.inlet_prewhirl;
    z[EntropyGeneration] = request.operating.entropy_generation;
    z[ExitDiameter] = op_case.geometry.exit.diameter;
    z[ExitWidth] = op_case.geometry.exit.width;
    z[ExitBlockage] = op_case.geometry.exit.blockage;
    z[InletDiameter] = op_case.geometry.impeller.inlet_diameter;
    z[ExitAngleDegrees] =
        op_case.geometry.impeller.exit_blade_angle_deg;
    z[InletStagnationEnthalpy] =
        turbocompgpu::evaluate_stagnation_state(op_case.air,
                                                inlet_station)
            .stagnation_enthalpy;
    z[InletEntropy] = inlet.s;
    z[PressureScale] = op_case.scales.pressure;
    z[WorkScale] = op_case.scales.specific_work;
    z[EntropyScale] = op_case.scales.specific_entropy;
    z[InletStagnationPressure] = inlet_total.p;
    z[InletStagnationTemperature] = inlet_total.T;
    z[Cp] = op_case.air.parameters().cp;
    z[GasConstant] = op_case.air.parameters().R;
    return z;
}

OperatingCase solve_case(Point point, double slip, double beta_deg,
                         double exit_width, double blockage,
                         double entropy_generation)
{
    OperatingCase result;
    result.air = thermomodel::IdealGas{};
    result.geometry = make_geometry();
    result.geometry.exit.width = exit_width;
    result.geometry.exit.blockage = blockage;
    result.geometry.impeller.exit_blade_angle_deg = beta_deg;
    result.scales = {50000.0, 100.0, 300000.0};
    result.request = make_request(point, result.air, result.geometry,
                                  slip, entropy_generation);
    result.solution = turbocompgpu::evaluate_operating_point(
        result.air, result.geometry, result.request, result.scales);
    if (!(result.solution.converged &&
          result.solution.independently_validated &&
          result.solution.basic_physics_valid)) {
        throw std::runtime_error("Operating-point solve failed at " +
                                 std::to_string(point.mass_flow) +
                                 " kg/s: " + result.solution.message);
    }
    result.packed = pack(result);
    return result;
}

double parameter_value(const OperatingCase& op_case, Index parameter)
{
    switch (parameter) {
    case Slip: return op_case.request.operating.slip_factor;
    case ExitAngleDegrees:
        return op_case.geometry.impeller.exit_blade_angle_deg;
    case ExitWidth: return op_case.geometry.exit.width;
    case ExitBlockage: return op_case.geometry.exit.blockage;
    case EntropyGeneration:
        return op_case.request.operating.entropy_generation;
    default: throw std::invalid_argument("Unsupported active parameter");
    }
}

double finite_difference(ScalarKernel objective,
                         const OperatingCase& op_case,
                         Index parameter)
{
    const double step = parameter == ExitAngleDegrees ? 1.0e-3
        : parameter == ExitWidth ? 1.0e-6
        : parameter == EntropyGeneration ? 1.0e-3
        : parameter == ExitBlockage ? 1.0e-5
        : 1.0e-5;
    double plus_slip = op_case.request.operating.slip_factor;
    double minus_slip = plus_slip;
    double plus_beta = op_case.geometry.impeller.exit_blade_angle_deg;
    double minus_beta = plus_beta;
    double plus_width = op_case.geometry.exit.width;
    double minus_width = plus_width;
    double plus_blockage = op_case.geometry.exit.blockage;
    double minus_blockage = plus_blockage;
    double plus_entropy = op_case.request.operating.entropy_generation;
    double minus_entropy = plus_entropy;
    switch (parameter) {
    case Slip: plus_slip += step; minus_slip -= step; break;
    case ExitAngleDegrees: plus_beta += step; minus_beta -= step; break;
    case ExitWidth: plus_width += step; minus_width -= step; break;
    case ExitBlockage: plus_blockage += step; minus_blockage -= step; break;
    case EntropyGeneration: plus_entropy += step; minus_entropy -= step; break;
    default: throw std::invalid_argument("Unsupported finite-difference parameter");
    }
    const Point point{op_case.request.operating.mass_flow,
                      op_case.request.operating.rotational_speed_rpm};
    const auto plus = solve_case(point, plus_slip, plus_beta, plus_width,
                                 plus_blockage, plus_entropy);
    const auto minus = solve_case(point, minus_slip, minus_beta, minus_width,
                                  minus_blockage, minus_entropy);
    return (objective(plus.packed.data()) - objective(minus.packed.data())) /
        (2.0 * step);
}

Vector tangent_state_T(const thermomodel::Derivatives& derivatives)
{
    Vector tangent{};
    tangent[Temperature] = 1.0;
    tangent[Density] = derivatives.drho_dT_p;
    tangent[Enthalpy] = derivatives.dh_dT_p;
    tangent[Entropy] = derivatives.ds_dT_p;
    return tangent;
}

Vector tangent_state_logp(const Vector& z,
                          const thermomodel::Derivatives& derivatives)
{
    Vector tangent{};
    const double pressure = z[PressureScale] * std::exp(z[LogPressure]);
    tangent[LogPressure] = 1.0;
    tangent[Density] = derivatives.drho_dp_T * pressure;
    tangent[Enthalpy] = derivatives.dh_dp_T * pressure;
    tangent[Entropy] = derivatives.ds_dp_T * pressure;
    return tangent;
}

void calculate_sensitivities(const OperatingCase& op_case)
{
    const Vector& z = op_case.packed;
    const auto derivatives_result = op_case.air.derivatives_Tp(
        z[Temperature],
        z[PressureScale] * std::exp(z[LogPressure]));
    if (!derivatives_result)
        throw std::runtime_error("Ideal-gas derivative evaluation failed");
    const auto& eos = derivatives_result.value;
    const Vector dT = tangent_state_T(eos);
    const Vector dlogp = tangent_state_logp(z, eos);
    const std::array<std::array<double, 2>, 2> Fx{{
        {forward<normalized_energy>(z, dT),
         forward<normalized_energy>(z, dlogp)},
        {forward<normalized_entropy>(z, dT),
         forward<normalized_entropy>(z, dlogp)}
    }};

    const auto report_objective = [&]<ScalarKernel Objective,
                                      LagrangianKernel Lagrangian>(
        const char* name) {
        const std::array<double, 2> Jx{
            forward<Objective>(z, dT),
            forward<Objective>(z, dlogp)
        };
        std::array<double, 2> lambda{};
        if (!solve_transpose_2x2(Fx, Jx, lambda))
            throw std::runtime_error("Singular implicit adjoint system");

        const Vector gradient = reverse<Lagrangian>(
            z, lambda[0], lambda[1]);

        for (std::size_t k = 0; k < active_parameters.size(); ++k) {
            Vector dparameter{};
            dparameter[active_parameters[k]] = 1.0;
            const std::array<double, 2> Ftheta{
                forward<normalized_energy>(z, dparameter),
                forward<normalized_entropy>(z, dparameter)
            };
            std::array<double, 2> dx{};
            if (!solve_2x2(Fx, {-Ftheta[0], -Ftheta[1]}, dx))
                throw std::runtime_error("Singular implicit TLM system");

            Vector total_tangent{};
            for (std::size_t i = 0; i < total_tangent.size(); ++i) {
                total_tangent[i] = dparameter[i] +
                    dx[0] * dT[i] + dx[1] * dlogp[i];
            }
            const double tlm = forward<Objective>(z, total_tangent);
            const double adjoint = gradient[active_parameters[k]];
            const double finite_difference_value = finite_difference(
                Objective, op_case, active_parameters[k]);
            std::cout << std::setprecision(10)
                      << op_case.request.operating.mass_flow << ','
                      << name << ','
                      << active_parameter_names[k] << ','
                      << tlm << ',' << adjoint << ','
                      << finite_difference_value << ','
                      << (tlm - adjoint) << ','
                      << (tlm - finite_difference_value) << '\n';
            const double adjoint_tolerance = 2.0e-8 *
                std::max({1.0, std::abs(tlm), std::abs(adjoint)});
            const double finite_difference_tolerance = 2.0e-3 *
                std::max({1.0e-6, std::abs(tlm),
                          std::abs(finite_difference_value)});
            if (std::abs(tlm - adjoint) > adjoint_tolerance ||
                std::abs(tlm - finite_difference_value) >
                    finite_difference_tolerance) {
                throw std::runtime_error(
                    std::string("Sensitivity verification failed for ") +
                    name + " / " + active_parameter_names[k]);
            }
        }
    };
    report_objective.template operator()<total_to_total_pressure_ratio,
                                           lagrangian_pressure>("PR_tt");
    report_objective.template operator()<total_to_total_efficiency,
                                           lagrangian_efficiency>("eta_tt");
}

} // namespace

int main()
{
    try {
        constexpr double rpm = 14000.0;
        const double slip = 1.0 -
            std::sqrt(std::cos(30.0 * std::numbers::pi_v<double> / 180.0)) /
            std::pow(20.0, 0.7);
        std::cout << "flow_kg_s,objective,parameter,TLM,adjoint,finite_difference,"
                     "TLM_minus_adjoint,TLM_minus_FD\n";
        for (const double flow : {3.99, 4.56, 5.11}) {
            const auto op_case = solve_case(
                {flow, rpm}, slip, 60.0, 0.026,
                make_geometry().exit.blockage, 10.0);
            calculate_sensitivities(op_case);
        }
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "Enzyme implicit sensitivity failed: "
                  << error.what() << '\n';
        return 1;
    }
}
