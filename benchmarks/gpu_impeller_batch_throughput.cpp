// SPDX-License-Identifier: GPL-3.0-or-later
#include <cmath>
// ThermoGPU's CUDA PR-mixture header is also included by this host C++ file.
using std::isfinite;

#include <turbocompgpu/fast_nonlinear_solver.hpp>
#include <turbocompgpu/gpu_impeller_batch_cuda.hpp>
#include <turbocompgpu/vaneless_diffuser_friction_analytical.hpp>
#include "lin2026_natural_gas.hpp"

#include <friction/colebrook.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#ifdef TURBOCOMPGPU_BENCHMARK_OPENMP
#include <omp.h>
#endif

namespace {

using Clock = std::chrono::steady_clock;
using turbocompgpu::cuda::ImpellerBatchInput;
using turbocompgpu::cuda::ImpellerBatchOptions;
using turbocompgpu::cuda::ImpellerBatchOutput;
using turbocompgpu::cuda::ImpellerBatchStatus;

ImpellerBatchInput point(double flow)
{
    return {
        293.15, 5.0e5, flow, 14000.0, 0.90, 0.0, 0.0,
        0.030, 0.150, 0.03,
        0.200, 0.030, 0.0, 0.090, 50.0,
        5.0e4, 100.0, 5.0e5,
        0.15, 0.020, 1.1e-5, 0.0,
        310.0, 0.20,
        0.320, 512
    };
}

ImpellerBatchOutput cpu_solve(
    const ImpellerBatchInput& x,
    const thermomodel::PengRobinsonMixture& mixture,
    const turbocompgpu::transport::Iso20765Parameters& viscosity_model,
    const ImpellerBatchOptions& gpu_options)
{
    constexpr double pi = 3.141592653589793238462643383279502884;
    const auto inlet = mixture.evaluate_Tp(
        x.inlet_temperature_K, x.inlet_pressure_Pa);
    if (!inlet) throw std::runtime_error("CPU inlet EOS evaluation failed");

    const double inlet_area = (pi / 4.0) *
        (x.inlet_shroud_diameter_m * x.inlet_shroud_diameter_m -
         x.inlet_hub_diameter_m * x.inlet_hub_diameter_m) *
        (1.0 - x.inlet_blockage);
    const double vm1 = x.mass_flow_kg_s / (inlet.value.rho * inlet_area);
    const double relative_speed = std::hypot(vm1, x.inlet_prewhirl_m_s);
    const double inlet_mu = gpu_options.use_iso20765_viscosity
        ? turbocompgpu::transport::iso20765_5_dynamic_viscosity_Pa_s(
              viscosity_model, inlet.value.T, inlet.value.p, inlet.value.rho)
        : x.dynamic_viscosity_Pa_s;
    if (!(std::isfinite(inlet_mu) && inlet_mu > 0.0))
        throw std::runtime_error("CPU ISO 20765 inlet viscosity failed");
    const double reynolds = inlet.value.rho * relative_speed *
        x.hydraulic_diameter_m / inlet_mu;
    const auto friction = friction::colebrook_turbulent(
        reynolds, x.relative_roughness);
    if (!(std::isfinite(friction.darcy_friction_factor) &&
          friction.darcy_friction_factor > 0.0))
        throw std::runtime_error("CPU Colebrook evaluation failed");

    const double passage_dh = friction.darcy_friction_factor *
        (x.passage_length_m / x.hydraulic_diameter_m) *
        0.5 * relative_speed * relative_speed;
    const double passage_entropy = passage_dh / x.inlet_temperature_K;

    const turbocompgpu::ImpellerSystemGeometry geometry{
        {x.inlet_hub_diameter_m, x.inlet_shroud_diameter_m, x.inlet_blockage},
        {x.exit_diameter_m, x.exit_width_m, x.exit_blockage},
        {x.exit_diameter_m, x.exit_blade_angle_deg, x.inlet_impeller_diameter_m}
    };
    const turbocompgpu::NonlinearOperatingPoint operating{
        x.inlet_temperature_K, x.inlet_pressure_Pa, x.mass_flow_kg_s,
        x.speed_rpm, x.slip_factor, x.inlet_prewhirl_m_s,
        x.base_entropy_generation_J_kgK + passage_entropy
    };
    const turbocompgpu::NonlinearResidualScales scales{
        x.specific_work_scale_J_kg,
        x.specific_entropy_scale_J_kgK,
        x.pressure_scale_Pa
    };
    const turbocompgpu::NonlinearUnknowns initial{
        x.initial_temperature_K, x.initial_log_pressure
    };
    const auto context = turbocompgpu::make_fast_residual_context(
        mixture, geometry, operating, scales);

    turbocompgpu::NonlinearSolverOptions options;
    options.maximum_iterations = gpu_options.maximum_iterations;
    options.maximum_backtracks = gpu_options.maximum_backtracks;
    options.residual_tolerance = gpu_options.residual_tolerance;
    options.difference_step = gpu_options.difference_step;
    options.temperature_scale = gpu_options.temperature_scale_K;
    options.minimum_step = gpu_options.minimum_step;
    options.jacobian_singularity_tolerance =
        gpu_options.jacobian_singularity_tolerance;

    const auto solved = turbocompgpu::solve_fast_nonlinear_compressor(
        context, initial, options);
    if (!solved.converged() || !solved.solution)
        throw std::runtime_error("CPU impeller solve failed to converge");
    const auto& c = *solved.solution;
    const auto state2 = mixture.evaluate_Tp(
        c.unknowns.outlet_temperature, c.outlet_pressure);
    if (!state2) throw std::runtime_error("CPU impeller outlet EOS failed");

    const turbocompgpu::VanelessDiffuserGeometry diffuser_geometry{
        x.exit_diameter_m / 2.0, x.diffuser_outlet_radius_m,
        x.exit_width_m, 0.0
    };
    const turbocompgpu::VanelessDiffuserInlet diffuser_inlet{
        state2.value, c.outlet_meridional_velocity, c.outlet_whirl,
        x.mass_flow_kg_s
    };
    turbocompgpu::FrictionalDiffuserOptions diffuser_options;
    diffuser_options.radial_steps = x.diffuser_radial_steps;
    const auto d =
        turbocompgpu::evaluate_frictional_vaneless_diffuser_local_colebrook_analytical(
            mixture, mixture, diffuser_geometry, diffuser_inlet,
            viscosity_model, gpu_options.use_iso20765_viscosity,
            x.dynamic_viscosity_Pa_s, x.relative_roughness, diffuser_options);

    ImpellerBatchOutput out{};
    out.status = ImpellerBatchStatus::Converged;
    out.iterations = solved.iterations;
    out.outlet_temperature_K = c.unknowns.outlet_temperature;
    out.outlet_pressure_Pa = c.outlet_pressure;
    out.outlet_density_kg_m3 = state2.value.rho;
    out.meridional_velocity_m_s = c.outlet_meridional_velocity;
    out.whirl_velocity_m_s = c.outlet_whirl;
    out.euler_work_J_kg = c.euler_specific_work;
    out.shaft_torque_N_m = c.shaft_torque;
    out.residual_norm = c.residual_norm;
    out.darcy_friction_factor = friction.darcy_friction_factor;
    out.passage_entropy_generation_J_kgK = passage_entropy;
    out.diffuser_outlet_temperature_K = d.outlet_static.T;
    out.diffuser_outlet_pressure_Pa = d.outlet_static.p;
    out.diffuser_outlet_density_kg_m3 = d.outlet_static.rho;
    out.diffuser_outlet_meridional_velocity_m_s = d.outlet_radial_velocity;
    out.diffuser_outlet_whirl_m_s = d.outlet_tangential_velocity;
    out.diffuser_stagnation_enthalpy_residual_J_kg =
        d.stagnation_enthalpy_residual;
    out.diffuser_entropy_change_J_kgK = d.entropy_change;
    out.diffuser_integrated_entropy_J_kgK = d.integrated_entropy_production;
    out.diffuser_entropy_balance_residual_J_kgK = d.entropy_balance_residual;
    out.diffuser_angular_momentum_residual_m2_s =
        d.angular_momentum_balance_residual;
    out.diffuser_integrated_wall_torque_m2_s_per_kg =
        d.integrated_wall_torque_per_mass;
    out.diffuser_mass_residual_kg_s = d.mass_residual;
    return out;
}

std::vector<ImpellerBatchOutput> cpu_batch(
    const std::vector<ImpellerBatchInput>& inputs,
    const thermomodel::PengRobinsonMixture& mixture,
    const turbocompgpu::transport::Iso20765Parameters& viscosity_model,
    const ImpellerBatchOptions& options,
    int threads)
{
    std::vector<ImpellerBatchOutput> outputs(inputs.size());
    std::vector<std::string> errors(inputs.size());
    const auto solve_one = [&](std::size_t i) {
        try {
            outputs[i] = cpu_solve(inputs[i], mixture, viscosity_model, options);
        } catch (const std::exception& e) {
            errors[i] = e.what();
        }
    };
#ifdef TURBOCOMPGPU_BENCHMARK_OPENMP
#pragma omp parallel for schedule(static) num_threads(threads)
    for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(inputs.size()); ++i) {
        solve_one(static_cast<std::size_t>(i));
    }
#else
    const int worker_count = std::min<int>(threads, static_cast<int>(inputs.size()));
    if (worker_count <= 1) {
        for (std::size_t i=0;i<inputs.size();++i) solve_one(i);
    } else {
        std::vector<std::thread> workers;
        workers.reserve(static_cast<std::size_t>(worker_count));
        for (int worker=0;worker<worker_count;++worker) {
            const std::size_t begin=inputs.size()*static_cast<std::size_t>(worker)/worker_count;
            const std::size_t end=inputs.size()*static_cast<std::size_t>(worker+1)/worker_count;
            workers.emplace_back([&,begin,end]{
                for (std::size_t i=begin;i<end;++i) solve_one(i);
            });
        }
        for (auto& worker:workers) worker.join();
    }
#endif
    for (const auto& error : errors)
        if (!error.empty()) throw std::runtime_error(error);
    return outputs;
}

struct ParityStats {
    double max_tolerance_ratio = 0.0;
    double max_pressure_error_Pa = 0.0;
    double max_temperature_error_K = 0.0;
};

void compare_value(const char* name, double cpu, double gpu,
                   double atol, double rtol, ParityStats& stats)
{
    const double error = std::abs(cpu - gpu);
    const double tolerance = atol + rtol * std::max(std::abs(cpu), std::abs(gpu));
    if (!std::isfinite(cpu) || !std::isfinite(gpu) || error > tolerance) {
        throw std::runtime_error(std::string("batch CPU/GPU mismatch: ") + name +
            " cpu=" + std::to_string(cpu) + " gpu=" + std::to_string(gpu) +
            " abs_error=" + std::to_string(error));
    }
    stats.max_tolerance_ratio = std::max(stats.max_tolerance_ratio, error / tolerance);
    if (std::string(name).find("pressure") != std::string::npos)
        stats.max_pressure_error_Pa = std::max(stats.max_pressure_error_Pa, error);
    if (std::string(name).find("temperature") != std::string::npos)
        stats.max_temperature_error_K = std::max(stats.max_temperature_error_K, error);
}

ParityStats check_parity(const std::vector<ImpellerBatchOutput>& cpu,
                         const std::vector<ImpellerBatchOutput>& gpu)
{
    if (cpu.size() != gpu.size()) throw std::runtime_error("batch lengths differ");
    ParityStats stats;
    for (std::size_t i = 0; i < cpu.size(); ++i) {
        const auto& c = cpu[i];
        const auto& g = gpu[i];
        if (g.status != ImpellerBatchStatus::Converged ||
            !(g.euler_work_J_kg > 0.0 && g.whirl_velocity_m_s > 0.0))
            throw std::runtime_error("GPU batch point did not converge on compressor branch");
        compare_value("impeller pressure",c.outlet_pressure_Pa,g.outlet_pressure_Pa,2.0,2e-7,stats);
        compare_value("impeller temperature",c.outlet_temperature_K,g.outlet_temperature_K,2e-4,2e-7,stats);
        compare_value("impeller density",c.outlet_density_kg_m3,g.outlet_density_kg_m3,2e-5,2e-7,stats);
        compare_value("impeller meridional velocity",c.meridional_velocity_m_s,g.meridional_velocity_m_s,2e-4,2e-7,stats);
        compare_value("impeller whirl velocity",c.whirl_velocity_m_s,g.whirl_velocity_m_s,2e-4,2e-7,stats);
        compare_value("Euler work",c.euler_work_J_kg,g.euler_work_J_kg,2e-2,2e-7,stats);
        compare_value("shaft torque",c.shaft_torque_N_m,g.shaft_torque_N_m,2e-4,2e-7,stats);
        compare_value("Darcy friction factor",c.darcy_friction_factor,g.darcy_friction_factor,2e-13,2e-11,stats);
        compare_value("passage entropy generation",c.passage_entropy_generation_J_kgK,g.passage_entropy_generation_J_kgK,2e-13,2e-11,stats);
        compare_value("diffuser pressure",c.diffuser_outlet_pressure_Pa,g.diffuser_outlet_pressure_Pa,10.0,2e-6,stats);
        compare_value("diffuser temperature",c.diffuser_outlet_temperature_K,g.diffuser_outlet_temperature_K,2e-4,2e-6,stats);
        compare_value("diffuser density",c.diffuser_outlet_density_kg_m3,g.diffuser_outlet_density_kg_m3,2e-5,2e-6,stats);
        compare_value("diffuser meridional velocity",c.diffuser_outlet_meridional_velocity_m_s,g.diffuser_outlet_meridional_velocity_m_s,2e-4,2e-6,stats);
        compare_value("diffuser whirl velocity",c.diffuser_outlet_whirl_m_s,g.diffuser_outlet_whirl_m_s,2e-4,2e-6,stats);
        compare_value("diffuser energy residual",c.diffuser_stagnation_enthalpy_residual_J_kg,g.diffuser_stagnation_enthalpy_residual_J_kg,0.2,2e-5,stats);
        compare_value("diffuser entropy change",c.diffuser_entropy_change_J_kgK,g.diffuser_entropy_change_J_kgK,2e-6,2e-6,stats);
        compare_value("diffuser entropy production",c.diffuser_integrated_entropy_J_kgK,g.diffuser_integrated_entropy_J_kgK,2e-6,2e-6,stats);
        compare_value("diffuser entropy balance",c.diffuser_entropy_balance_residual_J_kgK,g.diffuser_entropy_balance_residual_J_kgK,2e-5,2e-5,stats);
        compare_value("diffuser torque integral",c.diffuser_integrated_wall_torque_m2_s_per_kg,g.diffuser_integrated_wall_torque_m2_s_per_kg,2e-6,2e-6,stats);
        compare_value("diffuser angular balance",c.diffuser_angular_momentum_residual_m2_s,g.diffuser_angular_momentum_residual_m2_s,2e-6,2e-5,stats);
        compare_value("diffuser mass balance",c.diffuser_mass_residual_kg_s,g.diffuser_mass_residual_kg_s,1e-9,2e-7,stats);
    }
    return stats;
}

std::vector<std::size_t> parse_sizes(const std::string& text)
{
    std::vector<std::size_t> sizes;
    std::stringstream stream(text);
    std::string item;
    while (std::getline(stream, item, ',')) {
        const auto n = std::stoull(item);
        if (n == 0) throw std::invalid_argument("batch sizes must be positive");
        sizes.push_back(static_cast<std::size_t>(n));
    }
    if (sizes.empty()) throw std::invalid_argument("no batch sizes supplied");
    std::sort(sizes.begin(), sizes.end());
    sizes.erase(std::unique(sizes.begin(), sizes.end()), sizes.end());
    return sizes;
}

template<class F>
double median_seconds(int repeats, F&& function)
{
    std::vector<double> samples;
    samples.reserve(static_cast<std::size_t>(repeats));
    for (int i = 0; i < repeats; ++i) {
        const auto start = Clock::now();
        function();
        samples.push_back(std::chrono::duration<double>(Clock::now()-start).count());
    }
    std::sort(samples.begin(), samples.end());
    return samples[samples.size()/2];
}

} // namespace

int main(int argc, char** argv)
{
    try {
        std::string sizes_text = "1,4,16,64,256";
        int repeats = 3;
#ifdef TURBOCOMPGPU_BENCHMARK_OPENMP
        int cpu_threads = std::max(1, omp_get_max_threads());
#else
        int cpu_threads = static_cast<int>(std::max(1u,std::thread::hardware_concurrency()));
#endif
        for (int i=1;i<argc;++i) {
            const std::string arg=argv[i];
            if (arg=="--sizes" && i+1<argc) sizes_text=argv[++i];
            else if (arg=="--repeats" && i+1<argc) repeats=std::stoi(argv[++i]);
            else if (arg=="--cpu-threads" && i+1<argc) cpu_threads=std::stoi(argv[++i]);
            else throw std::invalid_argument("usage: benchmark [--sizes 1,4,...] [--repeats N] [--cpu-threads N]");
        }
        if (repeats<1 || cpu_threads<1) throw std::invalid_argument("repeats and CPU threads must be positive");
        const auto sizes=parse_sizes(sizes_text);
        const auto mixture=turbocompgpu::benchmark_data::lin2026_natural_gas();
        const auto viscosity_model=
            turbocompgpu::transport::make_iso20765_parameters(mixture);
        const ImpellerBatchOptions options{};
        turbocompgpu::cuda::ImpellerBatchSolver gpu_solver(mixture,sizes.back());

        const auto warm_point=point(3.0);
        const auto cpu_warm=cpu_solve(warm_point,mixture,viscosity_model,options);
        const auto gpu_warm=gpu_solver.solve({warm_point},options);
        (void)check_parity({cpu_warm},gpu_warm);

        std::cout << "# Model: LIN-2026 screening mixture; PR EOS; ISO 20765-5 local viscosity; "
                     "inlet-Re impeller friction; stage-local Colebrook diffuser; "
                     "512-step constant-width diffuser\n"
                  << "# GPU buffers are retained at solver capacity; call time includes H2D/D2H copies and synchronization\n"
                  << "# CPU parallel path parallelizes independent operating points with OpenMP when available\n"
                  << "# CPU threads requested: " << cpu_threads << "; repeats: " << repeats << '\n'
                  << "batch_size,cpu_serial_s,cpu_serial_points_s,cpu_parallel_s,cpu_parallel_points_s,"
                     "gpu_end_to_end_s,gpu_points_s,gpu_vs_cpu_serial,gpu_vs_cpu_parallel,"
                     "max_parity_tolerance_ratio,max_abs_pressure_error_Pa,max_abs_temperature_error_K\n";

        for (const auto n : sizes) {
            std::vector<ImpellerBatchInput> inputs;
            inputs.reserve(n);
            for (std::size_t i=0;i<n;++i)
                inputs.push_back(point(std::array<double,4>{3.0,4.0,5.0,6.0}[i%4]));

            // Untimed warm-up at this size avoids including first-use effects.
            const auto cpu_serial_warm=cpu_batch(inputs,mixture,viscosity_model,options,1);
            const auto cpu_parallel_warm=cpu_batch(inputs,mixture,viscosity_model,options,cpu_threads);
            const auto gpu_batch_warm=gpu_solver.solve(inputs,options);
            const auto stats=check_parity(cpu_serial_warm,gpu_batch_warm);
            (void)check_parity(cpu_parallel_warm,gpu_batch_warm);

            const double serial_s=median_seconds(repeats,[&]{(void)cpu_batch(inputs,mixture,viscosity_model,options,1);});
            const double parallel_s=median_seconds(repeats,[&]{(void)cpu_batch(inputs,mixture,viscosity_model,options,cpu_threads);});
            const double gpu_s=median_seconds(repeats,[&]{(void)gpu_solver.solve(inputs,options);});
            const double count=static_cast<double>(n);
            std::cout << n << ',' << serial_s << ',' << count/serial_s << ','
                      << parallel_s << ',' << count/parallel_s << ','
                      << gpu_s << ',' << count/gpu_s << ','
                      << serial_s/gpu_s << ',' << parallel_s/gpu_s << ','
                      << stats.max_tolerance_ratio << ','
                      << stats.max_pressure_error_Pa << ','
                      << stats.max_temperature_error_K << '\n';
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "GPU impeller throughput benchmark failed: " << e.what() << '\n';
        return 1;
    }
}
