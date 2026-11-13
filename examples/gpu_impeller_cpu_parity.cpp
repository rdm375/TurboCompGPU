// SPDX-License-Identifier: GPL-3.0-or-later
#include <cmath>
// ThermoGPU's CUDA PR-mixture header uses ::isfinite in functions that are
// also parsed by host C++. Make the standard overload visible at global scope
// for this C++ parity driver (the CUDA translation units provide it directly).
using std::isfinite;

#include <turbocompgpu/fast_nonlinear_solver.hpp>
#include <turbocompgpu/gpu_impeller_batch_cuda.hpp>
#include <turbocompgpu/vaneless_diffuser_friction_analytical.hpp>
#include "lin2026_natural_gas.hpp"

#include <friction/colebrook.hpp>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require_close(const char* name, double cpu, double gpu,
                   double atol, double rtol)
{
    const double error = std::abs(cpu-gpu);
    if (!std::isfinite(cpu) || !std::isfinite(gpu) ||
        error > atol+rtol*std::max(std::abs(cpu),std::abs(gpu))) {
        throw std::runtime_error(std::string("CPU/GPU parity failed for ")+name+
            ": cpu="+std::to_string(cpu)+", gpu="+std::to_string(gpu)+
            ", abs_error="+std::to_string(error));
    }
}

} // namespace

int main()
{
    try {
        using namespace turbocompgpu;
        using namespace turbocompgpu::cuda;

        const auto mixture=benchmark_data::lin2026_natural_gas();
        const auto viscosity_model=transport::make_iso20765_parameters(mixture);
        const std::vector<double> speeds_rpm{10000.0,12000.0,14000.0,16000.0,18000.0};
        const std::vector<double> flows_kg_s{3.0,4.0,5.0,6.0};
        std::vector<ImpellerBatchInput> points;
        points.reserve(speeds_rpm.size()*flows_kg_s.size());
        for (double speed : speeds_rpm) {
            for (double flow : flows_kg_s) {
                points.push_back({
                    293.15,5.0e5,flow,speed,0.90,0.0,0.0,
                    0.030,0.150,0.03,
                    0.200,0.030,0.0,0.090,50.0,
                    5.0e4,100.0,5.0e5,
                    0.15,0.020,1.1e-5,0.0,
                    310.0,0.20,
                    0.320,512
                });
            }
        }
        ImpellerBatchSolver gpu_solver(mixture,points.size());
        const ImpellerBatchOptions gpu_options{};
        const auto gpu=gpu_solver.solve(points,gpu_options);

        NonlinearSolverOptions cpu_options;
        cpu_options.maximum_iterations=gpu_options.maximum_iterations;
        cpu_options.maximum_backtracks=gpu_options.maximum_backtracks;
        cpu_options.residual_tolerance=gpu_options.residual_tolerance;
        cpu_options.difference_step=gpu_options.difference_step;
        cpu_options.temperature_scale=gpu_options.temperature_scale_K;
        cpu_options.minimum_step=gpu_options.minimum_step;
        cpu_options.jacobian_singularity_tolerance=
            gpu_options.jacobian_singularity_tolerance;

        std::cout << std::setprecision(12)
                  << "speed_rpm,flow_kg_s,cpu_T_K,gpu_T_K,abs_dT_K,cpu_p_Pa,gpu_p_Pa,abs_dp_Pa,"
                     "cpu_rho,gpu_rho,cpu_vm,gpu_vm,cpu_whirl,gpu_whirl,"
                     "cpu_work_J_kg,gpu_work_J_kg,cpu_torque_Nm,gpu_torque_Nm,"
                     "cpu_f_D,gpu_f_D,cpu_passage_entropy_J_kgK,gpu_passage_entropy_J_kgK,"
                     "residual_norm_gpu,"
                     "cpu_diffuser_p3_Pa,gpu_diffuser_p3_Pa,"
                     "cpu_diffuser_T3_K,gpu_diffuser_T3_K,"
                     "cpu_diffuser_vm3_m_s,gpu_diffuser_vm3_m_s,"
                     "cpu_diffuser_vtheta3_m_s,gpu_diffuser_vtheta3_m_s,"
                     "cpu_diffuser_dh0_residual_J_kg,gpu_diffuser_dh0_residual_J_kg,"
                     "cpu_diffuser_entropy_change_J_kgK,gpu_diffuser_entropy_change_J_kgK,"
                     "cpu_diffuser_entropy_production_J_kgK,gpu_diffuser_entropy_production_J_kgK,"
                     "cpu_diffuser_entropy_residual_J_kgK,gpu_diffuser_entropy_residual_J_kgK,"
                     "cpu_diffuser_wall_torque_integral_m2_s_kg,gpu_diffuser_wall_torque_integral_m2_s_kg,"
                     "cpu_diffuser_angular_residual_m2_s,gpu_diffuser_angular_residual_m2_s,"
                     "cpu_diffuser_mass_residual_kg_s,gpu_diffuser_mass_residual_kg_s\n";

        for (std::size_t i=0;i<points.size();++i) {
            const auto& x=points[i];
            try {
                const auto& g=gpu[i];
                if (g.status!=ImpellerBatchStatus::Converged ||
                    !(g.whirl_velocity_m_s>0.0 && g.euler_work_J_kg>0.0))
                    throw std::runtime_error("GPU impeller did not converge to a compressor branch");

                constexpr double pi=3.14159265358979323846;
                const double inlet_area=(pi/4.0)*
                    (x.inlet_shroud_diameter_m*x.inlet_shroud_diameter_m-
                     x.inlet_hub_diameter_m*x.inlet_hub_diameter_m)*(1.0-x.inlet_blockage);
                const auto inlet=mixture.evaluate_Tp(x.inlet_temperature_K,x.inlet_pressure_Pa);
                if (!inlet) throw std::runtime_error("CPU inlet EOS evaluation failed");
                const double vm1=x.mass_flow_kg_s/(inlet.value.rho*inlet_area);
                const double relative_speed=std::hypot(vm1,x.inlet_prewhirl_m_s);
                const double inlet_mu=gpu_options.use_iso20765_viscosity
                    ? transport::iso20765_5_dynamic_viscosity_Pa_s(
                          viscosity_model,inlet.value.T,inlet.value.p,inlet.value.rho)
                    : x.dynamic_viscosity_Pa_s;
                if (!(std::isfinite(inlet_mu) && inlet_mu>0.0))
                    throw std::runtime_error("CPU inlet viscosity evaluation failed");
                const double reynolds=inlet.value.rho*relative_speed*x.hydraulic_diameter_m/
                    inlet_mu;
                const auto f_cpu=friction::colebrook_turbulent(reynolds,x.relative_roughness);
                const double passage_friction_work=f_cpu.darcy_friction_factor*
                    (x.passage_length_m/x.hydraulic_diameter_m)*
                    0.5*relative_speed*relative_speed;
                const double passage_entropy=passage_friction_work/x.inlet_temperature_K;
                require_close("Darcy friction factor",f_cpu.darcy_friction_factor,
                              g.darcy_friction_factor,2.0e-13,2.0e-11);
                require_close("impeller passage entropy generation",passage_entropy,
                              g.passage_entropy_generation_J_kgK,2.0e-13,2.0e-11);

                const ImpellerSystemGeometry geometry{
                    {x.inlet_hub_diameter_m,x.inlet_shroud_diameter_m,x.inlet_blockage},
                    {x.exit_diameter_m,x.exit_width_m,x.exit_blockage},
                    {x.exit_diameter_m,x.exit_blade_angle_deg,x.inlet_impeller_diameter_m}
                };
                const NonlinearOperatingPoint operating{
                    x.inlet_temperature_K,x.inlet_pressure_Pa,x.mass_flow_kg_s,
                    x.speed_rpm,x.slip_factor,x.inlet_prewhirl_m_s,
                    x.base_entropy_generation_J_kgK+passage_entropy
                };
                const NonlinearResidualScales scales{
                    x.specific_work_scale_J_kg,
                    x.specific_entropy_scale_J_kgK,
                    x.pressure_scale_Pa
                };
                const NonlinearUnknowns initial{
                    x.initial_temperature_K,x.initial_log_pressure
                };
                const auto context=make_fast_residual_context(
                    mixture,geometry,operating,scales);
                const auto cpu=solve_fast_nonlinear_compressor(
                    context,initial,cpu_options);
                if (!cpu.converged() || !cpu.solution)
                    throw std::runtime_error("CPU fast residual reference did not converge");
                const auto& c=*cpu.solution;
                const auto cstate=mixture.evaluate_Tp(
                    c.unknowns.outlet_temperature,c.outlet_pressure);
                if (!cstate) throw std::runtime_error("CPU outlet EOS evaluation failed");

                const VanelessDiffuserGeometry diffuser_geometry{
                    x.exit_diameter_m/2.0,x.diffuser_outlet_radius_m,
                    x.exit_width_m,0.0
                };
                const VanelessDiffuserInlet diffuser_inlet{
                    cstate.value,c.outlet_meridional_velocity,c.outlet_whirl,
                    x.mass_flow_kg_s
                };
                FrictionalDiffuserOptions diffuser_options;
                diffuser_options.radial_steps=x.diffuser_radial_steps;
                const auto cpu_diffuser=
                    evaluate_frictional_vaneless_diffuser_local_colebrook_analytical(
                        mixture,mixture,diffuser_geometry,diffuser_inlet,
                        viscosity_model,gpu_options.use_iso20765_viscosity,
                        x.dynamic_viscosity_Pa_s,x.relative_roughness,diffuser_options);

                require_close("outlet temperature",c.unknowns.outlet_temperature,
                              g.outlet_temperature_K,2.0e-4,2.0e-7);
                require_close("outlet pressure",c.outlet_pressure,
                              g.outlet_pressure_Pa,2.0,2.0e-7);
                require_close("outlet density",cstate.value.rho,
                              g.outlet_density_kg_m3,2.0e-5,2.0e-7);
                require_close("exit meridional velocity",c.outlet_meridional_velocity,
                              g.meridional_velocity_m_s,2.0e-4,2.0e-7);
                require_close("exit whirl velocity",c.outlet_whirl,
                              g.whirl_velocity_m_s,2.0e-4,2.0e-7);
                require_close("Euler work",c.euler_specific_work,
                              g.euler_work_J_kg,2.0e-2,2.0e-7);
                require_close("shaft torque",c.shaft_torque,
                              g.shaft_torque_N_m,2.0e-4,2.0e-7);
                require_close("diffuser outlet pressure",cpu_diffuser.outlet_static.p,
                              g.diffuser_outlet_pressure_Pa,10.0,2.0e-6);
                require_close("diffuser outlet temperature",cpu_diffuser.outlet_static.T,
                              g.diffuser_outlet_temperature_K,2.0e-4,2.0e-6);
                require_close("diffuser outlet density",cpu_diffuser.outlet_static.rho,
                              g.diffuser_outlet_density_kg_m3,2.0e-5,2.0e-6);
                require_close("diffuser outlet meridional velocity",
                              cpu_diffuser.outlet_radial_velocity,
                              g.diffuser_outlet_meridional_velocity_m_s,2.0e-4,2.0e-6);
                require_close("diffuser outlet whirl velocity",
                              cpu_diffuser.outlet_tangential_velocity,
                              g.diffuser_outlet_whirl_m_s,2.0e-4,2.0e-6);
                require_close("diffuser stagnation enthalpy residual",
                              cpu_diffuser.stagnation_enthalpy_residual,
                              g.diffuser_stagnation_enthalpy_residual_J_kg,0.2,2.0e-5);
                require_close("diffuser entropy change",
                              cpu_diffuser.entropy_change,
                              g.diffuser_entropy_change_J_kgK,2.0e-6,2.0e-6);
                require_close("diffuser integrated entropy production",
                              cpu_diffuser.integrated_entropy_production,
                              g.diffuser_integrated_entropy_J_kgK,2.0e-6,2.0e-6);
                require_close("diffuser entropy balance residual",
                              cpu_diffuser.entropy_balance_residual,
                              g.diffuser_entropy_balance_residual_J_kgK,2.0e-5,2.0e-5);
                require_close("diffuser integrated wall torque per mass",
                              cpu_diffuser.integrated_wall_torque_per_mass,
                              g.diffuser_integrated_wall_torque_m2_s_per_kg,2.0e-6,2.0e-6);
                require_close("diffuser angular momentum residual",
                              cpu_diffuser.angular_momentum_balance_residual,
                              g.diffuser_angular_momentum_residual_m2_s,2.0e-6,2.0e-5);
                require_close("diffuser mass residual",cpu_diffuser.mass_residual,
                              g.diffuser_mass_residual_kg_s,1.0e-9,2.0e-7);

                std::cout << x.speed_rpm << ',' << x.mass_flow_kg_s << ','
                          << c.unknowns.outlet_temperature << ',' << g.outlet_temperature_K << ','
                          << std::abs(c.unknowns.outlet_temperature-g.outlet_temperature_K) << ','
                          << c.outlet_pressure << ',' << g.outlet_pressure_Pa << ','
                          << std::abs(c.outlet_pressure-g.outlet_pressure_Pa) << ','
                          << cstate.value.rho << ',' << g.outlet_density_kg_m3 << ','
                          << c.outlet_meridional_velocity << ',' << g.meridional_velocity_m_s << ','
                          << c.outlet_whirl << ',' << g.whirl_velocity_m_s << ','
                          << c.euler_specific_work << ',' << g.euler_work_J_kg << ','
                          << c.shaft_torque << ',' << g.shaft_torque_N_m << ','
                          << f_cpu.darcy_friction_factor << ',' << g.darcy_friction_factor << ','
                          << passage_entropy << ',' << g.passage_entropy_generation_J_kgK << ','
                          << g.residual_norm << ','
                          << cpu_diffuser.outlet_static.p << ',' << g.diffuser_outlet_pressure_Pa << ','
                          << cpu_diffuser.outlet_static.T << ',' << g.diffuser_outlet_temperature_K << ','
                          << cpu_diffuser.outlet_radial_velocity << ','
                          << g.diffuser_outlet_meridional_velocity_m_s << ','
                          << cpu_diffuser.outlet_tangential_velocity << ','
                          << g.diffuser_outlet_whirl_m_s << ','
                          << cpu_diffuser.stagnation_enthalpy_residual << ','
                          << g.diffuser_stagnation_enthalpy_residual_J_kg << ','
                          << cpu_diffuser.entropy_change << ','
                          << g.diffuser_entropy_change_J_kgK << ','
                          << cpu_diffuser.integrated_entropy_production << ','
                          << g.diffuser_integrated_entropy_J_kgK << ','
                          << cpu_diffuser.entropy_balance_residual << ','
                          << g.diffuser_entropy_balance_residual_J_kgK << ','
                          << cpu_diffuser.integrated_wall_torque_per_mass << ','
                          << g.diffuser_integrated_wall_torque_m2_s_per_kg << ','
                          << cpu_diffuser.angular_momentum_balance_residual << ','
                          << g.diffuser_angular_momentum_residual_m2_s << ','
                          << cpu_diffuser.mass_residual << ','
                          << g.diffuser_mass_residual_kg_s << '\n';
            } catch (const std::exception& error) {
                throw std::runtime_error(
                    "Parity failed at " + std::to_string(x.speed_rpm) +
                    " rpm, " + std::to_string(x.mass_flow_kg_s) +
                    " kg/s: " + error.what());
            }
        }
        std::cerr << "CPU/GPU impeller, Colebrook, and frictional diffuser parity passed for "
                  << points.size() << " points across " << speeds_rpm.size()
                  << " speeds and " << flows_kg_s.size() << " flow rates.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
