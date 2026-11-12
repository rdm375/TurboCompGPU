// SPDX-License-Identifier: GPL-3.0-or-later
#include <turbocompgpu/gpu_impeller_batch_cuda.hpp>
#include "lin2026_natural_gas.hpp"

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

int main()
{
    try {
        const auto mixture = turbocompgpu::benchmark_data::lin2026_natural_gas();
        turbocompgpu::cuda::ImpellerBatchSolver solver(mixture, 4);

        std::vector<turbocompgpu::cuda::ImpellerBatchInput> points;
        for (double flow : {3.0, 4.0, 5.0, 6.0}) {
            points.push_back({
                293.15, 5.0e5, flow, 14000.0, 0.90, 0.0, 0.0,
                0.030, 0.150, 0.03,
                0.200, 0.030, 0.0, 0.090, 50.0,
                5.0e4, 100.0, 5.0e5,
                0.15, 0.020, 1.1e-5, 0.0,
                310.0, 0.20,
                0.320, 512
            });
        }
        const auto result = solver.solve(points);
        std::cout << std::setprecision(10)
                  << "flow_kg_s,status,iterations,p2_Pa,T2_K,vm2_m_s,Ctheta2_m_s,"
                     "Euler_work_J_kg,torque_N_m,Darcy_f,friction_entropy_J_kgK,residual_norm,"
                     "diffuser_p3_Pa,diffuser_T3_K,diffuser_vm3_m_s,diffuser_vtheta3_m_s,"
                     "diffuser_dh0_residual_J_kg,diffuser_entropy_balance_residual_J_kgK,"
                     "diffuser_angular_balance_residual_m2_s\n";
        for (std::size_t i=0;i<result.size();++i) {
            const auto& r=result[i];
            std::cout << points[i].mass_flow_kg_s << ','
                      << static_cast<int>(r.status) << ',' << r.iterations << ','
                      << r.outlet_pressure_Pa << ',' << r.outlet_temperature_K << ','
                      << r.meridional_velocity_m_s << ',' << r.whirl_velocity_m_s << ','
                      << r.euler_work_J_kg << ',' << r.shaft_torque_N_m << ','
                      << r.darcy_friction_factor << ','
                      << r.passage_entropy_generation_J_kgK << ',' << r.residual_norm << ','
                      << r.diffuser_outlet_pressure_Pa << ','
                      << r.diffuser_outlet_temperature_K << ','
                      << r.diffuser_outlet_meridional_velocity_m_s << ','
                      << r.diffuser_outlet_whirl_m_s << ','
                      << r.diffuser_stagnation_enthalpy_residual_J_kg << ','
                      << r.diffuser_entropy_balance_residual_J_kgK << ','
                      << r.diffuser_angular_momentum_residual_m2_s << '\n';
            if (r.status != turbocompgpu::cuda::ImpellerBatchStatus::Converged)
                throw std::runtime_error("A GPU impeller point did not converge");
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "GPU impeller batch failed: " << e.what() << '\n';
        return 1;
    }
}
