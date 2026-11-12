// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once

#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/detail/pr_mixture_state_cuda.hpp>
#include <thermogpu/thermo_mixture_batch_cuda.hpp>
#include <turbocompgpu/iso20765_5_viscosity.hpp>

#include <cstddef>
#include <memory>
#include <vector>

namespace turbocompgpu::cuda {

// One independent impeller operating point. SI units throughout.
// Impeller passage friction is currently a lumped Darcy-Weisbach screening
// closure using inlet relative speed and inlet Re. The diffuser instead
// recomputes local PR density, velocity, ISO 20765-5 viscosity, Re, and
// Colebrook f at each RK4 stage. Relative roughness is held constant.
// Other empirical loss terms are supplied as base_entropy_generation.
struct ImpellerBatchInput {
    double inlet_temperature_K;
    double inlet_pressure_Pa;
    double mass_flow_kg_s;
    double speed_rpm;
    double slip_factor;
    double inlet_prewhirl_m_s;
    double base_entropy_generation_J_kgK;

    double inlet_hub_diameter_m;
    double inlet_shroud_diameter_m;
    double inlet_blockage;
    double exit_diameter_m;
    double exit_width_m;
    double exit_blockage;
    double inlet_impeller_diameter_m;
    double exit_blade_angle_deg;

    double specific_work_scale_J_kg;
    double specific_entropy_scale_J_kgK;
    double pressure_scale_Pa;

    // Darcy-Weisbach passage-friction screening inputs.
    double passage_length_m;
    double hydraulic_diameter_m;
    // Fallback viscosity used only when ImpellerBatchOptions disables ISO
    // 20765-5; retained as a scalar input for API compatibility.
    double dynamic_viscosity_Pa_s;
    double relative_roughness;

    double initial_temperature_K;
    double initial_log_pressure;

    // Constant-width, zero-blockage vaneless diffuser exit radius. Its inlet
    // radius is fixed by the impeller exit diameter/2 and its width is the
    // impeller exit width, matching the currently validated CPU stage model.
    double diffuser_outlet_radius_m;
    int diffuser_radial_steps;
};

enum class ImpellerBatchStatus : int {
    Converged = 0,
    InvalidInput = 1,
    InletThermoFailed = 2,
    OutletThermoFailed = 3,
    FrictionFailed = 4,
    JacobianFailed = 5,
    SingularJacobian = 6,
    LineSearchFailed = 7,
    MaximumIterations = 8
};

struct ImpellerBatchOutput {
    ImpellerBatchStatus status;
    int iterations;
    double outlet_temperature_K;
    double outlet_pressure_Pa;
    double outlet_density_kg_m3;
    double meridional_velocity_m_s;
    double whirl_velocity_m_s;
    double euler_work_J_kg;
    double shaft_torque_N_m;
    double energy_residual_J_kg;
    double entropy_residual_J_kgK;
    double residual_norm;
    double darcy_friction_factor;
    double passage_entropy_generation_J_kgK;
    double diffuser_outlet_temperature_K;
    double diffuser_outlet_pressure_Pa;
    double diffuser_outlet_density_kg_m3;
    double diffuser_outlet_meridional_velocity_m_s;
    double diffuser_outlet_whirl_m_s;
    double diffuser_stagnation_enthalpy_residual_J_kg;
    double diffuser_entropy_change_J_kgK;
    double diffuser_integrated_entropy_J_kgK;
    double diffuser_entropy_balance_residual_J_kgK;
    double diffuser_angular_momentum_residual_m2_s;
    double diffuser_integrated_wall_torque_m2_s_per_kg;
    double diffuser_mass_residual_kg_s;
};

struct ImpellerBatchOptions {
    int maximum_iterations = 40;
    int maximum_backtracks = 24;
    double residual_tolerance = 1.0e-10;
    double difference_step = 1.0e-5;
    double temperature_scale_K = 300.0;
    double minimum_step = 1.0e-8;
    double jacobian_singularity_tolerance = 1.0e-12;
    // Use ISO 20765-5:2022 gas viscosity from local T, p, PR density, and
    // fixed composition. If false, use ImpellerBatchInput::dynamic_viscosity_Pa_s.
    bool use_iso20765_viscosity = true;
};

class ImpellerBatchSolver {
public:
    ImpellerBatchSolver(const thermomodel::PengRobinsonMixture& mixture,
                        std::size_t capacity);
    ~ImpellerBatchSolver();
    ImpellerBatchSolver(const ImpellerBatchSolver&) = delete;
    ImpellerBatchSolver& operator=(const ImpellerBatchSolver&) = delete;
    ImpellerBatchSolver(ImpellerBatchSolver&&) noexcept;
    ImpellerBatchSolver& operator=(ImpellerBatchSolver&&) noexcept;

    // This backend covers the implicit impeller equations and the constant-
    // width frictional vaneless diffuser. Inputs/outputs are copied and the
    // call synchronizes before return. Device input/output storage is retained
    // at capacity between calls; concurrent calls on one solver serialize.
    // Volute and mechanical losses are omitted.
    [[nodiscard]] std::vector<ImpellerBatchOutput> solve(
        const std::vector<ImpellerBatchInput>& inputs,
        const ImpellerBatchOptions& options = {}) const;
    [[nodiscard]] std::size_t capacity() const noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace turbocompgpu::cuda
