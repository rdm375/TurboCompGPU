cmake_minimum_required(VERSION 3.20)

project(TurboCompGPU
    VERSION 0.1.0
    DESCRIPTION "Physics-based centrifugal compressor modeling"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CTest)

find_package(ThermoGPU 2.0 REQUIRED CONFIG)

add_executable(turbocompgpu_thermo_smoke
    examples/thermo_smoke.cpp
)

target_link_libraries(turbocompgpu_thermo_smoke
    PRIVATE thermogpu::thermo
)

if(BUILD_TESTING)
    add_test(
        NAME thermogpu_integration
        COMMAND turbocompgpu_thermo_smoke
    )
endif()

add_library(turbocompgpu INTERFACE)

target_include_directories(turbocompgpu
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(turbocompgpu
    INTERFACE thermogpu::thermo
)

if(BUILD_TESTING)
    add_executable(turbocompgpu_compressor_test
        tests/compressor_ideal_gas.cpp
    )

    target_link_libraries(turbocompgpu_compressor_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME compressor_ideal_gas
        COMMAND turbocompgpu_compressor_test
    )
endif()

if(BUILD_TESTING)
    add_executable(turbocompgpu_peng_robinson_test
        tests/compressor_peng_robinson.cpp
    )

    target_link_libraries(turbocompgpu_peng_robinson_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME compressor_peng_robinson
        COMMAND turbocompgpu_peng_robinson_test
    )
endif()

if(BUILD_TESTING)
    add_executable(turbocompgpu_pr_derivatives_test
        tests/peng_robinson_derivatives.cpp
    )

    target_link_libraries(turbocompgpu_pr_derivatives_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME peng_robinson_derivatives
        COMMAND turbocompgpu_pr_derivatives_test
    )
endif()

if(BUILD_TESTING)
    add_executable(turbocompgpu_failure_test
        tests/compressor_failures.cpp
    )

    target_link_libraries(turbocompgpu_failure_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME compressor_failures
        COMMAND turbocompgpu_failure_test
    )
endif()

if(BUILD_TESTING)
    add_executable(turbocompgpu_reference_shift_test
        tests/compressor_reference_shift.cpp
    )

    target_link_libraries(turbocompgpu_reference_shift_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME compressor_reference_shift
        COMMAND turbocompgpu_reference_shift_test
    )
endif()

if(BUILD_TESTING)
    add_executable(turbocompgpu_impeller_test
        tests/impeller.cpp
    )

    target_link_libraries(turbocompgpu_impeller_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME impeller_aerodynamics
        COMMAND turbocompgpu_impeller_test
    )
endif()

if(BUILD_TESTING)
    add_executable(turbocompgpu_prewhirl_test
        tests/impeller_prewhirl.cpp
    )

    target_link_libraries(turbocompgpu_prewhirl_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME impeller_prewhirl
        COMMAND turbocompgpu_prewhirl_test
    )
endif()

if(BUILD_TESTING)
    add_executable(turbocompgpu_continuity_test
        tests/continuity.cpp
    )

    target_link_libraries(turbocompgpu_continuity_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME impeller_continuity
        COMMAND turbocompgpu_continuity_test
    )
endif()

if(BUILD_TESTING)
    add_executable(turbocompgpu_station_test
        tests/station.cpp
    )

    target_link_libraries(turbocompgpu_station_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME station_thermodynamics
        COMMAND turbocompgpu_station_test
    )
endif()

if(BUILD_TESTING)
    add_executable(turbocompgpu_station_reference_shift_test
        tests/station_reference_shift.cpp
    )

    target_link_libraries(
        turbocompgpu_station_reference_shift_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME station_reference_shift
        COMMAND turbocompgpu_station_reference_shift_test
    )
endif()

if(BUILD_TESTING)
    add_executable(
        turbocompgpu_conservation_test
        tests/conservation.cpp
    )

    target_link_libraries(
        turbocompgpu_conservation_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME impeller_conservation
        COMMAND turbocompgpu_conservation_test
    )
endif()

if(BUILD_TESTING)
    add_executable(
        turbocompgpu_continuity_independent_test
        tests/continuity_independent.cpp
    )

    target_link_libraries(
        turbocompgpu_continuity_independent_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME impeller_continuity_independent
        COMMAND turbocompgpu_continuity_independent_test
    )
endif()

if(BUILD_TESTING)
    add_executable(
        turbocompgpu_impeller_system_test
        tests/impeller_system.cpp
    )

    target_link_libraries(
        turbocompgpu_impeller_system_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME impeller_system
        COMMAND turbocompgpu_impeller_system_test
    )
endif()

if(BUILD_TESTING)
    add_executable(
        turbocompgpu_impeller_manufactured_test
        tests/impeller_system_manufactured.cpp
    )

    target_link_libraries(
        turbocompgpu_impeller_manufactured_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME impeller_manufactured
        COMMAND turbocompgpu_impeller_manufactured_test
    )
endif()

if(BUILD_TESTING)
    add_executable(
        turbocompgpu_loss_ideal_gas_test
        tests/loss_ideal_gas.cpp
    )

    target_link_libraries(
        turbocompgpu_loss_ideal_gas_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME loss_ideal_gas
        COMMAND turbocompgpu_loss_ideal_gas_test
    )
endif()

if(BUILD_TESTING)
    add_executable(
        turbocompgpu_impeller_system_loss_test
        tests/impeller_system_loss.cpp
    )

    target_link_libraries(
        turbocompgpu_impeller_system_loss_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME impeller_system_loss
        COMMAND turbocompgpu_impeller_system_loss_test
    )
endif()

if(BUILD_TESTING)
    add_executable(
        turbocompgpu_five_residual_test
        tests/impeller_system_five_residual.cpp
    )

    target_link_libraries(
        turbocompgpu_five_residual_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME impeller_five_residual
        COMMAND turbocompgpu_five_residual_test
    )
endif()

if(BUILD_TESTING)
    add_executable(
        turbocompgpu_nonlinear_residual_test
        tests/nonlinear_residual.cpp
    )

    target_link_libraries(
        turbocompgpu_nonlinear_residual_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME nonlinear_residual
        COMMAND turbocompgpu_nonlinear_residual_test
    )
endif()

if(BUILD_TESTING)
    add_executable(
        turbocompgpu_nonlinear_solver_test
        tests/nonlinear_solver.cpp
    )

    target_link_libraries(
        turbocompgpu_nonlinear_solver_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME nonlinear_solver
        COMMAND turbocompgpu_nonlinear_solver_test
    )
endif()

if(BUILD_TESTING)
    add_executable(
        turbocompgpu_nonlinear_solver_extended_test
        tests/nonlinear_solver_extended.cpp
    )

    target_link_libraries(
        turbocompgpu_nonlinear_solver_extended_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME nonlinear_solver_extended
        COMMAND turbocompgpu_nonlinear_solver_extended_test
    )
endif()

if(BUILD_TESTING)
    add_executable(
        turbocompgpu_nonlinear_evaluation_profile_test
        tests/nonlinear_evaluation_profile.cpp
    )

    target_link_libraries(
        turbocompgpu_nonlinear_evaluation_profile_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME nonlinear_evaluation_profile
        COMMAND turbocompgpu_nonlinear_evaluation_profile_test
    )
endif()

if(BUILD_TESTING)
    add_executable(
        turbocompgpu_fast_nonlinear_solver_test
        tests/fast_nonlinear_solver.cpp
    )

    target_link_libraries(
        turbocompgpu_fast_nonlinear_solver_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME fast_nonlinear_solver
        COMMAND turbocompgpu_fast_nonlinear_solver_test
    )
endif()

if(BUILD_TESTING)
    add_executable(
        turbocompgpu_nonlinear_solver_type_safety_test
        tests/nonlinear_solver_type_safety.cpp
    )

    target_link_libraries(
        turbocompgpu_nonlinear_solver_type_safety_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME nonlinear_solver_type_safety
        COMMAND turbocompgpu_nonlinear_solver_type_safety_test
    )
endif()

if(BUILD_TESTING)
    add_executable(
        turbocompgpu_fast_nonlinear_reference_invariance_test
        tests/fast_nonlinear_reference_invariance.cpp
    )

    target_link_libraries(
        turbocompgpu_fast_nonlinear_reference_invariance_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME fast_nonlinear_reference_invariance
        COMMAND turbocompgpu_fast_nonlinear_reference_invariance_test
    )
endif()

add_executable(
    turbocompgpu_convergence_basin
    benchmarks/convergence_basin.cpp
)

target_link_libraries(
    turbocompgpu_convergence_basin
    PRIVATE turbocompgpu
)

add_executable(
    turbocompgpu_nonlinear_solver_failure_paths_test
    tests/nonlinear_solver_failure_paths.cpp
)

target_link_libraries(
    turbocompgpu_nonlinear_solver_failure_paths_test
    PRIVATE turbocompgpu
)

add_test(
    NAME nonlinear_solver_failure_paths
    COMMAND turbocompgpu_nonlinear_solver_failure_paths_test
)

# M1.4d-3.5: Wall-clock performance benchmark.
# Not registered with CTest.
add_executable(
    turbocompgpu_nonlinear_wall_clock
    benchmarks/nonlinear_wall_clock.cpp
)

target_link_libraries(
    turbocompgpu_nonlinear_wall_clock
    PRIVATE turbocompgpu
)

# M1.5 operating-envelope infrastructure

add_executable(turbocompgpu_operating_envelope
    benchmarks/operating_envelope.cpp
)
target_link_libraries(turbocompgpu_operating_envelope
    PRIVATE turbocompgpu
)

if(BUILD_TESTING)
    add_executable(turbocompgpu_operating_point_test
        tests/operating_point.cpp
    )
    target_link_libraries(turbocompgpu_operating_point_test
        PRIVATE turbocompgpu
    )
    add_test(
        NAME operating_point
        COMMAND turbocompgpu_operating_point_test
    )
endif()

# M1.6a: Compressor-map coordinate conventions.
if(BUILD_TESTING)
    add_executable(turbocompgpu_map_coordinates_test
        tests/map_coordinates.cpp
    )

    target_link_libraries(turbocompgpu_map_coordinates_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME map_coordinates
        COMMAND turbocompgpu_map_coordinates_test
    )
endif()

# M1.6a: EOS-consistent compressor-map postprocessing.
if(BUILD_TESTING)
    add_executable(turbocompgpu_compressor_map_test
        tests/compressor_map.cpp
    )

    target_link_libraries(turbocompgpu_compressor_map_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME compressor_map
        COMMAND turbocompgpu_compressor_map_test
    )
endif()

# M1.6a: Corrected-coordinate compressor-map envelope.
add_executable(turbocompgpu_compressor_map_envelope
    benchmarks/compressor_map_envelope.cpp
)

target_link_libraries(turbocompgpu_compressor_map_envelope
    PRIVATE turbocompgpu
)

# M1.6b.1: Controlled thermodynamic-model comparison.
add_executable(
    turbocompgpu_thermodynamic_comparison
    benchmarks/thermodynamic_comparison.cpp
)

target_link_libraries(
    turbocompgpu_thermodynamic_comparison
    PRIVATE turbocompgpu
)

# M1.6b.1: Matched-methane ideal-gas limiting behavior.
if(BUILD_TESTING)
    add_executable(
        turbocompgpu_thermodynamic_limit_test
        tests/thermodynamic_limit.cpp
    )

    target_link_libraries(
        turbocompgpu_thermodynamic_limit_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME thermodynamic_limit
        COMMAND turbocompgpu_thermodynamic_limit_test
    )
endif()

# M1.6b.2: Independent aerodynamic incidence-loss kernel.
if(BUILD_TESTING)
    add_executable(
        turbocompgpu_incidence_loss_test
        tests/incidence_loss.cpp
    )

    target_link_libraries(
        turbocompgpu_incidence_loss_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME incidence_loss
        COMMAND turbocompgpu_incidence_loss_test
    )
endif()

# M1.6b.2 Stage 2: Operating-point incidence calibration.
if(BUILD_TESTING)
    add_executable(
        turbocompgpu_operating_incidence_test
        tests/operating_incidence.cpp
    )

    target_link_libraries(
        turbocompgpu_operating_incidence_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME operating_incidence
        COMMAND turbocompgpu_operating_incidence_test
    )
endif()

# M1.6b.2 Stage 3: Incidence-aware nonlinear operating points.
if(BUILD_TESTING)
    add_executable(
        turbocompgpu_incidence_operating_point_test
        tests/incidence_operating_point.cpp
    )

    target_link_libraries(
        turbocompgpu_incidence_operating_point_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME incidence_operating_point
        COMMAND turbocompgpu_incidence_operating_point_test
    )
endif()

# M1.6b.2 Stage 4: Paired incidence compressor-map benchmark.
add_executable(
    turbocompgpu_incidence_map_comparison
    benchmarks/incidence_map_comparison.cpp
)

target_link_libraries(
    turbocompgpu_incidence_map_comparison
    PRIVATE turbocompgpu
)

# M1.6b.2: Analytical zero-incidence locus validation.
if(BUILD_TESTING)
    add_executable(
        turbocompgpu_incidence_locus_test
        tests/incidence_locus.cpp
    )

    target_link_libraries(
        turbocompgpu_incidence_locus_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME incidence_locus
        COMMAND turbocompgpu_incidence_locus_test
    )
endif()

# M1.6b.2: Independent algebraic incidence-loss oracle.
if(BUILD_TESTING)
    add_executable(
        turbocompgpu_incidence_algebraic_test
        tests/incidence_algebraic.cpp
    )

    target_link_libraries(
        turbocompgpu_incidence_algebraic_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME incidence_algebraic
        COMMAND turbocompgpu_incidence_algebraic_test
    )
endif()

# M1.6b.3 Stage 1: Independent diffusion-loss kernel.
if(BUILD_TESTING)
    add_executable(
        turbocompgpu_diffusion_loss_test
        tests/diffusion_loss.cpp
    )

    target_link_libraries(
        turbocompgpu_diffusion_loss_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME diffusion_loss
        COMMAND turbocompgpu_diffusion_loss_test
    )
endif()

# M1.6b.3 Stage 2: Self-consistent diffusion coupling.
if(BUILD_TESTING)
    add_executable(
        turbocompgpu_diffusion_operating_point_test
        tests/diffusion_operating_point.cpp
    )

    target_link_libraries(
        turbocompgpu_diffusion_operating_point_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME diffusion_operating_point
        COMMAND turbocompgpu_diffusion_operating_point_test
    )
endif()

# M1.6b.3 Stage 3: Diffusion-aware compressor-map comparison.
add_executable(
    turbocompgpu_diffusion_map_comparison
    benchmarks/diffusion_map_comparison.cpp
)

target_link_libraries(
    turbocompgpu_diffusion_map_comparison
    PRIVATE turbocompgpu
)

# M1.6b.3 Stage 4A: Diffusion activation-boundary root finder.
if(BUILD_TESTING)
    add_executable(
        turbocompgpu_diffusion_boundary_test
        tests/diffusion_boundary.cpp
    )

    target_link_libraries(
        turbocompgpu_diffusion_boundary_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME diffusion_boundary
        COMMAND turbocompgpu_diffusion_boundary_test
    )
endif()

# M1.6b.3 Stage 4B: Operating diffusion-boundary integration.
if(BUILD_TESTING)
    add_executable(
        turbocompgpu_operating_diffusion_boundary_test
        tests/operating_diffusion_boundary.cpp
    )

    target_link_libraries(
        turbocompgpu_operating_diffusion_boundary_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME operating_diffusion_boundary
        COMMAND turbocompgpu_operating_diffusion_boundary_test
    )
endif()

# M1.6b.3 Stage 4B: Diffusion activation-boundary sweep.
add_executable(
    turbocompgpu_diffusion_boundary_sweep
    benchmarks/diffusion_boundary_sweep.cpp
)

target_link_libraries(
    turbocompgpu_diffusion_boundary_sweep
    PRIVATE turbocompgpu
)

# M1.6b.3: End-to-end diffusion-boundary regression.
if(BUILD_TESTING)
    add_test(
        NAME diffusion_boundary_sweep
        COMMAND turbocompgpu_diffusion_boundary_sweep
            "${CMAKE_CURRENT_BINARY_DIR}/diffusion_boundary_ctest.csv"
    )
endif()

# M1.6b.4 Stage 1: EOS-consistent isentropic sound speed.
if(BUILD_TESTING)
    add_executable(
        turbocompgpu_isentropic_sound_speed_test
        tests/isentropic_sound_speed.cpp
    )

    target_link_libraries(
        turbocompgpu_isentropic_sound_speed_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME isentropic_sound_speed
        COMMAND turbocompgpu_isentropic_sound_speed_test
    )
endif()

# M1.6b.4 Stage 2: EOS-consistent choking mass flux.
if(BUILD_TESTING)
    add_executable(
        turbocompgpu_choking_mass_flux_test
        tests/choking_mass_flux.cpp
    )

    target_link_libraries(
        turbocompgpu_choking_mass_flux_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME choking_mass_flux
        COMMAND turbocompgpu_choking_mass_flux_test
    )
endif()

# M1.6b.4 Stage 3: Independent thermodynamic validation.
if(BUILD_TESTING)
    add_executable(
        turbocompgpu_choking_independent_validation_test
        tests/choking_independent_validation.cpp
    )

    target_link_libraries(
        turbocompgpu_choking_independent_validation_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME choking_independent_validation
        COMMAND turbocompgpu_choking_independent_validation_test
    )
endif()

# M1.6b.4 Stage 4: Rotating-frame throat.
if(BUILD_TESTING)
    add_executable(
        turbocompgpu_rotating_throat_test
        tests/rotating_throat.cpp
    )

    target_link_libraries(
        turbocompgpu_rotating_throat_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME rotating_throat
        COMMAND turbocompgpu_rotating_throat_test
    )
endif()

# M1.6b.4 Stage 5: Inlet-coupled rotating-throat choking.
if(BUILD_TESTING)
    add_executable(
        turbocompgpu_inlet_coupled_choking_test
        tests/inlet_coupled_choking.cpp
    )

    target_link_libraries(
        turbocompgpu_inlet_coupled_choking_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME inlet_coupled_choking
        COMMAND turbocompgpu_inlet_coupled_choking_test
    )
endif()

# M1.6b.4 Stage 6A: Physical throat geometry.
if(BUILD_TESTING)
    add_executable(
        turbocompgpu_throat_geometry_test
        tests/throat_geometry.cpp
    )

    target_link_libraries(
        turbocompgpu_throat_geometry_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME throat_geometry
        COMMAND turbocompgpu_throat_geometry_test
    )
endif()

# M1.6b.4 Stage 6B: Throat-local velocity coordinates.
if(BUILD_TESTING)
    add_executable(
        turbocompgpu_throat_velocity_test
        tests/throat_velocity.cpp
    )

    target_link_libraries(
        turbocompgpu_throat_velocity_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME throat_velocity
        COMMAND turbocompgpu_throat_velocity_test
    )
endif()

# M1.6b.4 Stage 6C: Geometry-coupled inlet choking.
if(BUILD_TESTING)
    add_executable(
        turbocompgpu_geometry_coupled_choking_test
        tests/geometry_coupled_choking.cpp
    )

    target_link_libraries(
        turbocompgpu_geometry_coupled_choking_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME geometry_coupled_choking
        COMMAND turbocompgpu_geometry_coupled_choking_test
    )
endif()

# M1.6b.4 Stage 6D: Prescribed-incidence throat choking.
if(BUILD_TESTING)
    add_executable(
        turbocompgpu_incidence_coupled_choking_test
        tests/incidence_coupled_choking.cpp
    )

    target_link_libraries(
        turbocompgpu_incidence_coupled_choking_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME incidence_coupled_choking
        COMMAND turbocompgpu_incidence_coupled_choking_test
    )
endif()

# M1.6b.4 Stage 6E: Inlet-to-throat incidence kinematics.
if(BUILD_TESTING)
    add_executable(
        turbocompgpu_throat_incidence_test
        tests/throat_incidence.cpp
    )

    target_link_libraries(
        turbocompgpu_throat_incidence_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME throat_incidence
        COMMAND turbocompgpu_throat_incidence_test
    )
endif()

# M1.6b.4 Stage 6E: Mass-flow-coupled throat incidence.
if(BUILD_TESTING)
    add_executable(
        turbocompgpu_predictive_incidence_choking_test
        tests/predictive_incidence_choking.cpp
    )

    target_link_libraries(
        turbocompgpu_predictive_incidence_choking_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME predictive_incidence_choking
        COMMAND turbocompgpu_predictive_incidence_choking_test
    )
endif()

# M1.6b.4 Stage 6F: Parameterized passage turning.
if(BUILD_TESTING)
    add_executable(
        turbocompgpu_passage_turning_choking_test
        tests/passage_turning_choking.cpp
    )

    target_link_libraries(
        turbocompgpu_passage_turning_choking_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME passage_turning_choking
        COMMAND turbocompgpu_passage_turning_choking_test
    )
endif()

# M1.6b.4 Stage 6G: Implicit passage-turning sensitivity.
if(BUILD_TESTING)
    add_executable(
        turbocompgpu_passage_turning_sensitivity_test
        tests/passage_turning_sensitivity.cpp
    )

    target_link_libraries(
        turbocompgpu_passage_turning_sensitivity_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME passage_turning_sensitivity
        COMMAND turbocompgpu_passage_turning_sensitivity_test
    )
endif()

# M1.6b.4 Stage 6H: Choking-boundary diagnostic continuation.
add_executable(
    turbocompgpu_choking_continuation
    benchmarks/choking_continuation.cpp
)

target_link_libraries(
    turbocompgpu_choking_continuation
    PRIVATE turbocompgpu
)

if(BUILD_TESTING)
    add_executable(
        turbocompgpu_choking_continuation_test
        tests/choking_continuation.cpp
    )

    target_link_libraries(
        turbocompgpu_choking_continuation_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME choking_continuation
        COMMAND turbocompgpu_choking_continuation_test
    )
endif()

# M1.6b.4 Stage 6H.1: Constant-area blade-angle comparison.
add_executable(
    turbocompgpu_choking_area_comparison
    benchmarks/choking_area_comparison.cpp
)

target_link_libraries(
    turbocompgpu_choking_area_comparison
    PRIVATE turbocompgpu
)

if(BUILD_TESTING)
    add_executable(
        turbocompgpu_choking_area_comparison_test
        tests/choking_area_comparison.cpp
    )

    target_link_libraries(
        turbocompgpu_choking_area_comparison_test
        PRIVATE turbocompgpu
    )

    add_test(
        NAME choking_area_comparison
        COMMAND turbocompgpu_choking_area_comparison_test
    )
endif()

# Stage 7B: Low-flow compressor characteristic and choking audit.
# Diagnostic benchmark; not registered with CTest.
add_executable(
    turbocompgpu_low_flow_characteristic
    benchmarks/low_flow_characteristic.cpp
)

target_link_libraries(
    turbocompgpu_low_flow_characteristic
    PRIVATE turbocompgpu
)

# Stage 7D: Independently evaluated provisional choking endpoints.
add_executable(
    turbocompgpu_choking_characteristic_endpoints
    benchmarks/choking_characteristic_endpoints.cpp
)

target_link_libraries(
    turbocompgpu_choking_characteristic_endpoints
    PRIVATE turbocompgpu
)

# Stage 7E: Aerodynamic-loss compressor characteristic comparison.
add_executable(
    turbocompgpu_aerodynamic_characteristic
    benchmarks/aerodynamic_characteristic.cpp
)

target_link_libraries(
    turbocompgpu_aerodynamic_characteristic
    PRIVATE turbocompgpu
)

# Stage 7G: Refined characteristic extrema and capacity endpoints.
add_executable(
    turbocompgpu_refined_characteristic
    benchmarks/refined_characteristic.cpp
)

target_link_libraries(
    turbocompgpu_refined_characteristic
    PRIVATE turbocompgpu
)

# Stage 7H: Numerical convergence of compressor characteristic maxima.
add_executable(
    turbocompgpu_characteristic_convergence
    benchmarks/characteristic_convergence.cpp
)

target_link_libraries(
    turbocompgpu_characteristic_convergence
    PRIVATE turbocompgpu
)

# Stage 7I: Aerodynamic closure sensitivity.
add_executable(
    turbocompgpu_aerodynamic_sensitivity
    benchmarks/aerodynamic_sensitivity.cpp
)

target_link_libraries(
    turbocompgpu_aerodynamic_sensitivity
    PRIVATE turbocompgpu
)

# Stage 7I.2: Aerodynamic boundary transitions and peak sensitivities.
add_executable(
    turbocompgpu_aerodynamic_transitions
    benchmarks/aerodynamic_transitions.cpp
)

target_link_libraries(
    turbocompgpu_aerodynamic_transitions
    PRIVATE turbocompgpu
)

# Stage 7I.3: Transition-root and mixed-derivative convergence.
add_executable(
    turbocompgpu_aerodynamic_transition_convergence
    benchmarks/aerodynamic_transition_convergence.cpp
)

target_link_libraries(
    turbocompgpu_aerodynamic_transition_convergence
    PRIVATE turbocompgpu
)

# Stage 7J: Impeller outlet geometry sensitivity.
add_executable(
    turbocompgpu_geometry_sensitivity
    benchmarks/geometry_sensitivity.cpp
)

target_link_libraries(
    turbocompgpu_geometry_sensitivity
    PRIVATE turbocompgpu
)

# Stage 7J.1: Outlet-width capacity-boundary transition.
add_executable(
    turbocompgpu_geometry_width_transition
    benchmarks/geometry_width_transition.cpp
)

target_link_libraries(
    turbocompgpu_geometry_width_transition
    PRIVATE turbocompgpu
)

# Stage 7J.2: Constant-outlet-area geometry sensitivity.
add_executable(
    turbocompgpu_constant_area_sensitivity
    benchmarks/constant_area_sensitivity.cpp
)

target_link_libraries(
    turbocompgpu_constant_area_sensitivity
    PRIVATE turbocompgpu
)

# Stage 7J.3: Constant-area geometry boundary transition.
add_executable(
    turbocompgpu_constant_area_transition
    benchmarks/constant_area_transition.cpp
)

target_link_libraries(
    turbocompgpu_constant_area_transition
    PRIVATE turbocompgpu
)

# Stage 7J.3: Constant-area elasticity convergence.
add_executable(
    turbocompgpu_constant_area_convergence
    benchmarks/constant_area_convergence.cpp
)

target_link_libraries(
    turbocompgpu_constant_area_convergence
    PRIVATE turbocompgpu
)


# Stage 8B.1: lossless complete-stage reference
if(BUILD_TESTING)
    add_executable(turbocompgpu_complete_stage_test
        tests/complete_stage.cpp
    )
    target_link_libraries(turbocompgpu_complete_stage_test
        PRIVATE turbocompgpu
    )
    add_test(
        NAME complete_stage_reference
        COMMAND turbocompgpu_complete_stage_test
    )
endif()


# Stage 8B.2: fixed-geometry vaneless diffuser
if(BUILD_TESTING)
    add_executable(turbocompgpu_vaneless_diffuser_test
        tests/vaneless_diffuser.cpp
    )
    target_link_libraries(turbocompgpu_vaneless_diffuser_test
        PRIVATE turbocompgpu
    )
    add_test(
        NAME vaneless_diffuser
        COMMAND turbocompgpu_vaneless_diffuser_test
    )
endif()


# Stage 8B.3: impeller/diffuser integration
if(BUILD_TESTING)
    add_executable(turbocompgpu_impeller_diffuser_integration_test
        tests/impeller_diffuser_integration.cpp
    )
    target_link_libraries(turbocompgpu_impeller_diffuser_integration_test
        PRIVATE turbocompgpu
    )
    add_test(
        NAME impeller_diffuser_integration
        COMMAND turbocompgpu_impeller_diffuser_integration_test
    )
endif()

# Long-running validation tests are retained but disabled
# during routine development.
option(
    TURBOCOMP_ENABLE_SLOW_TESTS
    "Enable long-running validation tests"
    OFF
)

if(BUILD_TESTING AND NOT TURBOCOMP_ENABLE_SLOW_TESTS)
    set_tests_properties(
        passage_turning_sensitivity
        PROPERTIES DISABLED TRUE
    )
endif()


# Stage 8B.4a: diffuser characteristic analysis
add_executable(turbocompgpu_diffuser_characteristic
    benchmarks/diffuser_characteristic.cpp
)
target_link_libraries(turbocompgpu_diffuser_characteristic
    PRIVATE turbocompgpu
)
if(BUILD_TESTING)
    add_test(
        NAME diffuser_characteristic
        COMMAND turbocompgpu_diffuser_characteristic
    )
    set_tests_properties(
        diffuser_characteristic
        PROPERTIES
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()


# Stage 8B.4b: full characteristic and sonic diagnostics
if(BUILD_TESTING)
    add_test(
        NAME diffuser_sonic_characteristic
        COMMAND turbocompgpu_diffuser_characteristic
    )
    set_tests_properties(
        diffuser_sonic_characteristic
        PROPERTIES
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()


# Stage 8B.4c: branch-aware vaneless diffuser
add_executable(turbocompgpu_vaneless_diffuser_branches_test
    tests/vaneless_diffuser_branches.cpp
)
target_link_libraries(turbocompgpu_vaneless_diffuser_branches_test
    PRIVATE turbocompgpu
)
if(BUILD_TESTING)
    add_test(
        NAME vaneless_diffuser_branches
        COMMAND turbocompgpu_vaneless_diffuser_branches_test
    )
endif()

# Stage 8B.5b.1: fixed-geometry complete-stage characteristic.
add_executable(
    turbocompgpu_complete_stage_characteristic
    benchmarks/complete_stage_characteristic.cpp
)

target_link_libraries(
    turbocompgpu_complete_stage_characteristic
    PRIVATE turbocompgpu
)

# Stage 8B.5b.2: coupled operating-envelope boundaries.
add_executable(
    turbocompgpu_complete_stage_boundaries
    benchmarks/complete_stage_boundaries.cpp
)

target_link_libraries(
    turbocompgpu_complete_stage_boundaries
    PRIVATE turbocompgpu
)

# Stage 8B.5b.3: extended coupled diffuser capacity envelope.
add_executable(
    turbocompgpu_complete_stage_capacity_extended
    benchmarks/complete_stage_capacity_extended.cpp
)

target_link_libraries(
    turbocompgpu_complete_stage_capacity_extended
    PRIVATE turbocompgpu
)
