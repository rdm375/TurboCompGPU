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
