# Optional Enzyme-generated tangent-linear and adjoint sensitivity diagnostic.
# The differentiated residual/objective kernels live in one translation unit;
# thermodynamic state derivatives enter through ThermoGPU's analytic
# derivatives_Tp contract.
option(
    TURBOCOMP_ENABLE_ENZYME_SENSITIVITY
    "Build the optional Enzyme implicit-sensitivity diagnostic"
    OFF
)
set(
    TURBOCOMP_ENZYME_PLUGIN
    ""
    CACHE FILEPATH
    "Path to the LLVM Enzyme pass plugin matching the selected Clang"
)

add_executable(
    turbocompgpu_eckardt_full_stage_sensitivity_fd
    "${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/eckardt_full_stage_sensitivity_fd.cpp"
)
target_link_libraries(
    turbocompgpu_eckardt_full_stage_sensitivity_fd
    PRIVATE turbocompgpu
)

if(BUILD_TESTING)
    add_test(
        NAME eckardt_full_stage_sensitivity_fd
        COMMAND turbocompgpu_eckardt_full_stage_sensitivity_fd
    )
    set_tests_properties(
        eckardt_full_stage_sensitivity_fd
        PROPERTIES TIMEOUT 1800
    )
endif()

if(TURBOCOMP_ENABLE_ENZYME_SENSITIVITY)
    if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        message(FATAL_ERROR
            "TURBOCOMP_ENABLE_ENZYME_SENSITIVITY requires Clang; configure "
            "with -DCMAKE_CXX_COMPILER=clang++")
    endif()
    if(NOT EXISTS "${TURBOCOMP_ENZYME_PLUGIN}")
        message(FATAL_ERROR
            "Set TURBOCOMP_ENZYME_PLUGIN to the Enzyme LLVM pass plugin "
            "built for this Clang/LLVM installation")
    endif()

    add_executable(
        turbocompgpu_eckardt_implicit_sensitivity_enzyme
        "${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/eckardt_implicit_sensitivity_enzyme.cpp"
    )
    target_link_libraries(
        turbocompgpu_eckardt_implicit_sensitivity_enzyme
        PRIVATE turbocompgpu
    )
    target_compile_options(
        turbocompgpu_eckardt_implicit_sensitivity_enzyme
        PRIVATE "-fplugin=${TURBOCOMP_ENZYME_PLUGIN}"
    )

    if(BUILD_TESTING)
        add_test(
            NAME eckardt_implicit_sensitivity_enzyme
            COMMAND turbocompgpu_eckardt_implicit_sensitivity_enzyme
        )
        set_tests_properties(
            eckardt_implicit_sensitivity_enzyme
            PROPERTIES TIMEOUT 120
        )
    endif()
endif()
