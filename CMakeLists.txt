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
