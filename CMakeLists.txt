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

option(TURBOCOMPGPU_ENABLE_OPENMP
    "Enable OpenMP execution"
    OFF
)

option(TURBOCOMPGPU_ENABLE_CUDA
    "Enable CUDA execution"
    OFF
)

# ThermoGPU integration will be introduced in M1.
