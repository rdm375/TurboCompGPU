# TurboCompGPU

Physics-based, GPU-accelerated centrifugal compressor modeling.

TurboCompGPU is an open-source scientific computing project for
modeling centrifugal compressor performance using real-gas
thermodynamics, meanline aerodynamics, and high-performance
numerical methods.

The project builds upon ThermoGPU:

https://github.com/rdm375/ThermoGPU

## Objectives

- Real-gas centrifugal compressor modeling
- Isentropic and polytropic compression calculations
- Impeller and diffuser meanline aerodynamics
- Aerodynamic loss modeling
- Compressor performance-map generation
- Scalar C++, OpenMP, and CUDA execution
- Numerical verification and independent validation
- Integration with transient pipeline simulations

## Development philosophy

1. Establish a correct scalar numerical reference.
2. Validate against analytical and independent external references.
3. Implement and differentially validate parallel backends.
4. Characterize performance using reproducible benchmarks.

## Technology

- C++20
- CMake
- ThermoGPU 2.0
- OpenMP (planned)
- NVIDIA CUDA (planned)

## Status

Initial project development.

## License

GPL-3.0-or-later.
