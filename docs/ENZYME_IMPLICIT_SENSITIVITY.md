# Enzyme implicit sensitivity diagnostic

This optional diagnostic evaluates tangent-linear (TLM) and adjoint
sensitivities of the converged impeller residual at three Eckardt A operating
flows: 3.99, 4.56, and 5.11 kg/s at 14,000 rpm. It varies slip, exit blade
angle, exit width, exit blockage, and prescribed entropy generation. Its
reported objectives are ideal-gas impeller total-to-total pressure ratio and
efficiency.

The nonlinear solver is run normally to obtain each converged state. Enzyme
then differentiates only the smooth scalar residual and objective kernels.
The state-property chain rule uses the thermodynamic model's
`derivatives_Tp` values. The implicit state tangent solves
`F_x dx/dtheta = -F_theta`; the adjoint solves
`F_x^T lambda = J_x` and evaluates `J_theta - lambda^T F_theta`.

Each reported TLM and adjoint derivative is also compared with a centered
finite-difference rerun of the nonlinear solve. The executable returns an
error if the TLM and adjoint disagree beyond numerical tolerance or if the
finite-difference check differs by more than 0.2% of the derivative scale.

## Build and run

Use the Clang compiler and Enzyme plugin built for the same LLVM version. The
plugin is optional and is not required by the normal TurboCompGPU build.

```sh
cmake -S . -B build-enzyme-sensitivity \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=clang++ \
  -DCMAKE_PREFIX_PATH=/path/to/ThermoGPU/install \
  -DTURBOCOMP_ENABLE_ENZYME_SENSITIVITY=ON \
  -DTURBOCOMP_ENZYME_PLUGIN=/path/to/LLVMEnzyme.so

cmake --build build-enzyme-sensitivity \
  --target turbocompgpu_eckardt_implicit_sensitivity_enzyme --parallel 8

ctest --test-dir build-enzyme-sensitivity --output-on-failure \
  -R '^eckardt_implicit_sensitivity_enzyme$'
```

The output columns contain flow, objective, parameter, TLM derivative, adjoint
derivative, finite-difference derivative, and the two pairwise differences.
These are local derivatives of this impeller-only ideal-gas diagnostic. They
are not calibrated parameter values and do not include the downstream diffuser
or establish physical causality. The current three-point data remain too
sparse to identify several correlated parameters uniquely.

## Differentiability scope

The code deliberately does not differentiate Newton iterations, line searches,
branch selection, exceptions, or validation checks. It differentiates the
converged residual on its current smooth branch and applies the implicit
function theorem. Thermodynamic state evaluations remain outside the active
Enzyme kernel; their first derivatives come from the `DifferentiableModel`
contract. This keeps the differentiated kernel within one C++ translation
unit and avoids requiring Enzyme to differentiate the separately compiled
thermodynamics implementation.
