#include <friction/colebrook.hpp>
#include <friction/colebrook_cuda.cuh>
#include <thermogpu/thermo_mixture_batch_cuda.hpp>
#include "lin2026_natural_gas.hpp"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

__global__ void friction_batch_kernel(const double* reynolds,
                                      const double* relative_roughness,
                                      friction::cuda::ColebrookDeviceResult* output,
                                      std::size_t count)
{
    const std::size_t i = static_cast<std::size_t>(blockIdx.x)*blockDim.x+threadIdx.x;
    if (i < count)
        output[i] = friction::cuda::colebrook_turbulent(reynolds[i],relative_roughness[i]);
}

void check(cudaError_t status, const char* operation)
{
    if (status != cudaSuccess)
        throw std::runtime_error(std::string(operation)+": "+cudaGetErrorString(status));
}

template<class T>
T* allocate(std::size_t count)
{
    T* pointer = nullptr;
    check(cudaMalloc(reinterpret_cast<void**>(&pointer),count*sizeof(T)),"cudaMalloc");
    return pointer;
}

bool close(double a, double b, double atol, double rtol)
{
    return std::isfinite(a) && std::isfinite(b) &&
           std::abs(a-b) <= atol + rtol*std::max(std::abs(a),std::abs(b));
}

} // namespace

int main()
{
    try {
        const auto mixture = turbocompgpu::benchmark_data::lin2026_natural_gas();
        const std::vector<thermogpu::thermo::PRMixtureTPInput> states{
            {293.15,6.5e6}, {310.0,7.1e6}, {340.0,8.0e6}, {380.0,9.0e6}};
        thermogpu::thermo::CUDAPRMixtureWorkspace thermo_gpu(mixture,states.size());
        const auto gpu_states = thermo_gpu.evaluate(states);

        const double reynolds_host[]{2.5e4,8.4e5,3.0e6,1.2e7};
        const double roughness_host[]{0.0,1.0e-5,2.0e-4,1.0e-3};
        constexpr std::size_t count = sizeof(reynolds_host)/sizeof(reynolds_host[0]);
        auto* d_re = allocate<double>(count);
        auto* d_rr = allocate<double>(count);
        auto* d_f = allocate<friction::cuda::ColebrookDeviceResult>(count);
        check(cudaMemcpy(d_re,reynolds_host,sizeof(reynolds_host),cudaMemcpyHostToDevice),"copy Re");
        check(cudaMemcpy(d_rr,roughness_host,sizeof(roughness_host),cudaMemcpyHostToDevice),"copy roughness");
        friction_batch_kernel<<<1,128>>>(d_re,d_rr,d_f,count);
        check(cudaGetLastError(),"launch friction kernel");
        check(cudaDeviceSynchronize(),"synchronize friction kernel");
        friction::cuda::ColebrookDeviceResult gpu_f[count]{};
        check(cudaMemcpy(gpu_f,d_f,sizeof(gpu_f),cudaMemcpyDeviceToHost),"copy friction results");
        cudaFree(d_re); cudaFree(d_rr); cudaFree(d_f);

        std::cout << std::setprecision(12)
                  << "state,T_K,p_Pa,Z,rho_kg_m3,h_J_kg,s_J_kgK,max_abs_cpu_error\n";
        for (std::size_t i=0;i<states.size();++i) {
            const auto cpu = mixture.evaluate_Tp(states[i].T,states[i].p);
            const auto cpu_derivatives = mixture.derivatives_Tp(states[i].T,states[i].p);
            if (!cpu || gpu_states[i].status != thermomodel::detail::PRMixtureCudaStatus::Success)
                throw std::runtime_error("CPU or CUDA mixture state evaluation failed");
            if (!cpu_derivatives)
                throw std::runtime_error("CPU PR-mixture derivative evaluation failed");
            const double errors[]{
                std::abs(cpu.value.rho-gpu_states[i].rho),
                std::abs(cpu.value.h-gpu_states[i].h),
                std::abs(cpu.value.s-gpu_states[i].s),
                std::abs(cpu.value.cp-gpu_states[i].cp),
                std::abs(cpu.value.cv-gpu_states[i].cv),
                std::abs(cpu.value.Z-gpu_states[i].Z),
                std::abs(cpu_derivatives.value.drho_dp_T-gpu_states[i].drho_dp_T),
                std::abs(cpu_derivatives.value.drho_dT_p-gpu_states[i].drho_dT_p),
                std::abs(cpu_derivatives.value.dh_dp_T-gpu_states[i].dh_dp_T),
                std::abs(cpu_derivatives.value.dh_dT_p-gpu_states[i].dh_dT_p),
                std::abs(cpu_derivatives.value.ds_dp_T-gpu_states[i].ds_dp_T),
                std::abs(cpu_derivatives.value.ds_dT_p-gpu_states[i].ds_dT_p)};
            double max_error=0.0;
            for(double e:errors) max_error=std::max(max_error,e);
            if (!close(cpu.value.rho,gpu_states[i].rho,1e-9,2e-11) ||
                !close(cpu.value.h,gpu_states[i].h,1e-7,2e-11) ||
                !close(cpu.value.s,gpu_states[i].s,1e-9,2e-11) ||
                !close(cpu.value.cp,gpu_states[i].cp,1e-7,2e-11) ||
                !close(cpu.value.cv,gpu_states[i].cv,1e-7,2e-11) ||
                !close(cpu.value.Z,gpu_states[i].Z,1e-12,2e-11) ||
                !close(cpu_derivatives.value.drho_dp_T,gpu_states[i].drho_dp_T,1e-14,2e-11) ||
                !close(cpu_derivatives.value.drho_dT_p,gpu_states[i].drho_dT_p,1e-12,2e-11) ||
                !close(cpu_derivatives.value.dh_dp_T,gpu_states[i].dh_dp_T,1e-12,2e-11) ||
                !close(cpu_derivatives.value.dh_dT_p,gpu_states[i].dh_dT_p,1e-7,2e-11) ||
                !close(cpu_derivatives.value.ds_dp_T,gpu_states[i].ds_dp_T,1e-14,2e-11) ||
                !close(cpu_derivatives.value.ds_dT_p,gpu_states[i].ds_dT_p,1e-12,2e-11))
                throw std::runtime_error("ThermoGPU CUDA PR-mixture parity check failed");
            std::cout << i << ',' << states[i].T << ',' << states[i].p << ','
                      << gpu_states[i].Z << ',' << gpu_states[i].rho << ','
                      << gpu_states[i].h << ',' << gpu_states[i].s << ','
                      << max_error << '\n';
        }

        std::cout << "Re,relative_roughness,Darcy_f\n";
        for (std::size_t i=0;i<count;++i) {
            const auto cpu=friction::colebrook_turbulent(reynolds_host[i],roughness_host[i]);
            if (gpu_f[i].status != friction::cuda::ColebrookStatus::Success ||
                !close(cpu.darcy_friction_factor,gpu_f[i].darcy_friction_factor,2e-14,2e-11))
                throw std::runtime_error("CUDA Colebrook parity check failed");
            std::cout << reynolds_host[i] << ',' << roughness_host[i] << ','
                      << gpu_f[i].darcy_friction_factor << '\n';
        }
        std::cout << "CUDA PR-mixture and Colebrook parity checks passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
