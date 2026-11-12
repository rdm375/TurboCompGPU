// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Richard Myers

#include "thermogpu/thermo_mixture_batch_cuda.hpp"

#include <cuda_runtime.h>

#include <stdexcept>
#include <string>
#include <utility>

namespace thermogpu::thermo {
namespace {

void cuda_check(cudaError_t status, const char* operation)
{
    if (status != cudaSuccess)
        throw std::runtime_error(std::string(operation) + ": " +
                                 cudaGetErrorString(status));
}

__global__ void mixture_state_kernel(
    thermomodel::detail::PRMixtureCudaParameters parameters,
    const PRMixtureTPInput* inputs,
    thermomodel::detail::PRMixtureCudaState* outputs,
    std::size_t count)
{
    const std::size_t i = static_cast<std::size_t>(blockIdx.x)*blockDim.x+threadIdx.x;
    if (i >= count) return;
    outputs[i] = thermomodel::detail::pr_mixture_evaluate_Tp_cuda(
        parameters, inputs[i].T, inputs[i].p);
}

template<class T>
T* device_allocate(std::size_t count, const char* label)
{
    if (count == 0) return nullptr;
    T* pointer = nullptr;
    cuda_check(cudaMalloc(reinterpret_cast<void**>(&pointer), count*sizeof(T)), label);
    return pointer;
}

} // namespace

struct CUDAPRMixtureWorkspace::Impl {
    thermomodel::detail::PRMixtureCudaParameters parameters{};
    thermomodel::detail::PRMixtureCudaComponent* components = nullptr;
    double* kij = nullptr;
    std::size_t capacity = 0;

    ~Impl() {
        if (components) cudaFree(components);
        if (kij) cudaFree(kij);
    }
};

CUDAPRMixtureWorkspace::CUDAPRMixtureWorkspace(
    const thermomodel::PengRobinsonMixture& model, std::size_t capacity)
    : impl_(std::make_unique<Impl>())
{
    constexpr std::size_t max_components = 32;
    const auto& source = model.components();
    if (source.empty() || source.size() > max_components)
        throw std::invalid_argument("CUDA PR mixture supports 1 to 32 components");
    if (capacity == 0)
        throw std::invalid_argument("CUDA PR mixture capacity must be positive");

    std::vector<thermomodel::detail::PRMixtureCudaComponent> components;
    components.reserve(source.size());
    for (const auto& c : source) {
        components.push_back({c.mole_fraction, c.molar_mass, c.Tc, c.pc,
                              c.acentric, c.ideal_cp});
    }
    const auto& kij = model.binary_interactions();
    if (kij.size() != source.size()*source.size())
        throw std::invalid_argument("PR mixture kij matrix has an invalid size");

    impl_->components = device_allocate<thermomodel::detail::PRMixtureCudaComponent>(
        components.size(), "allocate CUDA PR mixture components");
    impl_->kij = device_allocate<double>(kij.size(), "allocate CUDA PR mixture kij");
    cuda_check(cudaMemcpy(impl_->components, components.data(),
                          components.size()*sizeof(components[0]), cudaMemcpyHostToDevice),
               "copy CUDA PR mixture components");
    cuda_check(cudaMemcpy(impl_->kij, kij.data(), kij.size()*sizeof(double),
                          cudaMemcpyHostToDevice), "copy CUDA PR mixture kij");

    double ideal_cp = 0.0;
    for (const auto& c : source)
        ideal_cp += c.mole_fraction*c.molar_mass*c.ideal_cp;
    ideal_cp /= model.mixture_molar_mass();

    impl_->parameters.components = impl_->components;
    impl_->parameters.binary_interactions = impl_->kij;
    impl_->parameters.component_count = static_cast<int>(source.size());
    impl_->parameters.molar_mass = model.mixture_molar_mass();
    impl_->parameters.gas_constant = thermomodel::detail::pr_Ru /
                                     model.mixture_molar_mass();
    impl_->parameters.ideal_cp = ideal_cp;
    impl_->parameters.Tref = model.reference_temperature();
    impl_->parameters.pref = model.reference_pressure();
    impl_->parameters.Tmin = model.minimum_temperature();
    impl_->parameters.Tmax = model.maximum_temperature();
    impl_->parameters.pmin = model.minimum_pressure();
    impl_->parameters.pmax = model.maximum_pressure();
    impl_->parameters.root_policy = static_cast<int>(model.root_policy());
    impl_->capacity = capacity;
}

CUDAPRMixtureWorkspace::~CUDAPRMixtureWorkspace() = default;
CUDAPRMixtureWorkspace::CUDAPRMixtureWorkspace(CUDAPRMixtureWorkspace&&) noexcept = default;
CUDAPRMixtureWorkspace& CUDAPRMixtureWorkspace::operator=(CUDAPRMixtureWorkspace&&) noexcept = default;

void CUDAPRMixtureWorkspace::evaluate_device(
    const PRMixtureTPInput* inputs,
    thermomodel::detail::PRMixtureCudaState* outputs,
    std::size_t count, cudaStream_t stream) const
{
    if (!impl_ || count > impl_->capacity || (count && (!inputs || !outputs)))
        throw std::invalid_argument("Invalid CUDA PR mixture device batch");
    if (count == 0) return;
    constexpr int threads = 128;
    const int blocks = static_cast<int>((count+threads-1)/threads);
    mixture_state_kernel<<<blocks,threads,0,stream>>>(impl_->parameters,inputs,outputs,count);
    cuda_check(cudaGetLastError(), "launch CUDA PR mixture state kernel");
}

std::vector<thermomodel::detail::PRMixtureCudaState>
CUDAPRMixtureWorkspace::evaluate(const std::vector<PRMixtureTPInput>& inputs) const
{
    if (inputs.size() > capacity())
        throw std::invalid_argument("CUDA PR mixture batch exceeds workspace capacity");
    if (inputs.empty()) return {};
    auto* d_inputs = device_allocate<PRMixtureTPInput>(inputs.size(), "allocate CUDA PR mixture inputs");
    auto* d_outputs = device_allocate<thermomodel::detail::PRMixtureCudaState>(
        inputs.size(), "allocate CUDA PR mixture outputs");
    try {
        cuda_check(cudaMemcpy(d_inputs, inputs.data(), inputs.size()*sizeof(inputs[0]),
                              cudaMemcpyHostToDevice), "copy CUDA PR mixture inputs");
        evaluate_device(d_inputs,d_outputs,inputs.size());
        cuda_check(cudaDeviceSynchronize(), "execute CUDA PR mixture batch");
        std::vector<thermomodel::detail::PRMixtureCudaState> outputs(inputs.size());
        cuda_check(cudaMemcpy(outputs.data(),d_outputs,outputs.size()*sizeof(outputs[0]),
                              cudaMemcpyDeviceToHost), "copy CUDA PR mixture outputs");
        cudaFree(d_inputs); cudaFree(d_outputs);
        return outputs;
    } catch (...) {
        cudaFree(d_inputs); cudaFree(d_outputs);
        throw;
    }
}

thermomodel::detail::PRMixtureCudaParameters
CUDAPRMixtureWorkspace::device_parameters() const noexcept
{
    return impl_ ? impl_->parameters : thermomodel::detail::PRMixtureCudaParameters{};
}

std::size_t CUDAPRMixtureWorkspace::capacity() const noexcept
{
    return impl_ ? impl_->capacity : 0;
}

} // namespace thermogpu::thermo
