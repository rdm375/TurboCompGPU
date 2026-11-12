// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Richard Myers
#pragma once

#include "thermomodel/detail/pr_mixture_state_cuda.hpp"
#include "thermomodel/peng_robinson.hpp"

#include <cuda_runtime_api.h>

#include <cstddef>
#include <memory>
#include <vector>

namespace thermogpu::thermo {

struct PRMixtureTPInput { double T; double p; };

// Owns fixed mixture parameters and can evaluate states either through a
// host convenience method or into caller-owned device arrays. Device calls
// enqueue on the supplied stream and do not synchronize.
class CUDAPRMixtureWorkspace {
public:
    CUDAPRMixtureWorkspace(const thermomodel::PengRobinsonMixture& model,
                           std::size_t capacity);
    ~CUDAPRMixtureWorkspace();
    CUDAPRMixtureWorkspace(const CUDAPRMixtureWorkspace&) = delete;
    CUDAPRMixtureWorkspace& operator=(const CUDAPRMixtureWorkspace&) = delete;
    CUDAPRMixtureWorkspace(CUDAPRMixtureWorkspace&&) noexcept;
    CUDAPRMixtureWorkspace& operator=(CUDAPRMixtureWorkspace&&) noexcept;

    void evaluate_device(const PRMixtureTPInput* device_inputs,
                         thermomodel::detail::PRMixtureCudaState* device_outputs,
                         std::size_t count,
                         cudaStream_t stream = nullptr) const;
    [[nodiscard]] std::vector<thermomodel::detail::PRMixtureCudaState>
    evaluate(const std::vector<PRMixtureTPInput>& inputs) const;
    [[nodiscard]] thermomodel::detail::PRMixtureCudaParameters
    device_parameters() const noexcept;
    [[nodiscard]] std::size_t capacity() const noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace thermogpu::thermo
