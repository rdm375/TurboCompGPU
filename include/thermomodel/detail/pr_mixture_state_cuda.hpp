// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Richard Myers
#pragma once

#include "thermomodel/detail/cubic_kernel.hpp"
#include "thermomodel/detail/pr_coefficients_kernel.hpp"

#include <cmath>

#if defined(__CUDACC__)
#define THERMOMODEL_PR_MIXTURE_CUDA_HD __host__ __device__
#else
#define THERMOMODEL_PR_MIXTURE_CUDA_HD
#endif

namespace thermomodel::detail {

// Flat, pointer-based data contract for CUDA kernels. The pointed-to arrays
// are device resident when called from a kernel. Mole fractions and kij are
// fixed for a compressor map; kij is dense, symmetric and row-major.
struct PRMixtureCudaComponent {
    double mole_fraction;
    double molar_mass;
    double critical_temperature;
    double critical_pressure;
    double acentric_factor;
    double ideal_cp;
};

enum class PRMixtureCudaStatus : int {
    Success = 0,
    InvalidInput = 1,
    OutsideDomain = 2,
    NoPhysicalRoot = 3,
    MultipleRoots = 4,
    UnstableRoot = 5,
    NonphysicalProperties = 6
};

struct PRMixtureCudaParameters {
    const PRMixtureCudaComponent* components = nullptr;
    const double* binary_interactions = nullptr;
    int component_count = 0;
    double molar_mass = 0.0;
    double gas_constant = 0.0;
    double ideal_cp = 0.0;
    double Tref = 298.15;
    double pref = 101325.0;
    double Tmin = 150.0;
    double Tmax = 1000.0;
    double pmin = 100.0;
    double pmax = 100.0e6;
    int root_policy = 0; // 0 vapor, 1 liquid, 2 reject multiple
};

struct PRMixtureCudaState {
    double T = 0.0;
    double p = 0.0;
    double rho = 0.0;
    double h = 0.0;
    double s = 0.0;
    double u = 0.0;
    double cp = 0.0;
    double cv = 0.0;
    double Z = 0.0;
    double sound_speed = 0.0;
    double drho_dp_T = 0.0;
    double drho_dT_p = 0.0;
    double dh_dp_T = 0.0;
    double dh_dT_p = 0.0;
    double ds_dp_T = 0.0;
    double ds_dT_p = 0.0;
    PRMixtureCudaStatus status = PRMixtureCudaStatus::InvalidInput;
};

THERMOMODEL_PR_MIXTURE_CUDA_HD inline PRMixtureCudaState
pr_mixture_evaluate_Tp_cuda(const PRMixtureCudaParameters& q,
                            double T, double p)
{
    constexpr double sqrt2 = 1.41421356237309504880168872420969808;
    constexpr int max_components = 32;
    PRMixtureCudaState out{};
    if (q.components == nullptr || q.binary_interactions == nullptr ||
        q.component_count <= 0 || q.component_count > max_components ||
        !(q.root_policy == 0 || q.root_policy == 1 || q.root_policy == 2) ||
        !(::isfinite(T) && ::isfinite(p) && T > 0.0 && p > 0.0 &&
          ::isfinite(q.molar_mass) && q.molar_mass > 0.0 &&
          ::isfinite(q.gas_constant) && q.gas_constant > 0.0 &&
          ::isfinite(q.ideal_cp) && q.ideal_cp > q.gas_constant)) {
        out.status = PRMixtureCudaStatus::InvalidInput;
        return out;
    }
    if (!(T >= q.Tmin && T <= q.Tmax && p >= q.pmin && p <= q.pmax)) {
        out.status = PRMixtureCudaStatus::OutsideDomain;
        return out;
    }

    double ai[max_components], aip[max_components], aipp[max_components];
    double bi[max_components];
    double a = 0.0, ap = 0.0, app = 0.0, b = 0.0;
    for (int i = 0; i < q.component_count; ++i) {
        const auto& c = q.components[i];
        if (!(::isfinite(c.mole_fraction) && c.mole_fraction >= 0.0 &&
              ::isfinite(c.molar_mass) && c.molar_mass > 0.0 &&
              ::isfinite(c.critical_temperature) && c.critical_temperature > 0.0 &&
              ::isfinite(c.critical_pressure) && c.critical_pressure > 0.0 &&
              ::isfinite(c.acentric_factor) && ::isfinite(c.ideal_cp) && c.ideal_cp > 0.0)) {
            out.status = PRMixtureCudaStatus::InvalidInput;
            return out;
        }
        const double kappa = 0.37464 + 1.54226*c.acentric_factor -
            0.26992*c.acentric_factor*c.acentric_factor;
        const double root_tc_t = ::sqrt(c.critical_temperature*T);
        const double root_t_tc = ::sqrt(T/c.critical_temperature);
        const double qalpha = 1.0 + kappa*(1.0-root_t_tc);
        const double qalpha_p = -kappa/(2.0*root_tc_t);
        const double qalpha_pp = kappa/(4.0*::sqrt(c.critical_temperature)*::pow(T,1.5));
        const double a0 = pr_omega_a*pr_Ru*pr_Ru*
            c.critical_temperature*c.critical_temperature/c.critical_pressure;
        ai[i] = a0*qalpha*qalpha;
        aip[i] = 2.0*a0*qalpha*qalpha_p;
        aipp[i] = 2.0*a0*(qalpha_p*qalpha_p+qalpha*qalpha_pp);
        bi[i] = pr_omega_b*pr_Ru*c.critical_temperature/c.critical_pressure;
        b += c.mole_fraction*bi[i];
    }
    for (int i = 0; i < q.component_count; ++i) {
        for (int j = 0; j < q.component_count; ++j) {
            const double xij = q.components[i].mole_fraction*q.components[j].mole_fraction;
            const double kij = q.binary_interactions[i*q.component_count+j];
            const double aij = (1.0-kij)*::sqrt(ai[i]*ai[j]);
            const double log_deriv = 0.5*(aip[i]/ai[i]+aip[j]/ai[j]);
            const double aij_p = aij*log_deriv;
            const double aij_pp = aij*(0.5*(aipp[i]/ai[i]+aipp[j]/ai[j])-
                0.25*(aip[i]/ai[i]-aip[j]/ai[j])*(aip[i]/ai[i]-aip[j]/ai[j]));
            a += xij*aij;
            ap += xij*aij_p;
            app += xij*aij_pp;
        }
    }
    if (!(::isfinite(a) && a > 0.0 && ::isfinite(b) && b > 0.0)) {
        out.status = PRMixtureCudaStatus::NonphysicalProperties;
        return out;
    }

    const double A = a*p/(pr_Ru*pr_Ru*T*T);
    const double B = b*p/(pr_Ru*T);
    const auto roots = solve_monic_cubic(B-1.0,
        A-3.0*B*B-2.0*B, -A*B+B*B+B*B*B);
    double valid_roots[3]{};
    int valid_count = 0;
    for (int i = 0; i < roots.count; ++i)
        if (::isfinite(roots.values[i]) && roots.values[i] > B*(1.0+1.0e-12))
            valid_roots[valid_count++] = roots.values[i];
    if (valid_count == 0) {
        out.status = PRMixtureCudaStatus::NoPhysicalRoot;
        return out;
    }
    if (valid_count > 1 && q.root_policy == 2) {
        out.status = PRMixtureCudaStatus::MultipleRoots;
        return out;
    }
    const double Z = q.root_policy == 1 ? valid_roots[0] : valid_roots[valid_count-1];
    const double v = Z*pr_Ru*T/p;
    const double D = v*v+2.0*b*v-b*b;
    if (!(::isfinite(v) && v > b && ::isfinite(D) && D > 0.0)) {
        out.status = PRMixtureCudaStatus::NoPhysicalRoot;
        return out;
    }
    const double pv = -pr_Ru*T/((v-b)*(v-b))+a*(2.0*v+2.0*b)/(D*D);
    const double pt = pr_Ru/(v-b)-ap/D;
    if (!(::isfinite(pv) && pv < 0.0 && ::isfinite(pt))) {
        out.status = PRMixtureCudaStatus::UnstableRoot;
        return out;
    }
    const double logq = ::log((v+(1.0+sqrt2)*b)/(v+(1.0-sqrt2)*b));
    const double ur_molar = (T*ap-a)*logq/(2.0*sqrt2*b);
    const double sres_molar = pr_Ru*::log(p*v/(pr_Ru*T)*(v-b)/v)+
        ap*logq/(2.0*sqrt2*b);
    const double cv = q.ideal_cp-q.gas_constant+
        T*app*logq/(2.0*sqrt2*b)/q.molar_mass;
    const double cp = cv-T*(pt/q.molar_mass)*(pt/q.molar_mass)/(pv/q.molar_mass);
    const double c2 = -v*v*pv*cp/(q.molar_mass*cv);
    if (!(::isfinite(cv) && cv > 0.0 && ::isfinite(cp) && cp > cv &&
          ::isfinite(c2) && c2 > 0.0)) {
        out.status = PRMixtureCudaStatus::NonphysicalProperties;
        return out;
    }
    const double rho = q.molar_mass/v;
    const double h = q.ideal_cp*(T-q.Tref)+ur_molar/q.molar_mass+
        (Z-1.0)*q.gas_constant*T;
    const double s = q.ideal_cp*::log(T/q.Tref)-q.gas_constant*::log(p/q.pref)+
        sres_molar/q.molar_mass;
    out.T=T; out.p=p; out.rho=rho; out.h=h; out.s=s; out.u=h-p/rho;
    out.cp=cp; out.cv=cv; out.Z=Z; out.sound_speed=::sqrt(c2);
    out.drho_dp_T=-q.molar_mass/(v*v*pv);
    out.drho_dT_p=q.molar_mass*pt/(v*v*pv);
    out.dh_dp_T=v/q.molar_mass+T*pt/(q.molar_mass*pv);
    out.dh_dT_p=cp;
    out.ds_dp_T=pt/(q.molar_mass*pv);
    out.ds_dT_p=cp/T;
    out.status=PRMixtureCudaStatus::Success;
    return out;
}

} // namespace thermomodel::detail

#undef THERMOMODEL_PR_MIXTURE_CUDA_HD
