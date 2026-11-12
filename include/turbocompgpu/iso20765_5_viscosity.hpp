// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once

#include <thermomodel/peng_robinson.hpp>

#include <cmath>
#include <stdexcept>
#include <string>

#if defined(__CUDACC__)
#define TURBO_ISO20765_HD __host__ __device__
#else
#define TURBO_ISO20765_HD
#endif

namespace turbocompgpu::transport {

TURBO_ISO20765_HD inline bool iso20765_isfinite(double x)
{
#if defined(__CUDA_ARCH__)
    return ::isfinite(x);
#else
    return std::isfinite(x);
#endif
}

// Component critical molar densities are taken from ISO 20765-2:2015 Table B.1.
// ISO 20765-5:2022 LBC viscosity mixing uses these to form the mixture critical
// volume. Keep this fixed-size POD suitable for passing by value to CUDA kernels.
inline constexpr int max_components = 32;

struct Iso20765Component {
    double x;
    double molar_mass_g_mol;
    double Tc_K;
    double pc_MPa;
    double rhoc_mol_dm3;
    double eta_prefactor_mPa_s;
    double mixture_weight;
};

struct Iso20765Parameters {
    int count = 0;
    double molar_mass_kg_mol = 0.0;
    double Mmix_g_mol = 0.0;
    double Tcmix_K = 0.0;
    double pcmix_MPa = 0.0;
    double Vcmix_dm3_mol = 0.0;
    double weighted_sqrt_M = 0.0;
    double xi_mPa_s = 0.0;
    Iso20765Component components[max_components]{};
};

inline double critical_molar_density_iso20765(const std::string& name)
{
    if (name == "methane" || name == "CH4") return 10.139342719;
    if (name == "nitrogen" || name == "N2") return 11.1839;
    if (name == "carbon_dioxide" || name == "carbon dioxide" || name == "CO2") return 10.624978698;
    if (name == "ethane" || name == "C2H6") return 6.870854540;
    if (name == "propane" || name == "C3H8") return 5.000043088;
    if (name == "n-butane" || name == "n_butane" || name == "butane") return 3.920016792;
    if (name == "isobutane" || name == "i-butane") return 3.860142940;
    if (name == "n-pentane" || name == "n_pentane") return 3.215577588;
    if (name == "isopentane" || name == "i-pentane") return 3.271;
    if (name == "n-hexane" || name == "n_hexane") return 2.705877875;
    if (name == "n-heptane" || name == "n_heptane") return 2.315324434;
    if (name == "n-octane" || name == "n_octane") return 2.056404127;
    if (name == "n-nonane" || name == "n_nonane") return 1.81;
    if (name == "n-decane" || name == "n_decane") return 1.64;
    if (name == "hydrogen" || name == "H2") return 14.94;
    if (name == "oxygen" || name == "O2") return 13.63;
    if (name == "carbon_monoxide" || name == "carbon monoxide" || name == "CO") return 10.85;
    if (name == "water" || name == "H2O") return 17.873716090;
    if (name == "hydrogen_sulfide" || name == "hydrogen sulfide" || name == "H2S") return 10.19;
    if (name == "helium" || name == "He") return 17.399;
    if (name == "argon" || name == "Ar") return 13.407429659;
    throw std::invalid_argument("ISO 20765 viscosity has no critical-density entry for component: " + name);
}

inline Iso20765Parameters make_iso20765_parameters(
    const thermomodel::PengRobinsonMixture& mixture)
{
    const auto& components = mixture.components();
    if (components.empty() || components.size() > max_components)
        throw std::invalid_argument("ISO 20765 viscosity supports 1..32 components");

    Iso20765Parameters q{};
    q.count = static_cast<int>(components.size());
    q.molar_mass_kg_mol = mixture.mixture_molar_mass();
    for (int i = 0; i < q.count; ++i) {
        const auto& c = components[static_cast<std::size_t>(i)];
        const double M = c.molar_mass * 1000.0;
        const double pc_MPa = c.pc / 1.0e6;
        const double rhoc = critical_molar_density_iso20765(c.name);
        if (!(std::isfinite(c.mole_fraction) && c.mole_fraction >= 0.0 &&
              std::isfinite(M) && M > 0.0 && std::isfinite(c.Tc) && c.Tc > 0.0 &&
              std::isfinite(pc_MPa) && pc_MPa > 0.0 && rhoc > 0.0))
            throw std::invalid_argument("Invalid ISO 20765 viscosity component data");
        const double sqrt_M = std::sqrt(M);
        const double eta_prefactor = 1.0e-4 * sqrt_M *
            std::pow(c.Tc, -1.0 / 6.0) *
            std::pow(pc_MPa / 0.101325, 2.0 / 3.0);
        if (!(std::isfinite(eta_prefactor) && eta_prefactor > 0.0))
            throw std::invalid_argument("Invalid ISO 20765 component viscosity factor");
        q.components[i] = {
            c.mole_fraction, M, c.Tc, pc_MPa, rhoc,
            eta_prefactor, c.mole_fraction * sqrt_M
        };
        q.Mmix_g_mol += c.mole_fraction * M;
        q.Tcmix_K += c.mole_fraction * c.Tc;
        q.pcmix_MPa += c.mole_fraction * pc_MPa;
        q.Vcmix_dm3_mol += c.mole_fraction / rhoc;
        q.weighted_sqrt_M += c.mole_fraction * sqrt_M;
    }
    if (!(q.Mmix_g_mol > 0.0 && q.Tcmix_K > 0.0 &&
          q.pcmix_MPa > 0.0 && q.Vcmix_dm3_mol > 0.0 &&
          q.weighted_sqrt_M > 0.0 &&
          q.molar_mass_kg_mol > 0.0))
        throw std::invalid_argument("Invalid ISO 20765 mixture parameters");
    q.xi_mPa_s = 1.0e-4 * std::sqrt(q.Mmix_g_mol) *
        std::pow(q.Tcmix_K, -1.0 / 6.0) *
        std::pow(q.pcmix_MPa / 0.101325, 2.0 / 3.0);
    if (!(std::isfinite(q.xi_mPa_s) && q.xi_mPa_s > 0.0))
        throw std::invalid_argument("Invalid ISO 20765 mixture viscosity factor");
    return q;
}

// ISO 20765-5:2022, §5.1 LBC method. Inputs are mass density in kg/m^3,
// mixture molar mass in kg/mol, T in K, and p in Pa. Returns dynamic viscosity
// in Pa s. This method is intended for gas-phase natural gas in its stated
// metering range; it is not a general liquid or near-critical transport model.
TURBO_ISO20765_HD inline double iso20765_5_dynamic_viscosity_Pa_s(
    const Iso20765Parameters& q, double T, double p, double rho_mass)
{
    if (!(q.count > 0 && q.count <= max_components &&
          iso20765_isfinite(T) && T >= 253.15 && T <= 313.15 &&
          iso20765_isfinite(p) && p >= 1.0e5 && p <= 1.0e7 &&
          iso20765_isfinite(rho_mass) && rho_mass > 0.0 &&
          iso20765_isfinite(q.molar_mass_kg_mol) && q.molar_mass_kg_mol > 0.0 &&
          q.Mmix_g_mol > 0.0 && q.Tcmix_K > 0.0 &&
          q.pcmix_MPa > 0.0 && q.Vcmix_dm3_mol > 0.0 &&
          q.weighted_sqrt_M > 0.0 && iso20765_isfinite(q.xi_mPa_s) &&
          q.xi_mPa_s > 0.0)) return NAN;

    double eta_mix_mPa_s = 0.0;
    for (int i = 0; i < q.count; ++i) {
        const auto& c = q.components[i];
        if (!(iso20765_isfinite(c.x) && c.x >= 0.0 && c.molar_mass_g_mol > 0.0 &&
              c.Tc_K > 0.0 && c.pc_MPa > 0.0 && c.rhoc_mol_dm3 > 0.0 &&
              iso20765_isfinite(c.eta_prefactor_mPa_s) &&
              c.eta_prefactor_mPa_s > 0.0 &&
              iso20765_isfinite(c.mixture_weight) && c.mixture_weight >= 0.0)) return NAN;
        const double Tr = T / c.Tc_K;
        const double alpha = Tr <= 1.5
            ? 3.4 * ::pow(Tr, 0.94)
            : 1.778 * ::pow(4.58 * Tr - 1.67, 0.625);
        eta_mix_mPa_s += c.mixture_weight * c.eta_prefactor_mPa_s * alpha;
    }
    eta_mix_mPa_s /= q.weighted_sqrt_M;
    const double molar_density_mol_dm3 =
        (rho_mass / q.molar_mass_kg_mol) / 1000.0;
    const double rho_reduced = q.Vcmix_dm3_mol * molar_density_mol_dm3;
    const double delta = 1.023 + 0.23364 * rho_reduced +
        0.58533 * rho_reduced * rho_reduced -
        0.40758 * rho_reduced * rho_reduced * rho_reduced +
        0.093324 * rho_reduced * rho_reduced * rho_reduced * rho_reduced;
    const double eta_mPa_s = eta_mix_mPa_s + q.xi_mPa_s *
        (delta * delta * delta * delta - 1.0);
    const double eta_Pa_s = eta_mPa_s * 1.0e-3;
    return (iso20765_isfinite(eta_Pa_s) && eta_Pa_s > 0.0) ? eta_Pa_s : NAN;
}

} // namespace turbocompgpu::transport

#undef TURBO_ISO20765_HD
