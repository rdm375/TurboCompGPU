#pragma once

#include <thermomodel/peng_robinson.hpp>

#include <vector>

namespace turbocompgpu::benchmark_data {

inline thermomodel::PengRobinsonMixture lin2026_natural_gas()
{
    // Lin et al. (2026), Table 1. Component cp values are screening constants;
    // the paper does not publish ideal-gas cp fits or binary kij values.
    return thermomodel::PengRobinsonMixture({
        {"methane", 0.9254, 0.016043, 190.564, 4.5992e6, 0.01142, 2220.0},
        {"nitrogen", 0.0084, 0.0280134, 126.192, 3.3958e6, 0.0372, 1040.0},
        {"carbon_dioxide", 0.0189, 0.0440095, 304.1282, 7.3773e6, 0.22394, 850.0},
        {"ethane", 0.0396, 0.030069, 305.322, 4.8722e6, 0.0995, 1730.0},
        {"propane", 0.0034, 0.0440956, 369.83, 4.248e6, 0.1523, 1670.0},
        {"n-butane", 0.0009, 0.0581222, 425.12, 3.796e6, 0.2002, 1670.0},
        {"isobutane", 0.0012, 0.0581222, 407.81, 3.629e6, 0.1840, 1670.0},
        {"isopentane", 0.0022, 0.0721488, 460.35, 3.378e6, 0.2274, 1700.0}
    }, {}, 298.15, 101325.0, 200.0, 1000.0, 100.0, 100.0e6);
}

} // namespace turbocompgpu::benchmark_data
