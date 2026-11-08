#pragma once

#include <turbocompgpu/pipeline_compressor.hpp>

#include <functional>

namespace turbocompgpu {

// A0/A1: the external pipe solver owns the pipe equations and reports the
// resulting station measurements. The evaluator must perform a complete
// steady pipe-compressor-network solve for the supplied driver speed.
using SteadyPlantEvaluator =
    std::function<SteadyPlantMeasurements(double speed_rpm)>;

} // namespace turbocompgpu
