#pragma once

#include <thermomodel/model.hpp>

#include <cstddef>

namespace turbocompgpu {

struct ThermodynamicEvaluationCounts {
    std::size_t evaluate_Tp = 0;
    std::size_t evaluate_Trho = 0;

    std::size_t total() const noexcept
    {
        return evaluate_Tp + evaluate_Trho;
    }
};

class CountingModel final : public thermomodel::Model {
public:
    explicit CountingModel(const thermomodel::Model& underlying)
        : underlying_(underlying)
    {
    }

    thermomodel::Result<thermomodel::State>
    evaluate_Tp(
        double temperature,
        double pressure,
        const thermomodel::Composition& composition = {})
        const override
    {
        ++counts_.evaluate_Tp;

        return underlying_.evaluate_Tp(
            temperature, pressure, composition);
    }

    thermomodel::Result<thermomodel::State>
    evaluate_Trho(
        double temperature,
        double density,
        const thermomodel::Composition& composition = {})
        const override
    {
        ++counts_.evaluate_Trho;

        return underlying_.evaluate_Trho(
            temperature, density, composition);
    }

    thermomodel::Capabilities capabilities()
        const noexcept override
    {
        return underlying_.capabilities();
    }

    thermomodel::Metadata metadata() const override
    {
        return underlying_.metadata();
    }

    ThermodynamicEvaluationCounts counts() const noexcept
    {
        return counts_;
    }

    void reset() const noexcept
    {
        counts_ = {};
    }

private:
    const thermomodel::Model& underlying_;

    mutable ThermodynamicEvaluationCounts counts_;
};

} // namespace turbocompgpu
