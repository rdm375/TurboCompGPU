#pragma once

#include <turbocompgpu/pipeline_compressor.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <functional>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace turbocompgpu {

// A complete, rectilinear calibration grid. Values are stored with the last
// axis varying fastest. Queries on the grid boundary are accepted; any
// extrapolation is rejected. CSV loading requires one row at every grid node.
template<std::size_t Dimensions, std::size_t Outputs = 1>
class RectilinearCalibrationTable {
    static_assert(Dimensions > 0 && Dimensions <= 8);
    static_assert(Outputs > 0);

public:
    using Axes = std::array<std::vector<double>, Dimensions>;
    using Coordinate = std::array<double, Dimensions>;
    using Value = std::array<double, Outputs>;

    RectilinearCalibrationTable(Axes axes, std::vector<Value> values)
        : axes_(std::move(axes)), values_(std::move(values))
    {
        validate();
    }

    const Axes& axes() const noexcept { return axes_; }
    const std::vector<Value>& values() const noexcept { return values_; }
    std::pair<double, double> range(std::size_t dimension) const
    {
        if (dimension >= Dimensions)
            throw std::out_of_range("Calibration-map dimension is invalid");
        return {axes_[dimension].front(), axes_[dimension].back()};
    }

    Value evaluate(const Coordinate& coordinate) const
    {
        std::array<std::size_t, Dimensions> lower{};
        std::array<double, Dimensions> fraction{};
        for (std::size_t d = 0; d < Dimensions; ++d) {
            const auto& axis = axes_[d];
            const double x = coordinate[d];
            if (!std::isfinite(x) || x < axis.front() || x > axis.back())
                throw std::out_of_range(
                    "Calibration-map query is outside its tabulated domain");
            if (axis.size() == 1) {
                lower[d] = 0;
                fraction[d] = 0.0;
                continue;
            }
            if (x == axis.back()) {
                lower[d] = axis.size() - 2;
                fraction[d] = 1.0;
                continue;
            }
            const auto upper = std::upper_bound(axis.begin(), axis.end(), x);
            lower[d] = static_cast<std::size_t>(upper - axis.begin() - 1);
            const double x0 = axis[lower[d]];
            fraction[d] = (x - x0) / (axis[lower[d] + 1] - x0);
        }

        Value result{};
        constexpr std::size_t corners = std::size_t{1} << Dimensions;
        for (std::size_t mask = 0; mask < corners; ++mask) {
            std::size_t flat = 0;
            double weight = 1.0;
            for (std::size_t d = 0; d < Dimensions; ++d) {
                const bool upper = ((mask >> d) & 1U) != 0;
                const auto& axis = axes_[d];
                std::size_t index = lower[d];
                if (axis.size() > 1 && upper) ++index;
                weight *= axis.size() == 1
                    ? (upper ? 0.0 : 1.0)
                    : (upper ? fraction[d] : 1.0 - fraction[d]);
                flat = flat * axis.size() + index;
            }
            if (weight == 0.0) continue;
            for (std::size_t output = 0; output < Outputs; ++output)
                result[output] += weight * values_[flat][output];
        }
        return result;
    }

    static RectilinearCalibrationTable load_csv(
        const std::string& path,
        const std::array<std::string, Dimensions>& axis_columns,
        const std::array<std::string, Outputs>& value_columns)
    {
        std::ifstream input(path);
        if (!input)
            throw std::runtime_error("Cannot open calibration CSV: " + path);
        return load_csv(input, axis_columns, value_columns, path);
    }

    static RectilinearCalibrationTable load_csv(
        std::istream& input,
        const std::array<std::string, Dimensions>& axis_columns,
        const std::array<std::string, Outputs>& value_columns,
        const std::string& source_name = "calibration stream")
    {
        std::string line;
        if (!std::getline(input, line))
            throw std::runtime_error("Calibration CSV has no header: " + source_name);
        auto header = split_csv(line);
        if (!header.empty() && header[0].size() >= 3 &&
            static_cast<unsigned char>(header[0][0]) == 0xef &&
            static_cast<unsigned char>(header[0][1]) == 0xbb &&
            static_cast<unsigned char>(header[0][2]) == 0xbf) {
            header[0].erase(0, 3);
        }
        std::array<std::size_t, Dimensions> axis_indices{};
        std::array<std::size_t, Outputs> value_indices{};
        for (std::size_t d = 0; d < Dimensions; ++d)
            axis_indices[d] = find_column(header, axis_columns[d], source_name);
        for (std::size_t o = 0; o < Outputs; ++o)
            value_indices[o] = find_column(header, value_columns[o], source_name);

        struct Row {
            Coordinate coordinate{};
            Value value{};
        };
        std::vector<Row> rows;
        std::array<std::vector<double>, Dimensions> axes;
        std::size_t line_number = 1;
        while (std::getline(input, line)) {
            ++line_number;
            if (trim(line).empty()) continue;
            const auto fields = split_csv(line);
            Row row;
            for (std::size_t d = 0; d < Dimensions; ++d) {
                if (axis_indices[d] >= fields.size())
                    throw std::runtime_error(csv_error(source_name, line_number,
                                                       "missing axis field"));
                row.coordinate[d] = parse_number(
                    fields[axis_indices[d]], source_name, line_number);
                axes[d].push_back(row.coordinate[d]);
            }
            for (std::size_t o = 0; o < Outputs; ++o) {
                if (value_indices[o] >= fields.size())
                    throw std::runtime_error(csv_error(source_name, line_number,
                                                       "missing output field"));
                row.value[o] = parse_number(
                    fields[value_indices[o]], source_name, line_number);
            }
            rows.push_back(row);
        }
        if (rows.empty())
            throw std::runtime_error("Calibration CSV has no data rows: " + source_name);

        for (auto& axis : axes) {
            std::sort(axis.begin(), axis.end());
            axis.erase(std::unique(axis.begin(), axis.end()), axis.end());
        }
        std::size_t grid_size = 1;
        for (const auto& axis : axes) {
            if (axis.size() > std::numeric_limits<std::size_t>::max() / grid_size)
                throw std::runtime_error("Calibration grid is too large: " + source_name);
            grid_size *= axis.size();
        }
        if (rows.size() != grid_size) {
            std::ostringstream message;
            message << "Calibration CSV is not a complete rectilinear grid: "
                    << source_name << " has " << rows.size()
                    << " rows for " << grid_size << " grid nodes";
            throw std::runtime_error(message.str());
        }
        std::vector<Value> values(grid_size);
        std::vector<bool> filled(grid_size, false);
        for (const auto& row : rows) {
            std::size_t flat = 0;
            for (std::size_t d = 0; d < Dimensions; ++d) {
                const auto it = std::lower_bound(
                    axes[d].begin(), axes[d].end(), row.coordinate[d]);
                const auto index = static_cast<std::size_t>(it - axes[d].begin());
                flat = flat * axes[d].size() + index;
            }
            if (filled[flat])
                throw std::runtime_error(
                    "Calibration CSV has duplicate grid coordinates: " + source_name);
            filled[flat] = true;
            values[flat] = row.value;
        }
        if (std::find(filled.begin(), filled.end(), false) != filled.end())
            throw std::runtime_error(
                "Calibration CSV is missing a grid node: " + source_name);
        return RectilinearCalibrationTable(std::move(axes), std::move(values));
    }

private:
    static std::string trim(std::string value)
    {
        const auto first = value.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) return {};
        const auto last = value.find_last_not_of(" \t\r\n");
        return value.substr(first, last - first + 1);
    }

    static std::vector<std::string> split_csv(const std::string& line)
    {
        std::vector<std::string> fields;
        std::string field;
        bool quoted = false;
        for (std::size_t i = 0; i < line.size(); ++i) {
            const char c = line[i];
            if (c == '"') {
                if (quoted && i + 1 < line.size() && line[i + 1] == '"') {
                    field.push_back('"');
                    ++i;
                } else {
                    quoted = !quoted;
                }
            } else if (c == ',' && !quoted) {
                fields.push_back(trim(std::move(field)));
                field.clear();
            } else {
                field.push_back(c);
            }
        }
        if (quoted) throw std::runtime_error("Malformed quoted CSV field");
        fields.push_back(trim(std::move(field)));
        return fields;
    }

    static std::size_t find_column(
        const std::vector<std::string>& header,
        const std::string& name,
        const std::string& source)
    {
        const auto it = std::find(header.begin(), header.end(), name);
        if (it == header.end())
            throw std::runtime_error(
                "Calibration CSV " + source + " is missing column " + name);
        if (std::find(it + 1, header.end(), name) != header.end())
            throw std::runtime_error(
                "Calibration CSV " + source + " has duplicate column " + name);
        return static_cast<std::size_t>(it - header.begin());
    }

    static double parse_number(
        const std::string& field,
        const std::string& source,
        std::size_t line_number)
    {
        std::size_t used = 0;
        double value = 0.0;
        try {
            value = std::stod(field, &used);
        } catch (...) {
            throw std::runtime_error(csv_error(
                source, line_number, "expected a numeric calibration value"));
        }
        if (used != field.size() || !std::isfinite(value))
            throw std::runtime_error(csv_error(
                source, line_number, "calibration values must be finite numbers"));
        return value;
    }

    static std::string csv_error(
        const std::string& source,
        std::size_t line_number,
        const std::string& problem)
    {
        return "Calibration CSV " + source + " line " +
            std::to_string(line_number) + ": " + problem;
    }

    void validate() const
    {
        std::size_t expected = 1;
        for (const auto& axis : axes_) {
            if (axis.empty())
                throw std::invalid_argument("Calibration-map axes cannot be empty");
            for (std::size_t i = 0; i < axis.size(); ++i) {
                if (!std::isfinite(axis[i]) ||
                    (i > 0 && !(axis[i] > axis[i - 1])))
                    throw std::invalid_argument(
                        "Calibration-map axes must be finite and strictly increasing");
            }
            if (axis.size() > std::numeric_limits<std::size_t>::max() / expected)
                throw std::invalid_argument("Calibration-map grid is too large");
            expected *= axis.size();
        }
        if (values_.size() != expected)
            throw std::invalid_argument(
                "Calibration-map values do not fill the rectilinear grid");
        for (const auto& value : values_)
            for (double component : value)
                if (!std::isfinite(component))
                    throw std::invalid_argument(
                        "Calibration-map outputs must be finite");
    }

    Axes axes_;
    std::vector<Value> values_;
};

struct SurgeLineEstimate {
    double surge_flow_kg_s = 0.0;
    double surge_flow_uncertainty_kg_s = 0.0;
};

using SurgeLineTable = RectilinearCalibrationTable<2, 2>;

class TabulatedSurgeLineMap {
public:
    explicit TabulatedSurgeLineMap(SurgeLineTable table)
        : table_(std::move(table))
    {
        for (const auto& value : table_.values())
            if (!(value[0] > 0.0 && value[1] >= 0.0))
                throw std::invalid_argument(
                    "Surge map nodes need positive flow and nonnegative uncertainty");
    }

    SurgeLineEstimate estimate(double corrected_speed, double vane_angle_deg) const
    {
        const auto values = table_.evaluate({corrected_speed, vane_angle_deg});
        if (!(values[0] > 0.0 && values[1] >= 0.0))
            throw std::runtime_error("Surge map returned invalid flow or uncertainty");
        return {values[0], values[1]};
    }

    double surge_flow_kg_s(double corrected_speed, double vane_angle_deg) const
    {
        return estimate(corrected_speed, vane_angle_deg).surge_flow_kg_s;
    }

    double margin_fraction(
        double corrected_speed,
        double vane_angle_deg,
        double corrected_flow_kg_s) const
    {
        if (!(std::isfinite(corrected_flow_kg_s) && corrected_flow_kg_s > 0.0))
            throw std::invalid_argument("Surge-margin flow must be positive and finite");
        const auto line = estimate(corrected_speed, vane_angle_deg);
        return (corrected_flow_kg_s - line.surge_flow_kg_s) /
            line.surge_flow_kg_s;
    }

    double conservative_margin_fraction(
        double corrected_speed,
        double vane_angle_deg,
        double corrected_flow_kg_s,
        double uncertainty_multiplier = 1.0) const
    {
        if (!(std::isfinite(corrected_flow_kg_s) && corrected_flow_kg_s > 0.0 &&
              std::isfinite(uncertainty_multiplier) && uncertainty_multiplier >= 0.0))
            throw std::invalid_argument("Invalid conservative surge-margin query");
        const auto line = estimate(corrected_speed, vane_angle_deg);
        const double conservative_surge = line.surge_flow_kg_s +
            uncertainty_multiplier * line.surge_flow_uncertainty_kg_s;
        return (corrected_flow_kg_s - conservative_surge) / conservative_surge;
    }

private:
    SurgeLineTable table_;
};

using GuideVaneEffectsTable = RectilinearCalibrationTable<3, 4>;

class TabulatedGuideVaneCalibration {
public:
    explicit TabulatedGuideVaneCalibration(GuideVaneEffectsTable table)
        : table_(std::move(table))
    {
        for (const auto& value : table_.values())
            if (!(value[1] >= 0.0 && value[2] >= 0.0 && value[3] >= 0.0))
                throw std::invalid_argument(
                    "Guide-vane nodes need nonnegative entropy and uncertainty bounds");
    }

    GuideVaneMap bind(GuideVaneOperatingDomain domain) const
    {
        const auto angle_range = table_.range(0);
        const auto speed_range = table_.range(1);
        const auto flow_range = table_.range(2);
        if (domain.minimum_angle_deg < angle_range.first ||
            domain.maximum_angle_deg > angle_range.second ||
            domain.minimum_speed_rpm < speed_range.first ||
            domain.maximum_speed_rpm > speed_range.second ||
            domain.minimum_flow_kg_s < flow_range.first ||
            domain.maximum_flow_kg_s > flow_range.second) {
            throw std::invalid_argument(
                "Guide-vane domain exceeds its tabulated calibration grid");
        }
        const auto table = table_;
        GuideVaneMap result;
        result.domain = domain;
        result.calibrated_effects = [table] (
            double angle_deg, double speed_rpm, double flow_kg_s,
            const thermomodel::State&) {
            const auto values = table.evaluate({angle_deg, speed_rpm, flow_kg_s});
            return GuideVaneEffects{values[0], values[1], values[2], values[3]};
        };
        result.validate();
        return result;
    }

private:
    GuideVaneEffectsTable table_;
};

using RecycleValveFlowTable = RectilinearCalibrationTable<4, 2>;

class TabulatedRecycleValveCalibration {
public:
    explicit TabulatedRecycleValveCalibration(RecycleValveFlowTable table)
        : table_(std::move(table))
    {
        for (const auto& value : table_.values())
            if (!(value[0] >= 0.0 && value[1] >= 0.0))
                throw std::invalid_argument(
                    "Recycle-valve nodes need nonnegative flow and uncertainty bound");
    }

    RecycleValveMap bind(RecycleValveDomain domain) const
    {
        const auto suction_range = table_.range(1);
        const auto discharge_range = table_.range(2);
        const auto opening_range = table_.range(0);
        if (opening_range.first < 0.0 || opening_range.second > 1.0)
            throw std::invalid_argument(
                "Recycle-valve opening calibration must lie in [0,1]");
        if (domain.minimum_suction_pressure_pa < suction_range.first ||
            domain.maximum_suction_pressure_pa > suction_range.second ||
            domain.minimum_discharge_pressure_pa < discharge_range.first ||
            domain.maximum_discharge_pressure_pa > discharge_range.second) {
            throw std::invalid_argument(
                "Recycle-valve domain exceeds its tabulated calibration grid");
        }
        const auto table = table_;
        RecycleValveMap result;
        result.domain = domain;
        result.calibrated_flow_kg_s = [table] (
            double opening, double suction_pressure,
            const thermomodel::State& discharge) {
            return table.evaluate({opening, suction_pressure,
                                   discharge.p, discharge.T})[0];
        };
        result.calibrated_flow_uncertainty_kg_s = [table] (
            double opening, double suction_pressure,
            const thermomodel::State& discharge) {
            return table.evaluate({opening, suction_pressure,
                                   discharge.p, discharge.T})[1];
        };
        result.validate();
        return result;
    }

private:
    RecycleValveFlowTable table_;
};

using RecycleCoolerReturnTable = RectilinearCalibrationTable<4, 2>;

class TabulatedRecycleCoolerCalibration {
public:
    explicit TabulatedRecycleCoolerCalibration(RecycleCoolerReturnTable table)
        : table_(std::move(table))
    {
        for (const auto& value : table_.values())
            if (!(value[0] > 0.0 && value[1] >= 0.0))
                throw std::invalid_argument(
                    "Recycle-cooler nodes need positive temperature and nonnegative uncertainty bound");
    }

    double return_temperature_uncertainty_k(
        const thermomodel::State& discharge,
        double suction_pressure,
        double recycle_flow_kg_s) const
    {
        const auto values = table_.evaluate({
            discharge.p, discharge.T, suction_pressure, recycle_flow_kg_s});
        if (!(values[1] >= 0.0))
            throw std::runtime_error(
                "Recycle cooler map returned invalid temperature uncertainty");
        return values[1];
    }

    auto bind(const thermomodel::Model& eos) const
    {
        const auto table = table_;
        const auto* model = &eos;
        return [table, model] (
            const thermomodel::State& discharge,
            double suction_pressure,
            double recycle_flow_kg_s) {
            const double return_temperature = table.evaluate({
                discharge.p, discharge.T, suction_pressure,
                recycle_flow_kg_s})[0];
            if (!(return_temperature > 0.0))
                throw std::runtime_error(
                    "Recycle cooler map returned a nonpositive temperature");
            return evaluate(*model, return_temperature, suction_pressure);
        };
    }

private:
    RecycleCoolerReturnTable table_;
};

} // namespace turbocompgpu
