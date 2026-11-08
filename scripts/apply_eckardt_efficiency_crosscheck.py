#!/usr/bin/env python3
"""Add the digitized eta_tt comparison to the Eckardt friction benchmark."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "benchmarks/eckardt_diffuser_friction_sweep.cpp"
source = SOURCE.read_text()
if "struct EfficiencyPoint" in source or "total_to_total_efficiency" in source:
    raise SystemExit("Efficiency comparison code already appears to be installed")

edits = []

def replace_once(old: str, new: str, description: str) -> None:
    global source
    count = source.count(old)
    if count != 1:
        raise SystemExit(
            f"Expected one source anchor for {description}; found {count}. "
            "No files were changed."
        )
    source = source.replace(old, new, 1)
    edits.append(description)

replace_once(
    "struct MapPoint {\n    double rpm;\n    double mass_flow;\n    double measured_pressure_ratio;\n};\n\n",
    """struct MapPoint {
    double rpm;
    double mass_flow;
    double measured_pressure_ratio;
};

struct EfficiencyPoint {
    double rpm;
    double mass_flow;
    double measured_efficiency;
};

std::vector<EfficiencyPoint> read_efficiency_map(
    const std::filesystem::path& path)
{
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open efficiency data: " + path.string());
    std::vector<EfficiencyPoint> points;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line.front() == '#')
            continue;
        std::istringstream row(line);
        EfficiencyPoint point{};
        char comma1 = '\\0';
        char comma2 = '\\0';
        if (!(row >> point.rpm >> comma1 >> point.mass_flow >> comma2
                  >> point.measured_efficiency) ||
            comma1 != ',' || comma2 != ',')
            throw std::runtime_error("Invalid efficiency row: " + line);
        points.push_back(point);
    }
    if (points.empty())
        throw std::runtime_error("Efficiency map contains no points");
    return points;
}

""",
    "efficiency input parser",
)
replace_once(
    "    double inlet_reynolds;\n};\n\ndouble air_viscosity_sutherland",
    "    double inlet_reynolds;\n    double total_to_total_efficiency;\n};\n\ndouble air_viscosity_sutherland",
    "efficiency prediction field",
)
replace_once(
    "    const double p03 = evaluate_stagnation_state(air, station3)\n        .stagnation_state.p;",
    """    const double p03 = evaluate_stagnation_state(air, station3)
        .stagnation_state.p;
    const auto total1 = evaluate_stagnation_state(air, station1).stagnation_state;
    const auto total3 = evaluate_stagnation_state(air, station3).stagnation_state;
    const double gamma = air.parameters().cp /
        (air.parameters().cp - air.parameters().R);
    const double T03s = total1.T * std::pow(p03 / total1.p,
        (gamma - 1.0) / gamma);
    const double eta_tt = (T03s - total1.T) / (total3.T - total1.T);

""",
    "total-to-total efficiency calculation",
)
replace_once(
    "            100.0 * (1.0 - p03 / p02), used_darcy_factor, reynolds};",
    "            100.0 * (1.0 - p03 / p02), used_darcy_factor, reynolds,\n            eta_tt};",
    "prediction aggregate return",
)
replace_once(
    "            : \"benchmarks/results-eckardt-impeller-a/diffuser_friction_sweep.csv\";\n",
    """            : "benchmarks/results-eckardt-impeller-a/diffuser_friction_sweep.csv";
        const std::optional<std::filesystem::path> efficiency_input = argc > 3
            ? std::optional<std::filesystem::path>{argv[3]} : std::nullopt;
        const std::optional<std::filesystem::path> efficiency_output = argc > 4
            ? std::optional<std::filesystem::path>{argv[4]} : std::nullopt;
""",
    "optional efficiency CLI arguments",
)
replace_once(
    '        std::cout << "Per-point CSV: " << output << \'\\n\';\n',
    """        std::cout << "Per-point CSV: " << output << '\\n';

        if (efficiency_input) {
            if (!efficiency_output)
                throw std::invalid_argument(
                    "Efficiency CSV output path is required when efficiency input is given");
            const auto efficiency_points = read_efficiency_map(*efficiency_input);
            std::ofstream efficiency_csv(*efficiency_output);
            if (!efficiency_csv)
                throw std::runtime_error("Cannot write: " + efficiency_output->string());
            efficiency_csv << std::setprecision(12)
                << "model,darcy_factor,inlet_reynolds,rpm,mass_flow_kg_s,"
                   "measured_eta_tt,predicted_eta_tt,error_percentage_points,"
                   "predicted_total_pressure_ratio,diffuser_loss_percent\\n";
            std::cout << "\\nEfficiency digitization comparison (percentage points):\\n"
                         "Darcy factor, complete, bias_pp, MAE_pp, max_abs_pp\\n";
            for (const double f : friction_factors) {
                double bias = 0.0;
                double mae = 0.0;
                double max_abs = 0.0;
                for (const auto& point : efficiency_points) {
                    const MapPoint request{point.rpm, point.mass_flow, 1.0};
                    const auto result = predict(request, air, g, f);
                    const double error_pp = 100.0 *
                        (result.total_to_total_efficiency - point.measured_efficiency);
                    efficiency_csv << "fixed," << f << ',' << result.inlet_reynolds
                        << ',' << point.rpm << ',' << point.mass_flow << ','
                        << point.measured_efficiency << ','
                        << result.total_to_total_efficiency << ',' << error_pp << ','
                        << result.stage_pr << ',' << result.diffuser_loss_percent << '\\n';
                    bias += error_pp;
                    mae += std::abs(error_pp);
                    max_abs = std::max(max_abs, std::abs(error_pp));
                }
                const double n = static_cast<double>(efficiency_points.size());
                std::cout << f << ", " << n << ", " << bias / n << ", "
                          << mae / n << ", " << max_abs << '\\n';
            }
            double bias = 0.0;
            double mae = 0.0;
            double max_abs = 0.0;
            for (const auto& point : efficiency_points) {
                const MapPoint request{point.rpm, point.mass_flow, 1.0};
                const auto result = predict(request, air, g, -1.0);
                const double error_pp = 100.0 *
                    (result.total_to_total_efficiency - point.measured_efficiency);
                efficiency_csv << "smooth_wall_inlet_Re," << result.darcy_factor
                    << ',' << result.inlet_reynolds << ',' << point.rpm << ','
                    << point.mass_flow << ',' << point.measured_efficiency << ','
                    << result.total_to_total_efficiency << ',' << error_pp << ','
                    << result.stage_pr << ',' << result.diffuser_loss_percent << '\\n';
                bias += error_pp;
                mae += std::abs(error_pp);
                max_abs = std::max(max_abs, std::abs(error_pp));
            }
            const double n = static_cast<double>(efficiency_points.size());
            std::cout << "smooth_wall_inlet_Re, " << n << ", " << bias / n
                      << ", " << mae / n << ", " << max_abs << '\\n'
                      << "Efficiency comparison CSV: " << *efficiency_output << '\\n';
        }
""",
    "efficiency comparison loop",
)

SOURCE.write_text(source)
print("Updated:", SOURCE)
print("Applied anchors:")
for edit in edits:
    print(" -", edit)
