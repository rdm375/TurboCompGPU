#include <turbocompgpu/operating_point.hpp>
#include <turbocompgpu/continuity.hpp>

#include <thermomodel/peng_robinson.hpp>
#ifdef TURBOCOMPGPU_WITH_HEOS
#include <thermomodel/derivatives.hpp>
#include <thermomodel/model.hpp>
#include <Python.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace turbocompgpu;

namespace {

// NIST Shomate ideal-gas cp correlation for CO2, 298–1200 K.
// PR supplies the residual real-gas contribution. This is a reproducible
// caloric input, not a Span–Wagner CO2 backend.
class CarbonDioxideIdealCp final : public thermomodel::CaloricModel {
public:
    double minimum_temperature() const noexcept override { return 200.0; }
    double maximum_temperature() const noexcept override { return 1200.0; }

    double cp(double temperature) const noexcept override
    {
        const double t = temperature / 1000.0;
        return (A + B*t + C*t*t + D*t*t*t + E/(t*t)) / molar_mass;
    }

    double enthalpy_delta(double temperature, double reference) const noexcept override
    {
        const double t = temperature / 1000.0;
        const double r = reference / 1000.0;
        const auto h = [](double x) {
            return A*x + B*x*x/2.0 + C*x*x*x/3.0 +
                   D*x*x*x*x/4.0 - E/x;
        };
        return 1000.0 * (h(t) - h(r)) / molar_mass;
    }

    double entropy_thermal_delta(double temperature, double reference) const noexcept override
    {
        const double t = temperature / 1000.0;
        const double r = reference / 1000.0;
        const auto s = [](double x) {
            return A*std::log(x) + B*x + C*x*x/2.0 +
                   D*x*x*x/3.0 - E/(2.0*x*x);
        };
        return (s(t) - s(r)) / molar_mass;
    }

    double minimum_cp(double tmin, double tmax) const noexcept override
    {
        double result = std::min(cp(tmin), cp(tmax));
        for (int i = 1; i < 512; ++i)
            result = std::min(result, cp(tmin + (tmax-tmin)*i/512.0));
        return result;
    }

private:
    static constexpr double molar_mass = 0.0440095;
    static constexpr double A = 24.99735;
    static constexpr double B = 55.18696;
    static constexpr double C = -33.69137;
    static constexpr double D = 7.948387;
    static constexpr double E = -0.136638;
};

struct Point {
    double rpm, T0, p0, mass, pressure_ratio, p3_bar, p05_bar;
    std::optional<double> T05;
};

std::vector<std::string> split_csv(const std::string& line)
{
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) fields.push_back(field);
    if (!line.empty() && line.back() == ',') fields.emplace_back();
    return fields;
}

std::vector<Point> read_points(const std::filesystem::path& path)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open " + path.string());
    std::vector<Point> points;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line.front() == '#') continue;
        const auto f = split_csv(line);
        if (f.size() != 8) throw std::runtime_error("Expected 8 CSV columns: " + line);
        auto number = [&](std::size_t i) { return std::stod(f.at(i)); };
        Point p{number(0), number(1), number(2), number(3), number(4),
                number(5), number(6), std::nullopt};
        if (!f[7].empty()) p.T05 = number(7);
        points.push_back(p);
    }
    if (points.empty()) throw std::runtime_error("No measurement rows in " + path.string());
    return points;
}

using PRModel = thermomodel::PengRobinson;

#ifdef TURBOCOMPGPU_WITH_HEOS
class CoolPropCO2Model final : public thermomodel::Model,
                               public thermomodel::DifferentiableModel {
public:
    CoolPropCO2Model()
    {
        if (!Py_IsInitialized()) Py_Initialize();
#ifdef TURBOCOMPGPU_COOLPROP_SITEARCH
        PyObject* site_path = PyUnicode_FromString(TURBOCOMPGPU_COOLPROP_SITEARCH);
        if (!site_path || PyList_Append(PySys_GetObject("path"), site_path) != 0) {
            Py_XDECREF(site_path);
            throw std::runtime_error("Could not add CoolProp Python site-packages path");
        }
        Py_DECREF(site_path);
#endif
        constexpr const char* setup = R"PY(
import CoolProp.CoolProp as _cp
_turbocompgpu_co2 = _cp.AbstractState("HEOS", "CO2")
def _turbocompgpu_co2_tp(T, p):
    s = _turbocompgpu_co2
    try:
        s.update(_cp.PT_INPUTS, p, T)
    except ValueError:
        # CoolProp deliberately rejects PT states numerically indistinguishable
        # from saturation. During the inlet entropy bracket the root finder can
        # probe that boundary even when its converged state lies on one side.
        # Retry just off the boundary on the side indicated by the requested p.
        p_sat = _cp.PropsSI("P", "T", T, "Q", 1, "HEOS::CO2")
        if abs(p - p_sat) / p_sat > 2e-6:
            raise
        p_side = p_sat * (1 + 2e-6 if p > p_sat else 1 - 2e-6)
        s.update(_cp.PT_INPUTS, p_side, T)
    return (s.rhomass(), s.hmass(), s.smass(), s.umass(), s.cpmass(),
            s.cvmass(), s.speed_sound(),
            s.first_partial_deriv(_cp.iDmass, _cp.iP, _cp.iT),
            s.first_partial_deriv(_cp.iDmass, _cp.iT, _cp.iP),
            s.first_partial_deriv(_cp.iHmass, _cp.iP, _cp.iT),
            s.first_partial_deriv(_cp.iHmass, _cp.iT, _cp.iP),
            s.first_partial_deriv(_cp.iSmass, _cp.iP, _cp.iT),
            s.first_partial_deriv(_cp.iSmass, _cp.iT, _cp.iP),
            (0 if s.phase() == _cp.iphase_gas else
             3 if s.phase() in (_cp.iphase_supercritical,
                                _cp.iphase_supercritical_gas) else
             1 if s.phase() in (_cp.iphase_liquid,
                                _cp.iphase_supercritical_liquid) else
             2 if s.phase() == _cp.iphase_twophase else 4))
)PY";
        if (PyRun_SimpleString(setup) != 0)
            throw std::runtime_error(python_error("Initializing CoolProp HEOS::CO2"));
        PyObject* main_module = PyImport_AddModule("__main__");
        PyObject* globals = main_module ? PyModule_GetDict(main_module) : nullptr;
        tp_function_ = globals ? PyDict_GetItemString(globals, "_turbocompgpu_co2_tp") : nullptr;
        if (!tp_function_)
            throw std::runtime_error(python_error("Finding CoolProp CO2 state function"));
        Py_INCREF(tp_function_);
    }

    ~CoolPropCO2Model() override
    {
        if (tp_function_ && Py_IsInitialized()) Py_DECREF(tp_function_);
    }

    thermomodel::Result<thermomodel::State> evaluate_Tp(
        double T, double p, const thermomodel::Composition& composition = {}) const override
    {
        if (!composition.mole_fractions.empty())
            return thermomodel::Result<thermomodel::State>::failure(
                thermomodel::ErrorCode::Unsupported, "CO2 backend accepts a pure fluid only");
        const auto values = query(T, p);
        if (!values)
            return thermomodel::Result<thermomodel::State>::failure(
                thermomodel::ErrorCode::OutOfRange, "CoolProp HEOS CO2 state evaluation failed");
        thermomodel::State state;
        state.T = T;
        state.p = p;
        state.rho = (*values)[0];
        state.h = (*values)[1];
        state.s = (*values)[2];
        state.u = (*values)[3];
        state.cp = (*values)[4];
        state.cv = (*values)[5];
        state.sound_speed = (*values)[6];
        state.Z = p / (state.rho * specific_gas_constant * T);
        switch (static_cast<int>((*values)[13])) {
        case 0: state.phase = thermomodel::Phase::Vapor; break;
        case 1: state.phase = thermomodel::Phase::Liquid; break;
        case 2: state.phase = thermomodel::Phase::TwoPhase; break;
        case 3: state.phase = thermomodel::Phase::Supercritical; break;
        default: state.phase = thermomodel::Phase::Unknown; break;
        }
        return thermomodel::Result<thermomodel::State>::success(state);
    }

    thermomodel::Result<thermomodel::State> evaluate_Trho(
        double, double, const thermomodel::Composition& = {}) const override
    {
        return thermomodel::Result<thermomodel::State>::failure(
            thermomodel::ErrorCode::Unsupported,
            "CoolProp CO2 benchmark adapter currently supports T,p states only");
    }

    thermomodel::Result<thermomodel::Derivatives> derivatives_Tp(
        double T, double p, const thermomodel::Composition& composition = {}) const override
    {
        if (!composition.mole_fractions.empty())
            return thermomodel::Result<thermomodel::Derivatives>::failure(
                thermomodel::ErrorCode::Unsupported, "CO2 backend accepts a pure fluid only");
        const auto values = query(T, p);
        if (!values)
            return thermomodel::Result<thermomodel::Derivatives>::failure(
                thermomodel::ErrorCode::OutOfRange, "CoolProp HEOS CO2 derivative evaluation failed");
        return thermomodel::Result<thermomodel::Derivatives>::success({
            (*values)[7], (*values)[8], (*values)[9], (*values)[10],
            (*values)[11], (*values)[12]});
    }

    thermomodel::Capabilities capabilities() const noexcept override
    {
        return {true, false, false, false, true};
    }

    thermomodel::Metadata metadata() const override
    {
        return {"CoolProp HEOS CO2", "CoolProp 8", "CoolProp reference convention",
                216.592, 1100.0, 517950.0, 800.0e6};
    }

private:
    static constexpr double specific_gas_constant = 8.31446261815324 / 0.0440095;
    PyObject* tp_function_ = nullptr;

    static std::string python_error(const char* context)
    {
        if (PyErr_Occurred()) PyErr_Print();
        return std::string(context) + ": Python CoolProp call failed";
    }

    std::optional<std::array<double, 14>> query(double T, double p) const
    {
        const PyGILState_STATE gil = PyGILState_Ensure();
        PyObject* tuple = PyObject_CallFunction(tp_function_, "dd", T, p);
        if (!tuple) {
            PyErr_Print();
            PyGILState_Release(gil);
            return std::nullopt;
        }
        if (!PyTuple_Check(tuple) || PyTuple_Size(tuple) != 14) {
            Py_DECREF(tuple);
            PyGILState_Release(gil);
            return std::nullopt;
        }
        std::array<double, 14> values{};
        for (Py_ssize_t i = 0; i < static_cast<Py_ssize_t>(values.size()); ++i) {
            values[static_cast<std::size_t>(i)] =
                PyFloat_AsDouble(PyTuple_GetItem(tuple, i));
            if (PyErr_Occurred()) {
                PyErr_Print();
                Py_DECREF(tuple);
                PyGILState_Release(gil);
                return std::nullopt;
            }
        }
        Py_DECREF(tuple);
        PyGILState_Release(gil);
        return values;
    }
};

std::unique_ptr<thermomodel::Model> make_co2()
{
    return std::make_unique<CoolPropCO2Model>();
}

const char* co2_backend_name() { return "CoolProp HEOS"; }

#else
std::unique_ptr<thermomodel::Model> make_co2()
{
    thermomodel::PRParameters p;
    p.molar_mass = 0.0440095;
    p.Tc = 304.1282;
    p.pc = 7.3773e6;
    p.acentric = 0.22394;
    p.Tref = 298.15;
    p.pref = 101325.0;
    p.Tmin = 200.0;
    p.Tmax = 1200.0;
    p.pmin = 1.0e5;
    p.pmax = 100.0e6;
    p.root_policy = thermomodel::PRRootPolicy::Vapor;
    return std::make_unique<PRModel>(std::make_shared<CarbonDioxideIdealCp>(), p);
}

const char* co2_backend_name() { return "Peng-Robinson"; }
#endif

struct StaticState { thermomodel::State state; double velocity; };

// Convert measured inlet total state to a static state by enforcing h0,
// s0, continuity, and the documented annular area. This avoids treating
// total measurements as static boundary conditions.
StaticState inlet_static(const thermomodel::Model& eos, const Point& point,
                         double area)
{
#ifdef TURBOCOMPGPU_WITH_HEOS
    constexpr int entropy_iterations = 40;
    constexpr int temperature_iterations = 40;
    constexpr int velocity_iterations = 35;
#else
    constexpr int entropy_iterations = 60;
    constexpr int temperature_iterations = 55;
    constexpr int velocity_iterations = 45;
#endif
    const auto total_result = eos.evaluate_Tp(point.T0, point.p0);
    if (!total_result) throw std::runtime_error("CO2 total-state EOS failed: " + total_result.message);
    const auto total = total_result.value;

    auto state_at_entropy = [&](double temperature) {
        double lo = 1.0e5, hi = point.p0;
        const auto low_state = eos.evaluate_Tp(temperature, lo);
        const auto high_state = eos.evaluate_Tp(temperature, hi);
        if (!low_state || !high_state || low_state.value.s < total.s || high_state.value.s > total.s)
            throw std::runtime_error("Could not bracket static pressure at inlet entropy");
        for (int i = 0; i < entropy_iterations; ++i) {
            const double mid = 0.5 * (lo + hi);
            const auto candidate = eos.evaluate_Tp(temperature, mid);
            if (!candidate) throw std::runtime_error("PR failed in inlet entropy inversion");
            if (candidate.value.s > total.s) lo = mid;
            else hi = mid;
        }
        const auto result = eos.evaluate_Tp(temperature, 0.5 * (lo + hi));
        if (!result) throw std::runtime_error("PR failed at static inlet state");
        return result.value;
    };

    auto static_for_velocity = [&](double velocity) {
        const double target_h = total.h - 0.5 * velocity * velocity;
        double hi = point.T0;
        double interval = 1.0e-6;
        double lo = hi - interval;
        auto low_state = state_at_entropy(lo);
        while (low_state.h > target_h && interval < 20.0) {
            interval = std::min(20.0, 2.0 * interval);
            lo = point.T0 - interval;
            low_state = state_at_entropy(lo);
        }
        if (low_state.h > target_h)
            throw std::runtime_error("Could not bracket inlet static temperature (h_low=" +
                std::to_string(low_state.h) + ", target=" + std::to_string(target_h) +
                ", h0=" + std::to_string(total.h) + ", T0=" +
                std::to_string(point.T0) + ", p0=" + std::to_string(point.p0) +
                ", V=" + std::to_string(velocity) + ")");
        for (int i = 0; i < temperature_iterations; ++i) {
            const double mid = 0.5 * (lo + hi);
            const auto candidate = state_at_entropy(mid);
            if (candidate.h < target_h) lo = mid;
            else hi = mid;
        }
        return state_at_entropy(0.5 * (lo + hi));
    };

#ifdef TURBOCOMPGPU_WITH_HEOS
    // A fixed 0--100 m/s bisection evaluates its first midpoint at 50 m/s.
    // Near-critical HEOS states can leave the single-phase inlet branch at
    // that trial speed even when the continuity solution is at a much lower
    // speed. Grow a valid upper bracket from low velocity instead.
    double lo = 0.0;
    double hi = 1.0;
    bool flow_bracketed = false;
    std::string boundary_failure;
    while (hi <= 100.0) {
        try {
            const auto state = static_for_velocity(hi);
            const double flow = state.rho * area * hi;
            if (flow >= point.mass) {
                flow_bracketed = true;
                break;
            }
            lo = hi;
            hi *= 2.0;
        }
        catch (const std::exception& error) {
            boundary_failure = error.what();
            break;
        }
    }
    if (!flow_bracketed) {
        throw std::runtime_error(
            "Could not bracket inlet velocity on a valid HEOS branch"
            " (last valid V=" + std::to_string(lo) + " m/s, trial V=" +
            std::to_string(hi) + " m/s): " + boundary_failure);
    }
#else
    double lo = 0.0, hi = 100.0;
#endif
    for (int i = 0; i < velocity_iterations; ++i) {
        const double velocity = 0.5 * (lo + hi);
        const auto state = static_for_velocity(velocity);
        const double flow = state.rho * area * velocity;
        if (flow < point.mass) lo = velocity;
        else hi = velocity;
    }
    const double velocity = 0.5 * (lo + hi);
    return {static_for_velocity(velocity), velocity};
}

void write_eos_property_samples(
    const thermomodel::Model& eos,
    const std::vector<Point>& points,
    const std::filesystem::path& output)
{
    if (output.empty()) return;
    if (!output.parent_path().empty())
        std::filesystem::create_directories(output.parent_path());
    std::ofstream csv(output);
    if (!csv) throw std::runtime_error("Cannot write " + output.string());
    csv << std::setprecision(12)
        << "label,T_K,p_Pa,rho_kg_m3,cp_J_kgK,h_J_kg,s_J_kgK,status\n";
    auto write_state = [&](const std::string& label, double T, double p) {
        const auto state = eos.evaluate_Tp(T, p);
        csv << label << ',' << T << ',' << p << ',';
        if (!state) {
            csv << ",,,,ERROR: " << state.message << '\n';
            return;
        }
        csv << state.value.rho << ',' << state.value.cp << ','
            << state.value.h << ',' << state.value.s << ",OK\n";
    };

    // Preserve exact Sandia inlet states so the reference audit covers every
    // compressor-map operating point, not only the two original spot checks.
    for (const auto& point : points) {
        std::ostringstream label;
        label << "sandia_rpm_" << std::fixed << std::setprecision(0) << point.rpm;
        write_state(label.str(), point.T0, point.p0);
    }

    // Dense thermodynamic screen around the Sandia near-critical inlet region.
    write_state("datum_anchor", 300.0, 1.0e6);
    for (double T : {302.0, 304.0, 305.0, 306.0, 308.0, 310.0})
        for (double p : {7.4e6, 7.7e6, 8.0e6, 8.3e6, 8.5e6}) {
            std::ostringstream label;
            label << "grid_T" << std::fixed << std::setprecision(1) << T
                  << "_p" << std::setprecision(2) << p / 1.0e6 << "MPa";
            write_state(label.str(), T, p);
        }
}

ImpellerSystemGeometry sandia_geometry()
{
    constexpr double hub_r = 2.537585e-3;
    constexpr double shroud_r = 9.372047e-3;
    constexpr double exit_r = 18.68170e-3;
    constexpr double blade_thickness = 0.762e-3;
    const double mean_inlet_r = std::sqrt(0.5 * (hub_r*hub_r + shroud_r*shroud_r));
    const double inlet_blockage = 6.0 * blade_thickness /
        (2.0 * std::numbers::pi * mean_inlet_r);
    const double exit_blockage = 12.0 * blade_thickness /
        (2.0 * std::numbers::pi * exit_r);
    return {
        {2.0*hub_r, 2.0*shroud_r, inlet_blockage},
        {2.0*exit_r, 1.8e-3, exit_blockage},
        {2.0*exit_r, 50.0, 2.0*mean_inlet_r}
    };
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const std::filesystem::path input = argc > 1 ? argv[1] :
            "data/validation/sandia_sco2_report_table.csv";
        const std::filesystem::path output = argc > 2 ? argv[2] :
            "sandia_sco2_impeller_screen.csv";
        const std::filesystem::path property_output = argc > 3 ? argv[3] :
            std::filesystem::path{};
        const auto points = read_points(input);
        auto eos = make_co2();
        write_eos_property_samples(*eos, points, property_output);
        if (!property_output.empty())
            std::cout << co2_backend_name() << " property samples: "
                      << property_output << '\n';
        const auto geometry = sandia_geometry();
        const std::array<std::array<double, 4>, 2> co2_reference{{
            {{305.0, 656.77, 7312.5, 8.0e6}},
            {{310.0, 327.71, 9586.4, 8.0e6}}
        }};
        std::cout << "CO2 " << co2_backend_name() << " spot checks against Span-Wagner "
                     "reference (screen only; no acceptance gate):\n";
        for (const auto& reference : co2_reference) {
            const auto state = eos->evaluate_Tp(reference[0], reference[3]);
            if (!state) throw std::runtime_error(
                "CO2 " + std::string(co2_backend_name()) +
                " reference spot check failed: " + state.message);
            std::cout << "  " << reference[0] << " K, " << reference[3]/1.0e6
                      << " MPa: rho=" << state.value.rho << " kg/m3 (ref "
                      << reference[1] << ", error "
                      << 100.0*(state.value.rho/reference[1]-1.0) << "%), cp="
                      << state.value.cp << " J/(kg K) (ref " << reference[2]
                      << ", error " << 100.0*(state.value.cp/reference[2]-1.0)
                      << "%)\n";
        }
        const double inlet_area = effective_flow_area(
            inlet_annulus_area(geometry.inlet), geometry.inlet.blockage);
        constexpr double pressure_scale = 10.0e6;
        const NonlinearResidualScales scales{100000.0, 100.0, pressure_scale};
        NonlinearSolverOptions solver_options;
        solver_options.residual_tolerance = 1.0e-10;

        if (!output.parent_path().empty())
            std::filesystem::create_directories(output.parent_path());
        std::ofstream csv(output);
        if (!csv) throw std::runtime_error("Cannot write " + output.string());
        csv << std::setprecision(12)
            << "rpm,mass_flow_kg_s,inlet_total_T_K,inlet_total_p_Pa,"
               "inlet_static_T_K,inlet_static_p_Pa,inlet_velocity_m_s,"
               "measured_P3_static_Pa,predicted_P3_static_Pa,error_percent,"
               "converged,validation_message\n";

        int converged = 0;
        double abs_error_sum = 0.0;
        double max_abs_error = 0.0;
        constexpr double nan = std::numeric_limits<double>::quiet_NaN();
        for (const auto& point : points) {
            const double measured_p3 = point.p3_bar * 1.0e5;
            double inlet_T = nan, inlet_p = nan, inlet_velocity = nan;
            double predicted = nan, error = nan;
            bool point_valid = false;
            std::string message;
            try {
                const auto inlet = inlet_static(*eos, point, inlet_area);
                inlet_T = inlet.state.T;
                inlet_p = inlet.state.p;
                inlet_velocity = inlet.velocity;
                constexpr double slip = 0.84; // provisional correlation for 6 full blades
                OperatingPointRequest request{
                    {inlet.state.T, inlet.state.p, point.mass, point.rpm,
                     slip, 0.0, 0.0},
                    {point.T0 + 15.0, std::log(1.3 * point.p0 / pressure_scale)}
                };
                const auto result = evaluate_operating_point(
                    *eos, geometry, request, scales, solver_options);
                point_valid = result.independently_validated &&
                              result.basic_physics_valid;
                message = result.message;
                if (point_valid) {
                    predicted = result.outlet_pressure;
                    error = 100.0 * (predicted / measured_p3 - 1.0);
                    abs_error_sum += std::abs(error);
                    max_abs_error = std::max(max_abs_error, std::abs(error));
                    ++converged;
                }
            }
            catch (const std::exception& error_exception) {
                message = error_exception.what();
            }
            for (char& c : message) if (c == ',' || c == '\n' || c == '\r') c = ';';
            csv << point.rpm << ',' << point.mass << ',' << point.T0 << ','
                << point.p0 << ',' << inlet_T << ',' << inlet_p << ','
                << inlet_velocity << ',' << measured_p3 << ',' << predicted << ','
                << error << ',' << point_valid
                << ',' << message << '\n';
            std::cout << point.rpm << " rpm " << point.mass << " kg/s: "
                      << (std::isfinite(error) ? "P3 error=" + std::to_string(error) + "%" :
                          "no valid impeller result: " + message) << '\n';
        }
        std::cout << "Sandia CO2 " << co2_backend_name()
                  << " impeller-only screen: " << converged << '/' << points.size()
                  << " valid, MAPE=" << (converged ? abs_error_sum/converged : 0.0)
                  << "%, max_abs=" << max_abs_error << "%\n"
                  << "CSV: " << output << '\n';
#ifdef TURBOCOMPGPU_WITH_HEOS
        // This is a diagnostic screen without an accuracy gate. Preserve its
        // partial per-point results for invalid/saturated inlet reconstructions.
        return converged > 0 ? 0 : 2;
#else
        return converged == static_cast<int>(points.size()) ? 0 : 2;
#endif
    }
    catch (const std::exception& e) {
        std::cerr << "Sandia CO2 screen failed: " << e.what() << '\n';
        return 1;
    }
}
