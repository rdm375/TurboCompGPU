#include <turbocompgpu/vaneless_diffuser.hpp>
#include <turbocompgpu/isentropic_sound_speed.hpp>

#include <thermomodel/ideal_gas.hpp>
#include <thermomodel/peng_robinson.hpp>
#include <thermomodel/caloric.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace turbocompgpu;

namespace {

struct CharacteristicPoint {
    double pressure;
    double temperature;
    double density;
    double meridional_velocity;
    double tangential_velocity;
    double mass_capacity;
    double sound_speed;
    double meridional_mach;
    double energy_residual;
    double entropy_residual;
};

void analyze(
    const thermomodel::Model& model,
    const std::string& name,
    const std::string& filename)
{
    const ImpellerSystemGeometry geometry{
        {0.10, 0.20, 0.05},
        {0.40, 0.025, 0.10},
        {0.40, 45.0, 0.15}
    };

    const NonlinearResidualScales scales{
        50000.0, 100.0, 500000.0
    };

    const OperatingPointRequest request{
        {
            300.0, 500000.0, 5.0, 12000.0,
            0.90, 20.0, 20.0
        },
        {
            340.0,
            std::log(700000.0 / scales.pressure)
        }
    };

    NonlinearSolverOptions options;
    options.residual_tolerance = 1.0e-13;

    const auto impeller = evaluate_operating_point(
        model, geometry, request, scales, options);

    if (!(impeller.converged &&
          impeller.independently_validated &&
          impeller.basic_physics_valid)) {
        throw std::runtime_error(
            name + ": impeller validation failed");
    }

    const auto state2 = evaluate(
        model,
        impeller.outlet_temperature,
        impeller.outlet_pressure);

    const double pi = std::acos(-1.0);
    const double mass = request.operating.mass_flow;

    const double r2 = geometry.exit.diameter / 2.0;
    const double r3 = 1.5 * r2;

    const double area2 =
        pi * geometry.exit.diameter *
        geometry.exit.width *
        (1.0 - geometry.exit.blockage);

    const double vm2 =
        mass / (state2.rho * area2);

    const double omega =
        2.0 * pi *
        request.operating.rotational_speed_rpm / 60.0;

    const double beta =
        geometry.impeller.exit_blade_angle_deg *
        pi / 180.0;

    const double whirl2 =
        request.operating.slip_factor *
        (omega * r2 - vm2 / std::tan(beta));

    const double whirl3 = r2 * whirl2 / r3;

    const FlowStation station2{
        state2, {vm2, whirl2}
    };

    const auto stagnation2 =
        evaluate_stagnation_state(model, station2);

    const double h02 = stagnation2.stagnation_enthalpy;
    const double p02 = stagnation2.stagnation_state.p;

    const auto* differentiable =
        dynamic_cast<const thermomodel::DifferentiableModel*>(
            &model);

    if (!differentiable) {
        throw std::runtime_error(
            name + ": EOS does not provide analytical derivatives");
    }

    if (!(p02 > state2.p)) {
        throw std::runtime_error(
            name + ": invalid pressure interval");
    }

    // Reproduce the Stage 8B.3 manufactured geometry.
    const double p_reference =
        state2.p + 0.10 * (p02 - state2.p);

    const auto reference = invert_temperature(
        model, p_reference, state2.s,
        Property::Entropy, state2.T);

    const double reference_vm_squared =
        2.0 * (h02 - reference.state.h) -
        whirl3 * whirl3;

    if (!(reference_vm_squared > 0.0)) {
        throw std::runtime_error(
            name + ": reference velocity inadmissible");
    }

    const double reference_vm =
        std::sqrt(reference_vm_squared);

    const double area3 =
        mass / (reference.state.rho * reference_vm);

    const double width3 =
        area3 / (2.0 * pi * r3);

    // Confirm the manufactured geometry is accepted by
    // the existing predictive diffuser implementation.
    const VanelessDiffuserGeometry diffuser_geometry{
        r2, r3, width3, 0.0
    };

    const auto integrated =
        evaluate_vaneless_diffuser_from_operating_point(
            model, geometry, request, impeller,
            diffuser_geometry);

    constexpr int samples = 4000;

    std::vector<CharacteristicPoint> points;
    points.reserve(samples + 1);

    // Extend the pressure domain below impeller-exit pressure.
    // The logarithmic grid gives useful resolution at low pressure.
    //
    // EOS evaluation failures are NOT silently skipped: the scan
    // terminates at the first invalid state encountered while
    // moving from the valid impeller-exit region toward low pressure.

    const auto limits = model.metadata();

    const double p_min =
        std::max(limits.pmin * (1.0 + 1.0e-8),
                 state2.p * 0.01);

    if (!(p_min < state2.p))
        throw std::runtime_error(
            name + ": no accelerating pressure interval");

    auto point_at = [&](double pressure) -> CharacteristicPoint {
        const auto inversion = invert_temperature(
            model, pressure, state2.s,
            Property::Entropy, state2.T);

        const auto& state3 = inversion.state;

        const double vm_squared =
            2.0 * (h02 - state3.h) -
            whirl3 * whirl3;

        if (vm_squared < 0.0) {
            throw std::domain_error(
                "No real meridional velocity");
        }

        const double vm3 = std::sqrt(vm_squared);

        const auto acoustic =
            evaluate_isentropic_sound_speed(
                *differentiable, state3.T, pressure);

        const double capacity =
            state3.rho * area3 * vm3;

        const double energy_residual =
            state3.h +
            0.5 * (vm3 * vm3 + whirl3 * whirl3) -
            h02;

        return {
            pressure,
            state3.T,
            state3.rho,
            vm3,
            whirl3,
            capacity,
            acoustic.sound_speed,
            vm3 / acoustic.sound_speed,
            energy_residual,
            state3.s - state2.s
        };
    };

    // First scan the pressure-recovery region, stopping at the
    // physical kinetic-energy boundary rather than treating it
    // as an EOS failure.

    std::vector<CharacteristicPoint> recovery;

    for (int i = 0; i <= samples; ++i) {
        const double fraction =
            static_cast<double>(i) / samples;

        const double pressure =
            state2.p + fraction * (p02 - state2.p);

        const auto inversion = invert_temperature(
            model, pressure, state2.s,
            Property::Entropy, state2.T);

        const double vm_squared =
            2.0 * (h02 - inversion.state.h) -
            whirl3 * whirl3;

        if (vm_squared < 0.0)
            break;

        recovery.push_back(point_at(pressure));
    }

    // Scan toward lower pressure starting at p2. An EOS domain
    // failure is recorded explicitly; points beyond it are not
    // assumed to belong to the same thermodynamic branch.

    std::vector<CharacteristicPoint> accelerating;

    bool eos_boundary_encountered = false;
    std::string boundary_message;

    for (int i = 1; i <= samples; ++i) {
        const double fraction =
            static_cast<double>(i) / samples;

        const double pressure =
            state2.p *
            std::exp(fraction * std::log(p_min / state2.p));

        try {
            accelerating.push_back(point_at(pressure));
        }
        catch (const std::exception& error) {
            eos_boundary_encountered = true;
            boundary_message = error.what();
            break;
        }
    }

    std::reverse(
        accelerating.begin(), accelerating.end());

    points.insert(
        points.end(),
        accelerating.begin(), accelerating.end());

    points.insert(
        points.end(),
        recovery.begin(), recovery.end());

    if (points.size() < 3) {
        throw std::runtime_error(
            name + ": insufficient admissible samples");
    }

    const auto peak = std::max_element(
        points.begin(), points.end(),
        [](const auto& a, const auto& b) {
            return a.mass_capacity < b.mass_capacity;
        });

    // Count crossings with an explicit tolerance for exact
    // sampled roots. A run of near-zero samples counts once.
    //
    // This is still a sampled-root diagnostic, not a proof that
    // all mathematical roots have been isolated.

    const double root_tolerance =
        1.0e-9 * std::max(1.0, mass);

    auto root_sign = [&](double value) -> int {
        if (value > root_tolerance) return 1;
        if (value < -root_tolerance) return -1;
        return 0;
    };

    int crossings = 0;
    bool in_zero_run = false;
    int previous_nonzero_sign = 0;

    for (const auto& point : points) {
        const int sign =
            root_sign(point.mass_capacity - mass);

        if (sign == 0) {
            if (!in_zero_run) {
                ++crossings;
                in_zero_run = true;
            }
            continue;
        }

        if (!in_zero_run &&
            previous_nonzero_sign != 0 &&
            sign != previous_nonzero_sign) {
            ++crossings;
        }

        in_zero_run = false;
        previous_nonzero_sign = sign;
    }

    // Independently locate the sonic point by solving
    // Vm(p)/a(p) - 1 = 0 with bisection.

    bool sonic_found = false;
    double sonic_pressure = 0.0;
    double sonic_capacity = 0.0;
    double sonic_mach = 0.0;

    for (std::size_t i = 1; i < points.size(); ++i) {
        const double f0 =
            points[i - 1].meridional_mach - 1.0;

        const double f1 =
            points[i].meridional_mach - 1.0;

        if (f0 * f1 > 0.0)
            continue;

        double lo = points[i - 1].pressure;
        double hi = points[i].pressure;
        double flo = f0;

        for (int iteration = 0; iteration < 100; ++iteration) {
            const double mid = lo + 0.5 * (hi - lo);

            if (mid == lo || mid == hi)
                break;

            const double fm =
                point_at(mid).meridional_mach - 1.0;

            if ((fm < 0.0) == (flo < 0.0)) {
                lo = mid;
                flo = fm;
            } else {
                hi = mid;
            }
        }

        const auto sonic = point_at(
            lo + 0.5 * (hi - lo));

        sonic_pressure = sonic.pressure;
        sonic_capacity = sonic.mass_capacity;
        sonic_mach = sonic.meridional_mach;
        sonic_found = true;
        break;
    }

    if (sonic_found) {
        if (std::abs(sonic_mach - 1.0) > 1.0e-8) {
            throw std::runtime_error(
                name + ": sonic root did not converge");
        }

        if (std::abs(sonic_capacity - peak->mass_capacity) >
            2.0e-4 * std::max(1.0, sonic_capacity)) {
            throw std::runtime_error(
                name + ": sampled peak disagrees with sonic capacity");
        }
    }

    std::filesystem::create_directories(
        std::filesystem::path(filename).parent_path());

    std::ofstream csv(filename);

    if (!csv)
        throw std::runtime_error("Cannot open " + filename);

    csv << std::setprecision(17);

    csv
        << "pressure_pa,temperature_k,density_kg_m3,"
        << "vm_m_s,vtheta_m_s,mass_capacity_kg_s,"
        << "sound_speed_m_s,meridional_mach,"
        << "energy_residual_j_kg,"
        << "entropy_residual_j_kg_k\n";

    double maximum_energy_residual = 0.0;
    double maximum_entropy_residual = 0.0;

    for (const auto& point : points) {
        csv
            << point.pressure << ','
            << point.temperature << ','
            << point.density << ','
            << point.meridional_velocity << ','
            << point.tangential_velocity << ','
            << point.mass_capacity << ','
            << point.sound_speed << ','
            << point.meridional_mach << ','
            << point.energy_residual << ','
            << point.entropy_residual << '\n';

        maximum_energy_residual = std::max(
            maximum_energy_residual,
            std::abs(point.energy_residual));

        maximum_entropy_residual = std::max(
            maximum_entropy_residual,
            std::abs(point.entropy_residual));
    }

    if (!(maximum_energy_residual <
              1.0e-7 * std::max(1.0, std::abs(h02)) &&
          maximum_entropy_residual <
              1.0e-7 * std::max(1.0, std::abs(state2.s)))) {
        throw std::runtime_error(
            name + ": characteristic conservation failure");
    }

    // Verify that the predictive solution is represented
    // by the independently sampled characteristic.
    //
    // Interpolate the mass capacity between adjacent samples.
    const double selected_pressure =
        integrated.outlet_static.p;

    bool bracket_found = false;
    double interpolated_capacity = 0.0;

    for (std::size_t i = 1; i < points.size(); ++i) {
        const auto& left = points[i - 1];
        const auto& right = points[i];

        if (left.pressure <= selected_pressure &&
            selected_pressure <= right.pressure) {

            const double fraction =
                (selected_pressure - left.pressure) /
                (right.pressure - left.pressure);

            interpolated_capacity =
                left.mass_capacity +
                fraction * (
                    right.mass_capacity -
                    left.mass_capacity);

            bracket_found = true;
            break;
        }
    }

    if (!bracket_found ||
        std::abs(interpolated_capacity - mass) >
            1.0e-3 * mass) {
        throw std::runtime_error(
            name + ": predictive solution disagrees "
                   "with characteristic");
    }

    std::cout << std::setprecision(10)
              << name << '\n'
              << "  p2                 = " << state2.p << " Pa\n"
              << "  p02                = " << p02 << " Pa\n"
              << "  diffuser width     = " << width3 << " m\n"
              << "  operating mass     = " << mass << " kg/s\n"
              << "  peak sampled mass  = "
              << peak->mass_capacity << " kg/s\n"
              << "  peak pressure      = "
              << peak->pressure << " Pa\n"
              << "  peak at boundary   = "
              << (peak == points.begin() ||
                  peak == points.end() - 1) << '\n'
              << "  sampled crossings  = "
              << crossings << '\n'
              << "  sonic point found  = "
              << sonic_found << '\n'
              << "  sonic pressure     = "
              << sonic_pressure << " Pa\n"
              << "  sonic capacity     = "
              << sonic_capacity << " kg/s\n"
              << "  sonic Mach         = "
              << sonic_mach << '\n'
              << "  EOS boundary hit   = "
              << eos_boundary_encountered << '\n'
              << "  boundary message   = "
              << boundary_message << '\n'
              << "  solver pressure    = "
              << selected_pressure << " Pa\n"
              << "  max energy defect  = "
              << maximum_energy_residual << " J/kg\n"
              << "  max entropy defect = "
              << maximum_entropy_residual << " J/(kg K)\n"
              << "  CSV                = " << filename << '\n';
}

} // namespace

int main()
{
    try {
        const std::string directory =
            "benchmarks/results-diffuser-characteristic/";

        thermomodel::IdealGas ideal;

        analyze(
            ideal,
            "IdealGas",
            directory + "ideal-gas.csv");

        auto caloric =
            std::make_shared<thermomodel::ConstantCp>(2200.0);

        thermomodel::PRParameters parameters;
        parameters.root_policy =
            thermomodel::PRRootPolicy::Vapor;

        thermomodel::PengRobinson pr(
            caloric, parameters);

        analyze(
            pr,
            "PengRobinson",
            directory + "peng-robinson.csv");

        std::cout
            << "PASS: Stage 8B.4a diffuser characteristic\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
