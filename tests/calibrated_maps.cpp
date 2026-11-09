#include <turbocompgpu/calibrated_maps.hpp>

#include <thermomodel/ideal_gas.hpp>

#include <array>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

void require(bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}

void close(double actual, double expected, double tolerance,
           const std::string& message)
{
    if (!(std::abs(actual - expected) <= tolerance)) {
        throw std::runtime_error(message + ": expected " +
            std::to_string(expected) + ", got " + std::to_string(actual));
    }
}

template<class F>
void throws_out_of_range(F&& action, const std::string& message)
{
    bool threw = false;
    try { action(); }
    catch (const std::out_of_range&) { threw = true; }
    require(threw, message);
}

void test_csv_rectilinear_interpolation()
{
    using Table = turbocompgpu::RectilinearCalibrationTable<2, 1>;
    std::istringstream csv(
        "x,y,z,source\n"
        "0,0,0,\"published, point\"\n"
        "0,2,6,a\n"
        "1,0,2,a\n"
        "1,2,8,a\n");
    const auto table = Table::load_csv(
        csv, {"x", "y"}, {"z"}, "linear fixture");
    close(table.evaluate({0.25, 0.5})[0], 2.0, 1.0e-12,
          "bilinear interpolation did not reproduce an affine map");
    throws_out_of_range([&] { table.evaluate({1.01, 1.0}); },
                        "table extrapolated beyond the calibration grid");
}

void test_csv_requires_complete_grid()
{
    using Table = turbocompgpu::RectilinearCalibrationTable<2, 1>;
    std::istringstream csv(
        "x,y,z\n0,0,0\n0,1,1\n1,0,1\n1,0,1\n");
    bool threw = false;
    try { (void)Table::load_csv(csv, {"x", "y"}, {"z"}, "bad fixture"); }
    catch (const std::runtime_error&) { threw = true; }
    require(threw, "loader accepted a duplicate/missing grid node");
}

void test_surge_margin_lookup()
{
    using Table = turbocompgpu::SurgeLineTable;
    Table::Axes axes{{{1.0, 2.0}, {-10.0, 10.0}}};
    std::vector<Table::Value> values{
        {{2.0, 0.2}}, {{4.0, 0.4}}, {{4.0, 0.4}}, {{6.0, 0.6}}};
    const turbocompgpu::TabulatedSurgeLineMap surge(
        Table(std::move(axes), std::move(values)));
    const auto estimate = surge.estimate(1.5, 0.0);
    close(estimate.surge_flow_kg_s, 4.0, 1.0e-12,
          "surge-line bilinear interpolation failed");
    close(estimate.surge_flow_uncertainty_kg_s, 0.4, 1.0e-12,
          "surge-line uncertainty interpolation failed");
    close(surge.surge_flow_kg_s(1.5, 0.0), 4.0, 1.0e-12,
          "surge-line bilinear interpolation failed");
    close(surge.margin_fraction(1.5, 0.0, 5.0), 0.25, 1.0e-12,
          "surge-margin definition is incorrect");
    require(surge.margin_fraction(1.5, 0.0, 3.0) < 0.0,
            "operating point below the surge line did not have negative margin");
    close(surge.conservative_margin_fraction(1.5, 0.0, 5.0),
          0.6 / 4.4, 1.0e-12,
          "uncertainty-adjusted surge margin is incorrect");
    throws_out_of_range([&] { surge.margin_fraction(2.1, 0.0, 5.0); },
                        "surge map extrapolated beyond calibrated speed");
}

void test_guide_vane_table_adapter()
{
    using Table = turbocompgpu::GuideVaneEffectsTable;
    Table::Axes axes{{{0.0, 10.0}, {10000.0, 20000.0}, {1.0, 3.0}}};
    std::vector<Table::Value> values;
    for (double angle : axes[0])
        for (double speed : axes[1])
            for (double flow : axes[2])
                values.push_back({{angle + speed / 10000.0 + flow,
                                   0.001 * (angle + flow), 0.25, 0.0005}});
    turbocompgpu::TabulatedGuideVaneCalibration calibration(
        Table(std::move(axes), std::move(values)));
    auto map = calibration.bind({0.0, 10.0, 10000.0, 20000.0, 1.0, 3.0});
    thermomodel::IdealGas air;
    const auto inlet = turbocompgpu::evaluate(air, 300.0, 200000.0);
    const auto effects = map.evaluate(5.0, 15000.0, 2.0, inlet);
    close(effects.inlet_prewhirl_m_s, 8.5, 1.0e-12,
          "guide-vane table prewhirl interpolation failed");
    close(effects.entropy_generation_j_kg_k, 0.007, 1.0e-12,
          "guide-vane table entropy interpolation failed");
    close(effects.inlet_prewhirl_uncertainty_m_s, 0.25, 1.0e-12,
          "guide-vane prewhirl uncertainty interpolation failed");
    close(effects.entropy_generation_uncertainty_j_kg_k, 0.0005, 1.0e-12,
          "guide-vane entropy uncertainty interpolation failed");
    throws_out_of_range([&] { map.evaluate(11.0, 15000.0, 2.0, inlet); },
                        "guide-vane map extrapolated beyond vane travel");
}

void test_recycle_valve_table_adapter()
{
    using Table = turbocompgpu::RecycleValveFlowTable;
    Table::Axes axes{{{0.0, 1.0}, {100000.0, 200000.0},
                      {300000.0, 400000.0}, {300.0, 400.0}}};
    std::vector<Table::Value> values;
    for (double opening : axes[0])
        for (double suction : axes[1])
            for (double discharge : axes[2])
                for (double temperature : axes[3])
                    values.push_back({{2.0 * opening + suction / 100000.0 +
                                       discharge / 100000.0 + temperature / 100.0,
                                       0.1}});
    turbocompgpu::TabulatedRecycleValveCalibration calibration(
        Table(std::move(axes), std::move(values)));
    auto map = calibration.bind({100000.0, 200000.0,
                                 300000.0, 400000.0, 10.0});
    thermomodel::IdealGas air;
    const auto discharge = turbocompgpu::evaluate(air, 350.0, 350000.0);
    close(map.flow(0.5, 150000.0, discharge), 9.5, 1.0e-12,
          "recycle-valve table interpolation failed");
    close(map.flow_uncertainty(0.5, 150000.0, discharge), 0.1, 1.0e-12,
          "recycle-valve flow uncertainty interpolation failed");
    throws_out_of_range([&] { map.flow(0.5, 90000.0, discharge); },
                        "recycle-valve map extrapolated beyond suction pressure");
}

void test_recycle_cooler_flow_dependent_return()
{
    using Table = turbocompgpu::RecycleCoolerReturnTable;
    Table::Axes axes{{{300000.0, 400000.0}, {350.0, 400.0},
                      {100000.0, 200000.0}, {0.0, 2.0}}};
    std::vector<Table::Value> values;
    for (double p_discharge : axes[0])
        for (double t_discharge : axes[1])
            for (double p_suction : axes[2])
                for (double flow : axes[3])
                    values.push_back({{t_discharge - 40.0 - 2.0 * flow +
                                       p_discharge / 1.0e6 +
                                       p_suction / 1.0e6, 0.8}});
    turbocompgpu::TabulatedRecycleCoolerCalibration cooler(
        Table(std::move(axes), std::move(values)));
    thermomodel::IdealGas air;
    const auto return_model = cooler.bind(air);
    const auto discharge = turbocompgpu::evaluate(air, 375.0, 350000.0);
    const auto returned = return_model(discharge, 150000.0, 1.0);
    close(returned.T, 333.5, 1.0e-10,
          "cooler return map did not use recycle flow");
    close(returned.p, 150000.0, 1.0e-8,
          "cooler return state is not at suction pressure");
    close(cooler.return_temperature_uncertainty_k(discharge, 150000.0, 1.0),
          0.8, 1.0e-12,
          "cooler return-temperature uncertainty interpolation failed");
    static_assert(std::is_invocable_v<decltype(return_model),
                  const thermomodel::State&, double, double>);

    const auto flow_dependent = [&air](
        const thermomodel::State&, double pressure, double flow) {
        return turbocompgpu::evaluate(air, 310.0 + flow, pressure);
    };
    const auto flow_independent = [&air](
        const thermomodel::State&, double pressure) {
        return turbocompgpu::evaluate(air, 315.0, pressure);
    };
    close(turbocompgpu::evaluate_recycle_return_model(
              flow_dependent, discharge, 150000.0, 2.0).T,
          312.0, 1.0e-12,
          "recycle solver adapter did not pass the current valve flow");
    close(turbocompgpu::evaluate_recycle_return_model(
              flow_independent, discharge, 150000.0, 2.0).T,
          315.0, 1.0e-12,
          "two-argument recycle return callback compatibility failed");
}

} // namespace

int main()
{
    try {
        test_csv_rectilinear_interpolation();
        test_csv_requires_complete_grid();
        test_surge_margin_lookup();
        test_guide_vane_table_adapter();
        test_recycle_valve_table_adapter();
        test_recycle_cooler_flow_dependent_return();
        std::cout << "Calibrated-map tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Calibrated-map test failed: " << error.what() << '\n';
        return 1;
    }
}
