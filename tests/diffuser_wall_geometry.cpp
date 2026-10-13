#include <turbocompgpu/diffuser_wall_geometry.hpp>
#include <thermomodel/ideal_gas.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>

using namespace turbocompgpu;

namespace {

void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

bool near(double a, double b, double tol = 1e-11) {
    return std::abs(a-b) <=
        tol*std::max({1.0,std::abs(a),std::abs(b)});
}

template<class F>
void rejects(F&& f) {
    bool caught = false;
    try { f(); }
    catch (const std::invalid_argument&) { caught = true; }
    require(caught, "Expected invalid_argument");
}

void test_affine() {
    std::vector<WallPoint> hub, shroud;

    for (int i=0; i<=10; ++i) {
        const double r = 0.20 + i*0.01;
        hub.push_back({r, 0.10+0.02*(r-0.20)});
        shroud.push_back({r, 0.125-0.03*(r-0.20)});
    }

    DiffuserWallGeometry geometry(hub, shroud);

    for (int i=0; i<=100; ++i) {
        const double r = 0.20+i*0.001;
        const double exact = 0.025-0.05*(r-0.20);

        require(near(geometry.width(r), exact, 1e-13),
                "Affine width mismatch");
        require(near(geometry.width_derivative(r), -0.05, 1e-12),
                "Affine derivative mismatch");
    }

    // Check derivative continuity at interior knots.
    for (int i=1; i<10; ++i) {
        const double r = 0.20+i*0.01;
        const double eps = 1e-9;
        require(near(geometry.width_derivative(r-eps),
                     geometry.width_derivative(r+eps),1e-8),
                "Derivative discontinuity");
    }

    // Verify that the generated callbacks remain valid after
    // the original geometry object is destroyed.
    VariableWidthDiffuserGeometry callbacks;
    {
        DiffuserWallGeometry temporary(hub, shroud);
        callbacks = temporary.solver_geometry();
    }
    require(near(callbacks.width(0.25), 0.0225),
            "Geometry callback lifetime failure");

    // Verify compatibility with the existing frictional solver.
    thermomodel::IdealGas model;
    const auto state = evaluate(model, 340.0, 700000.0);
    const double mass =
        state.rho*2*std::numbers::pi*0.20*0.025*60.0;

    VanelessDiffuserInlet inlet{state,60.0,90.0,mass};

    FrictionalDiffuserOptions options;
    options.radial_steps = 512;
    options.darcy_friction_factor = 0.02;

    const auto result =
        evaluate_variable_width_vaneless_diffuser(
            model, model, callbacks, inlet, options);

    require(std::isfinite(result.outlet_static.p),
            "Nonfinite outlet pressure");
    require(std::abs(result.mass_residual) < 1e-8,
            "Mass conservation failure");
    require(std::abs(result.stagnation_enthalpy_residual)
                < 1e-6,
            "Energy conservation failure");

    std::cout << "Affine manufactured geometry: PASS\n";
}

void test_shape_preservation() {
    const std::vector<WallPoint> points{
        {0.20,0.0}, {0.21,0.001},
        {0.22,0.0011}, {0.25,0.006},
        {0.30,0.0061}
    };

    PchipWall wall(points);

    for (std::size_t i=0; i+1<points.size(); ++i) {
        for (int j=0; j<=100; ++j) {
            const double t = j/100.0;
            const double r =
                points[i].r+t*(points[i+1].r-points[i].r);

            const double x = wall.value(r);
            require(x >= points[i].x-1e-14 &&
                    x <= points[i+1].x+1e-14,
                    "PCHIP overshoot");

            // Check analytic derivative against centered FD,
            // away from knots.
            if (j>0 && j<100) {
                const double eps = 1e-7;
                const double fd =
                    (wall.value(r+eps)-wall.value(r-eps)) /
                    (2*eps);
                require(near(wall.derivative(r),fd,2e-6),
                        "PCHIP derivative mismatch");
            }
        }
    }

    std::cout << "Shape preservation and derivatives: PASS\n";
}

void test_rejections() {
    rejects([] {
        PchipWall wall({{0.2,0.0},{0.2,0.1}});
    });

    rejects([] {
        PchipWall wall({{0.3,0.0},{0.2,0.1}});
    });

    rejects([] {
        DiffuserWallGeometry geometry(
            {{0.2,0.02},{0.3,0.02}},
            {{0.2,0.01},{0.3,0.01}});
    });

    rejects([] {
        DiffuserWallGeometry geometry(
            {{0.2,0.0},{0.21,0.0}},
            {{0.3,0.02},{0.31,0.02}});
    });

    std::cout << "Invalid geometry rejection: PASS\n";
}

} // namespace

int main() {
    try {
        test_affine();
        test_shape_preservation();
        test_rejections();
        std::cout << "Diffuser wall geometry verification: PASS\n";
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
