#pragma once

#include <turbocompgpu/vaneless_diffuser_variable_width.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace turbocompgpu {

// Coordinates are SI: radius r [m], axial position x [m].
// This class does not infer units or coordinate conventions.
struct WallPoint {
    double r;
    double x;
};

class PchipWall {
public:
    explicit PchipWall(std::vector<WallPoint> points)
        : points_(std::move(points))
    {
        if (points_.size() < 2)
            throw std::invalid_argument("PCHIP needs at least two points");

        for (std::size_t i = 0; i < points_.size(); ++i) {
            if (!std::isfinite(points_[i].r) ||
                !std::isfinite(points_[i].x) ||
                points_[i].r <= 0.0)
                throw std::invalid_argument("Invalid wall coordinate");

            if (i && !(points_[i].r > points_[i-1].r))
                throw std::invalid_argument(
                    "Wall radii must be strictly increasing");
        }

        const std::size_t n = points_.size();
        std::vector<double> h(n-1), delta(n-1);

        for (std::size_t i = 0; i+1 < n; ++i) {
            h[i] = points_[i+1].r - points_[i].r;
            delta[i] = (points_[i+1].x - points_[i].x) / h[i];

            if (!std::isfinite(delta[i]))
                throw std::invalid_argument("Invalid wall secant");
        }

        slopes_.assign(n, 0.0);

        if (n == 2) {
            slopes_[0] = slopes_[1] = delta[0];
            return;
        }

        // Fritsch-Butland weighted harmonic interior slopes.
        for (std::size_t i = 1; i+1 < n; ++i) {
            if (delta[i-1] == 0.0 ||
                delta[i] == 0.0 ||
                std::signbit(delta[i-1]) != std::signbit(delta[i]))
                continue;

            const double w1 = 2.0*h[i] + h[i-1];
            const double w2 = h[i] + 2.0*h[i-1];

            // Normalized harmonic expression avoids large products.
            const double a = w1 / (w1+w2);
            const double b = w2 / (w1+w2);
            const double scale =
                std::max(std::abs(delta[i-1]), std::abs(delta[i]));

            const double d0 = delta[i-1] / scale;
            const double d1 = delta[i] / scale;

            slopes_[i] = scale / (a/d0 + b/d1);
        }

        slopes_[0] = endpoint(
            h[0], h[1], delta[0], delta[1]);

        slopes_[n-1] = endpoint(
            h[n-2], h[n-3], delta[n-2], delta[n-3]);

        for (double slope : slopes_)
            if (!std::isfinite(slope))
                throw std::invalid_argument("Nonfinite PCHIP slope");
    }

    double value(double r) const {
        const auto [i,t] = locate(r);
        const double h = points_[i+1].r - points_[i].r;
        const double y0 = points_[i].x;
        const double y1 = points_[i+1].x;

        return (2*t*t*t-3*t*t+1)*y0
             + (t*t*t-2*t*t+t)*h*slopes_[i]
             + (-2*t*t*t+3*t*t)*y1
             + (t*t*t-t*t)*h*slopes_[i+1];
    }

    double derivative(double r) const {
        const auto [i,t] = locate(r);
        const double h = points_[i+1].r - points_[i].r;
        const double y0 = points_[i].x;
        const double y1 = points_[i+1].x;

        return ((6*t*t-6*t)*y0
              + (3*t*t-4*t+1)*h*slopes_[i]
              + (-6*t*t+6*t)*y1
              + (3*t*t-2*t)*h*slopes_[i+1]) / h;
    }

    const std::vector<WallPoint>& points() const {
        return points_;
    }

private:
    static double endpoint(
        double h0, double h1, double d0, double d1)
    {
        double d = ((2*h0+h1)*d0-h0*d1)/(h0+h1);

        if (d0 == 0.0 ||
            std::signbit(d) != std::signbit(d0))
            return 0.0;

        if (std::signbit(d0) != std::signbit(d1) &&
            std::abs(d) > 3.0*std::abs(d0))
            return 3.0*d0;

        return d;
    }

    std::pair<std::size_t,double> locate(double r) const {
        if (!std::isfinite(r) ||
            r < points_.front().r ||
            r > points_.back().r)
            throw std::out_of_range("Radius outside wall profile");

        auto it = std::upper_bound(
            points_.begin(), points_.end(), r,
            [](double value, const WallPoint& p) {
                return value < p.r;
            });

        const std::size_t i =
            it == points_.begin() ? 0 :
            std::min<std::size_t>(
                static_cast<std::size_t>(it-points_.begin()-1),
                points_.size()-2);

        return {i, (r-points_[i].r) /
                   (points_[i+1].r-points_[i].r)};
    }

    std::vector<WallPoint> points_;
    std::vector<double> slopes_;
};

class DiffuserWallGeometry {
public:
    DiffuserWallGeometry(
        std::vector<WallPoint> hub,
        std::vector<WallPoint> shroud)
        : hub_(std::move(hub)),
          shroud_(std::move(shroud))
    {
        rmin_ = std::max(
            hub_.points().front().r,
            shroud_.points().front().r);

        rmax_ = std::min(
            hub_.points().back().r,
            shroud_.points().back().r);

        if (!(rmax_ > rmin_))
            throw std::invalid_argument(
                "Hub and shroud have no common radial interval");

        // Width is a piecewise cubic. Its extrema occur at
        // endpoints or at roots of its quadratic derivative.
        std::vector<double> knots{rmin_, rmax_};

        for (const auto& p : hub_.points())
            if (p.r > rmin_ && p.r < rmax_)
                knots.push_back(p.r);

        for (const auto& p : shroud_.points())
            if (p.r > rmin_ && p.r < rmax_)
                knots.push_back(p.r);

        std::sort(knots.begin(), knots.end());
        knots.erase(
            std::unique(knots.begin(), knots.end()),
            knots.end());

        for (std::size_t i = 0; i+1 < knots.size(); ++i) {
            const double lo = knots[i];
            const double hi = knots[i+1];
            const double h = hi-lo;

            const double w0 = width(lo);
            const double w1 = width(hi);
            const double d0 = width_derivative(lo);
            const double d1 = width_derivative(hi);

            if (!(w0 > 0.0 && w1 > 0.0))
                throw std::invalid_argument(
                    "Nonpositive passage width at knot");

            // w(t)=a*t^3+b*t^2+c*t+d, t in [0,1].
            const double a = 2*w0-2*w1+h*(d0+d1);
            const double b = -3*w0+3*w1-h*(2*d0+d1);
            const double c = h*d0;

            auto check = [&](double t) {
                if (t > 0.0 && t < 1.0) {
                    const double w =
                        ((a*t+b)*t+c)*t+w0;
                    if (!(std::isfinite(w) && w > 0.0))
                        throw std::invalid_argument(
                            "Nonpositive interior passage width");
                }
            };

            // Solve 3*a*t^2 + 2*b*t + c = 0.
            if (a == 0.0) {
                if (b != 0.0)
                    check(-c/(2*b));
            } else {
                const double discriminant = b*b-3*a*c;
                if (discriminant >= 0.0) {
                    const double root = std::sqrt(discriminant);
                    check((-b-root)/(3*a));
                    check((-b+root)/(3*a));
                }
            }
        }
    }

    double inlet_radius() const { return rmin_; }
    double outlet_radius() const { return rmax_; }

    double width(double r) const {
        return shroud_.value(r)-hub_.value(r);
    }

    double width_derivative(double r) const {
        return shroud_.derivative(r)-hub_.derivative(r);
    }

    double hub_axial(double r) const {
        return hub_.value(r);
    }

    double shroud_axial(double r) const {
        return shroud_.value(r);
    }

    // Own the geometry through shared_ptr so callbacks cannot dangle.
    VariableWidthDiffuserGeometry solver_geometry() const {
        auto owned = std::make_shared<DiffuserWallGeometry>(*this);

        return {
            rmin_, rmax_,
            [owned](double r) { return owned->width(r); },
            [owned](double r) {
                return owned->width_derivative(r);
            }
        };
    }

private:
    PchipWall hub_;
    PchipWall shroud_;
    double rmin_{};
    double rmax_{};
};

} // namespace turbocompgpu
