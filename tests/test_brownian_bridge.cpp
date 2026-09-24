#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

#include "riskengine/methods/montecarlo/brownian_bridge.hpp"

using namespace riskengine;

TEST_CASE("The Brownian bridge is an orthogonal map from normals to increments", "[bridge]") {
    // Increments must again be i.i.d. N(0, 1): the map A (columns = images of unit vectors)
    // must satisfy A A^T = I.
    for (std::uint32_t n : {1u, 2u, 3u, 7u, 12u, 64u, 252u}) {
        INFO("steps " << n);
        const BrownianBridge bridge(n);
        std::vector<std::vector<double>> columns(n, std::vector<double>(n));
        std::vector<double> e(n, 0.0);
        for (std::uint32_t i = 0; i < n; ++i) {
            e[i] = 1.0;
            bridge.transform(e, columns[i]);
            e[i] = 0.0;
        }
        double worst = 0.0;
        for (std::uint32_t r = 0; r < n; ++r)
            for (std::uint32_t c = 0; c < n; ++c) {
                double dot = 0.0;
                for (std::uint32_t i = 0; i < n; ++i) dot += columns[i][r] * columns[i][c];
                worst = std::max(worst, std::abs(dot - (r == c ? 1.0 : 0.0)));
            }
        CHECK(worst <= 1e-13);
    }
}

TEST_CASE("The first normal alone sets the terminal value", "[bridge]") {
    for (std::uint32_t n : {1u, 5u, 12u, 100u}) {
        const BrownianBridge bridge(n);
        std::vector<double> z(n, 0.0), inc(n);
        z[0] = 1.0;
        bridge.transform(z, inc);
        double terminal = 0.0;
        for (double x : inc) terminal += x;
        CHECK(std::abs(terminal - std::sqrt(double(n))) <= 1e-13 * std::sqrt(double(n)));
        // With z[0] alone, the path is the straight line to W(n): equal increments.
        for (double x : inc) CHECK(std::abs(x - 1.0 / std::sqrt(double(n))) <= 1e-14);
        // Every later normal leaves W(n) unchanged: its image sums to zero.
        for (std::uint32_t i = 1; i < n; ++i) {
            std::vector<double> e(n, 0.0);
            e[i] = 1.0;
            bridge.transform(e, inc);
            double sum = 0.0;
            for (double x : inc) sum += x;
            CHECK(std::abs(sum) <= 1e-13);
        }
    }
}
