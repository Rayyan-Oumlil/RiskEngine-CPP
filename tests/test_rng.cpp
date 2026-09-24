#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "riskengine/core/normal.hpp"
#include "riskengine/core/rng/random_stream.hpp"

using namespace riskengine;

// Known-answer vectors from Random123 (tests/kat_vectors, philox4x32 10). Checked at compile time.
static_assert(philox::generate({0, 0, 0, 0}, {0, 0}) ==
              philox::Counter{0x6627e8d5, 0xe169c58d, 0xbc57ac4c, 0x9b00dbd8});
static_assert(philox::generate({0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff}, {0xffffffff, 0xffffffff}) ==
              philox::Counter{0x408f276d, 0x41c83b0e, 0xa20bc7c6, 0x6d5451fd});
static_assert(philox::generate({0x243f6a88, 0x85a308d3, 0x13198a2e, 0x03707344}, {0xa4093822, 0x299f31d0}) ==
              philox::Counter{0xd16cfe09, 0x94fdcceb, 0x5001e420, 0x24126ea1});

namespace {

constexpr std::size_t kSamples = 1'000'000;

std::vector<double> draw(SeedKey key, std::uint32_t block, std::size_t n, bool normal) {
    RandomStream rng(key, block);
    std::vector<double> out(n);
    for (double& x : out) x = normal ? rng.normal() : rng.uniform();
    return out;
}

// Kolmogorov-Smirnov statistic sqrt(n) * sup |F_n - F|.
template <class Cdf>
double ks_statistic(std::vector<double> xs, Cdf cdf) {
    std::sort(xs.begin(), xs.end());
    const double n = static_cast<double>(xs.size());
    double d = 0.0;
    for (std::size_t i = 0; i < xs.size(); ++i) {
        const double f = cdf(xs[i]);
        d = std::max({d, f - static_cast<double>(i) / n, static_cast<double>(i + 1) / n - f});
    }
    return std::sqrt(n) * d;
}

double mean(const std::vector<double>& xs) {
    double s = 0.0;
    for (double x : xs) s += x;
    return s / static_cast<double>(xs.size());
}

double central_moment(const std::vector<double>& xs, double m, int k) {
    double s = 0.0;
    for (double x : xs) s += std::pow(x - m, k);
    return s / static_cast<double>(xs.size());
}

double correlation(const std::vector<double>& a, const std::vector<double>& b) {
    const double ma = mean(a), mb = mean(b);
    double sab = 0.0, saa = 0.0, sbb = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        sab += (a[i] - ma) * (b[i] - mb);
        saa += (a[i] - ma) * (a[i] - ma);
        sbb += (b[i] - mb) * (b[i] - mb);
    }
    return sab / std::sqrt(saa * sbb);
}

// 99.9 % critical value of the Kolmogorov distribution. The seeds are fixed, so these tests are
// deterministic: they either pass on every run and every compiler, or never.
constexpr double kKsCritical = 1.95;

} // namespace

TEST_CASE("Random streams are bit-identical across compilers (golden values)", "[rng][golden]") {
    // Integer arithmetic plus one exact int-to-double conversion: identical on every platform.
    RandomStream rng(SeedKey{42}, 0);
    CHECK(rng.uniform() == 0x1.39d5e0a6efea9p-1);
    CHECK(rng.uniform() == 0x1.2bf50ad5742b8p-4);
    CHECK(rng.uniform() == 0x1.f9b6424ea774dp-1);
    CHECK(rng.uniform() == 0x1.071eb4dce89c1p-1);

    // Normals also go through log and sqrt: sqrt is correctly rounded everywhere, log is not
    // guaranteed to be, so exact equality is only asserted where the libm is shared.
    RandomStream normals(SeedKey{42, 1}, 7);
    const double z[3] = {normals.normal(), normals.normal(), normals.normal()};
#if defined(_MSC_VER)
    constexpr double tol = 4 * std::numeric_limits<double>::epsilon();
    CHECK(std::abs(z[0] - 0x1.4676932726ac8p-4) <= tol * 0x1.4676932726ac8p-4);
    CHECK(std::abs(z[1] - 0x1.271345ed25b4ap+0) <= tol * 0x1.271345ed25b4ap+0);
    CHECK(std::abs(z[2] + 0x1.e97a095291d53p-4) <= tol * 0x1.e97a095291d53p-4);
#else
    CHECK(z[0] == 0x1.4676932726ac8p-4);
    CHECK(z[1] == 0x1.271345ed25b4ap+0);
    CHECK(z[2] == -0x1.e97a095291d53p-4);
#endif
}

TEST_CASE("Uniforms lie strictly inside (0, 1) on a symmetric grid", "[rng]") {
    RandomStream rng(SeedKey{7}, 3);
    for (int i = 0; i < 100'000; ++i) {
        const double u = rng.uniform();
        REQUIRE(u > 0.0);
        REQUIRE(u < 1.0);
        // (k + 1/2) 2^-52: the reflected point is exact and on the same grid.
        const double k = u * 0x1p52 - 0.5;
        REQUIRE(k == std::floor(k));
        REQUIRE((1.0 - u) * 0x1p52 - 0.5 == std::floor((1.0 - u) * 0x1p52 - 0.5));
    }
}

TEST_CASE("Same key, same numbers; any other key, block or stream, different numbers", "[rng]") {
    const auto base = draw(SeedKey{11}, 0, 1000, false);
    CHECK(draw(SeedKey{11}, 0, 1000, false) == base);
    CHECK(draw(SeedKey{12}, 0, 1000, false) != base);
    CHECK(draw(SeedKey{11, 1}, 0, 1000, false) != base);
    CHECK(draw(SeedKey{11}, 1, 1000, false) != base);
    // Seeds differing only in their high 32 bits map to different Philox keys.
    CHECK(draw(SeedKey{11 + (std::uint64_t{1} << 32)}, 0, 1000, false) != base);
}

TEST_CASE("Uniforms: moments, Kolmogorov-Smirnov and independence", "[rng][statistical]") {
    const auto u = draw(SeedKey{2026}, 0, kSamples, false);
    const double n = static_cast<double>(kSamples);
    const double m = mean(u);
    CHECK(std::abs(m - 0.5) <= 4.0 * std::sqrt(1.0 / 12.0 / n));
    CHECK(std::abs(central_moment(u, 0.5, 2) - 1.0 / 12.0) <= 4.0 * std::sqrt(1.0 / 180.0 / n));
    CHECK(ks_statistic(u, [](double x) { return x; }) < kKsCritical);

    // Lag-1 serial correlation, and correlation between two blocks of the same seed.
    const std::vector<double> head(u.begin(), u.end() - 1), tail(u.begin() + 1, u.end());
    CHECK(std::abs(correlation(head, tail)) <= 4.0 / std::sqrt(n));
    const auto other_block = draw(SeedKey{2026}, 1, kSamples, false);
    CHECK(std::abs(correlation(u, other_block)) <= 4.0 / std::sqrt(n));
}

TEST_CASE("Normals: moments and Kolmogorov-Smirnov", "[rng][statistical]") {
    const auto z = draw(SeedKey{2026}, 5, kSamples, true);
    const double n = static_cast<double>(kSamples);
    const double m = mean(z);
    const double var = central_moment(z, m, 2);
    CHECK(std::abs(m) <= 4.0 / std::sqrt(n));
    CHECK(std::abs(var - 1.0) <= 4.0 * std::sqrt(2.0 / n));
    CHECK(std::abs(central_moment(z, m, 3) / std::pow(var, 1.5)) <= 4.0 * std::sqrt(6.0 / n));
    CHECK(std::abs(central_moment(z, m, 4) / (var * var) - 3.0) <= 4.0 * std::sqrt(24.0 / n));
    CHECK(ks_statistic(z, [](double x) { return norm_cdf(x); }) < kKsCritical);
}

TEST_CASE("AS241 inverse normal matches a 60-digit reference", "[rng][icdf]") {
    struct Point {
        double p, x;
    };
    // Generated by tools/normal_icdf_reference.py (bisection at 60 digits, independent of AS241).
    const Point points[] = {
        {1e-300, -37.047096299361199237},
        {1e-100, -21.273453560965324294},
        {1e-20, -9.2623400897984075796},
        {1e-11, -6.7060231554951362961},
        {2e-11, -6.6040775904056342776},
        {1e-10, -6.3613409024040561991},
        {1e-05, -4.2648907939228246102},
        {0.001, -3.0902323061678135354},
        {0.02, -2.0537489106318230443},
        {0.075, -1.4395314709384559349},
        {0.07500000000000001, -1.4395314709384558369},
        {0.1, -1.2815515655446004353},
        {0.3, -0.52440051270804081597},
        {0.45, -0.12566134685507400616},
        {0.5, 0.0},
        {0.55, 0.12566134685507414641},
        {0.8, 0.8416212335729143638},
        {0.925, 1.4395314709384562291},
        {0.975, 1.9599639845400538556},
        {0.999, 3.0902323061678132778},
        {0.9999999999, 6.3613408896974218642},
    };
    for (const auto& pt : points) {
        INFO("p=" << pt.p);
        const double x = norm_icdf(pt.p);
        if (pt.x == 0.0) {
            CHECK(x == 0.0);
        } else {
            CHECK(std::abs(x - pt.x) <= 1e-15 * std::abs(pt.x));
        }
    }
}

TEST_CASE("AS241 inverse normal: monotone, antisymmetric on the grid, inverts the CDF", "[rng][icdf]") {
    CHECK(norm_icdf(0.0) == -std::numeric_limits<double>::infinity());
    CHECK(norm_icdf(1.0) == std::numeric_limits<double>::infinity());
    CHECK(std::isnan(norm_icdf(-0.1)));
    CHECK(std::isnan(norm_icdf(1.1)));
    CHECK(std::isnan(norm_icdf(std::numeric_limits<double>::quiet_NaN())));

    RandomStream rng(SeedKey{99}, 0);
    for (int i = 0; i < 100'000; ++i) {
        const double u = rng.uniform();
        // Exact antisymmetry is what makes antithetic variates exact.
        REQUIRE(norm_icdf(1.0 - u) == -norm_icdf(u));
    }
    double previous = -std::numeric_limits<double>::infinity();
    for (double p = 1e-12; p < 1.0; p += 1e-4) {
        const double x = norm_icdf(p);
        REQUIRE(x > previous);
        previous = x;
        // Compare tail probabilities (1 - p is exact for p >= 1/2), where the information is.
        const double tail = std::min(p, 1.0 - p);
        const double tail_of_x = x < 0.0 ? norm_cdf(x) : norm_cdf(-x);
        REQUIRE(std::abs(tail_of_x - tail) <= 1e-13 * tail);
    }
}
