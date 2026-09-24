#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "riskengine/core/rng/owen_scramble.hpp"
#include "riskengine/core/rng/sobol.hpp"

using namespace riskengine;

namespace {

// SciPy 1.17 scipy.stats.qmc.Sobol(d, scramble=False, bits=32), times 2^32: an independent
// implementation of the same Joe-Kuo direction numbers, in the same Gray-code order.
constexpr std::uint32_t kFirstPoints[16][8] = {
    {0x00000000u, 0x00000000u, 0x00000000u, 0x00000000u, 0x00000000u, 0x00000000u, 0x00000000u, 0x00000000u},
    {0x80000000u, 0x80000000u, 0x80000000u, 0x80000000u, 0x80000000u, 0x80000000u, 0x80000000u, 0x80000000u},
    {0xc0000000u, 0x40000000u, 0x40000000u, 0x40000000u, 0xc0000000u, 0xc0000000u, 0x40000000u, 0xc0000000u},
    {0x40000000u, 0xc0000000u, 0xc0000000u, 0xc0000000u, 0x40000000u, 0x40000000u, 0xc0000000u, 0x40000000u},
    {0x60000000u, 0x60000000u, 0xa0000000u, 0xe0000000u, 0x60000000u, 0x20000000u, 0x60000000u, 0xe0000000u},
    {0xe0000000u, 0xe0000000u, 0x20000000u, 0x60000000u, 0xe0000000u, 0xa0000000u, 0xe0000000u, 0x60000000u},
    {0xa0000000u, 0x20000000u, 0xe0000000u, 0xa0000000u, 0xa0000000u, 0xe0000000u, 0x20000000u, 0x20000000u},
    {0x20000000u, 0xa0000000u, 0x60000000u, 0x20000000u, 0x20000000u, 0x60000000u, 0xa0000000u, 0xa0000000u},
    {0x30000000u, 0x50000000u, 0xf0000000u, 0x70000000u, 0x90000000u, 0x50000000u, 0x70000000u, 0xf0000000u},
    {0xb0000000u, 0xd0000000u, 0x70000000u, 0xf0000000u, 0x10000000u, 0xd0000000u, 0xf0000000u, 0x70000000u},
    {0xf0000000u, 0x10000000u, 0xb0000000u, 0x30000000u, 0x50000000u, 0x90000000u, 0x30000000u, 0x30000000u},
    {0x70000000u, 0x90000000u, 0x30000000u, 0xb0000000u, 0xd0000000u, 0x10000000u, 0xb0000000u, 0xb0000000u},
    {0x50000000u, 0x30000000u, 0x50000000u, 0x90000000u, 0xf0000000u, 0x70000000u, 0x10000000u, 0x10000000u},
    {0xd0000000u, 0xb0000000u, 0xd0000000u, 0x10000000u, 0x70000000u, 0xf0000000u, 0x90000000u, 0x90000000u},
    {0x90000000u, 0x70000000u, 0x10000000u, 0xd0000000u, 0x30000000u, 0xb0000000u, 0x50000000u, 0xd0000000u},
    {0x10000000u, 0xf0000000u, 0x90000000u, 0x50000000u, 0xb0000000u, 0x30000000u, 0xd0000000u, 0x50000000u},
};

struct FarPoint {
    std::uint64_t index;
    std::uint32_t x[4]; // dimensions 10, 100, 511, 1024 (1-based)
};
constexpr FarPoint kFarPoints[] = {
    {1000u, {0x11c00000u, 0x2fc00000u, 0x46400000u, 0xb6c00000u}},
    {123456u, {0x3fc28000u, 0xecff8000u, 0x32e98000u, 0x726c8000u}},
    {1048583u, {0x63c77800u, 0x63f7d800u, 0x3e799800u, 0xffda8800u}},
};
constexpr std::uint32_t kFarDims[4] = {9, 99, 510, 1023}; // 0-based

} // namespace

TEST_CASE("Sobol matches SciPy's independent implementation", "[sobol]") {
    const Sobol sobol(8);
    SobolCursor cursor(sobol, 0);
    for (std::uint64_t i = 0; i < 16; ++i) {
        for (std::uint32_t j = 0; j < 8; ++j) {
            CHECK(sobol.point(i, j) == kFirstPoints[i][j]);
            CHECK(cursor.point()[j] == kFirstPoints[i][j]);
        }
        cursor.advance();
    }

    const Sobol wide(sobol_data::kMaxDimension);
    for (const auto& far : kFarPoints) {
        INFO("index " << far.index);
        SobolCursor before(wide, far.index - 1);
        before.advance();
        for (int d = 0; d < 4; ++d) {
            CHECK(wide.point(far.index, kFarDims[d]) == far.x[d]);
            CHECK(before.point()[kFarDims[d]] == far.x[d]);
        }
    }
}

TEST_CASE("Sobol cursor agrees with random access over long runs", "[sobol]") {
    const Sobol sobol(64);
    SobolCursor cursor(sobol, 777);
    for (std::uint64_t i = 777; i < 777 + 5000; ++i) {
        for (std::uint32_t j = 0; j < 64; j += 7) REQUIRE(cursor.point()[j] == sobol.point(i, j));
        cursor.advance();
    }
}

TEST_CASE("Every dimension is stratified: 2^m points, one per interval of width 2^-m", "[sobol]") {
    // The (0, m, 1)-net property of each one-dimensional projection, for m = 1..12.
    const Sobol sobol(sobol_data::kMaxDimension);
    for (int m = 1; m <= 12; ++m)
        for (std::uint32_t j = 0; j < sobol.dimensions(); j += 37) {
            std::vector<int> hits(std::size_t{1} << m, 0);
            for (std::uint64_t i = 0; i < (std::uint64_t{1} << m); ++i) ++hits[sobol.point(i, j) >> (32 - m)];
            for (int h : hits) REQUIRE(h == 1);
        }
}

TEST_CASE("Bit reversal is an involution with the expected images", "[owen]") {
    CHECK(reverse_bits(1u) == 0x80000000u);
    CHECK(reverse_bits(0x0000000fu) == 0xf0000000u);
    CHECK(reverse_bits(0x12345678u) == 0x1e6a2c48u);
    for (std::uint32_t x : {0u, 7u, 0xdeadbeefu, 0xffffffffu}) CHECK(reverse_bits(reverse_bits(x)) == x);
}

TEST_CASE("Owen scrambling preserves the net structure", "[owen]") {
    const Sobol sobol(64);
    for (std::uint32_t replication = 0; replication < 4; ++replication) {
        const auto seeds = scramble_seeds(SeedKey{42}, replication, sobol.dimensions());
        // One-dimensional projections stay (0, m, 1)-nets.
        for (int m = 1; m <= 12; ++m)
            for (std::uint32_t j = 0; j < sobol.dimensions(); j += 5) {
                std::vector<int> hits(std::size_t{1} << m, 0);
                for (std::uint64_t i = 0; i < (std::uint64_t{1} << m); ++i)
                    ++hits[owen_scramble(sobol.point(i, j), seeds[j]) >> (32 - m)];
                for (int h : hits) REQUIRE(h == 1);
            }
        // The first two dimensions stay a (0, m, 2)-net: every elementary box 2^-a x 2^-(m-a) of
        // area 2^-m holds exactly one of the first 2^m points.
        constexpr int m = 10;
        for (int a = 0; a <= m; ++a) {
            std::vector<int> hits(std::size_t{1} << m, 0);
            for (std::uint64_t i = 0; i < (std::uint64_t{1} << m); ++i) {
                const std::uint32_t x = owen_scramble(sobol.point(i, 0), seeds[0]);
                const std::uint32_t y = owen_scramble(sobol.point(i, 1), seeds[1]);
                const std::uint32_t bx = a == 0 ? 0 : x >> (32 - a);
                const std::uint32_t by = a == m ? 0 : y >> (32 - (m - a));
                ++hits[(std::size_t{bx} << (m - a)) | by];
            }
            for (int h : hits) REQUIRE(h == 1);
        }
    }
}

TEST_CASE("A scrambled point is uniform over seeds, and seeds differ", "[owen]") {
    // Randomization makes every point uniform on (0, 1), which is what makes the estimator
    // unbiased. Kolmogorov-Smirnov over 100,000 seeds for a fixed input point.
    constexpr int n = 100'000;
    RandomStream rng(SeedKey{3}, 0);
    for (std::uint32_t x : {0u, 0x80000000u, 0x12345678u}) {
        std::vector<double> u(n);
        for (double& v : u) v = to_open_unit(owen_scramble(x, static_cast<std::uint32_t>(rng.uniform() * 0x1p32)));
        std::sort(u.begin(), u.end());
        double d = 0.0;
        for (int i = 0; i < n; ++i) d = std::max({d, u[i] - double(i) / n, double(i + 1) / n - u[i]});
        CHECK(std::sqrt(double(n)) * d < 1.95);
    }
    const auto a = scramble_seeds(SeedKey{1}, 0, 16), b = scramble_seeds(SeedKey{1}, 1, 16), c = scramble_seeds(SeedKey{2}, 0, 16);
    CHECK(a != b);
    CHECK(a != c);
    CHECK(a == scramble_seeds(SeedKey{1}, 0, 16));
    CHECK(to_open_unit(0u) > 0.0);
    CHECK(to_open_unit(0xffffffffu) < 1.0);
}

TEST_CASE("The last 32-bit Sobol index is reachable by random access and by the cursor", "[sobol]") {
    // Index 2^32 - 1 has gray code bit 31 set: every direction number is used, none beyond.
    const Sobol sobol(16);
    constexpr std::uint64_t last = (std::uint64_t{1} << 32) - 1;
    SobolCursor cursor(sobol, last - 1);
    cursor.advance();
    for (std::uint32_t j = 0; j < 16; ++j) CHECK(cursor.point()[j] == sobol.point(last, j));
    // van der Corput: gray(2^32 - 1) = 2^31 uses only direction number 31, which is 1 (= 2^-32).
    CHECK(sobol.point(last, 0) == 1u);
}
