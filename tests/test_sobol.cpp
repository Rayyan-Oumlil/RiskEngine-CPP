#include <catch2/catch_test_macros.hpp>

#include <cstdint>
#include <vector>

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
