#pragma once

#include <array>
#include <cstdint>

namespace riskengine {

// Philox 4x32-10 (Salmon, Moraes, Dror, Shaw 2011, "Parallel random numbers: as easy as 1, 2, 3").
//
// A counter-based generator: a keyed bijection from a 128-bit counter to 128 random bits. There is
// no state to advance, so any draw of any stream can be computed directly from (key, counter).
// That is what makes results independent of how work is split across threads. Only integer
// operations are involved, so the output is bit-identical on every compiler and platform.
namespace philox {

using Counter = std::array<std::uint32_t, 4>;
using Key = std::array<std::uint32_t, 2>;

inline constexpr std::uint32_t kMultiplier0 = 0xD2511F53u;
inline constexpr std::uint32_t kMultiplier1 = 0xCD9E8D57u;
inline constexpr std::uint32_t kWeyl0 = 0x9E3779B9u; // golden ratio
inline constexpr std::uint32_t kWeyl1 = 0xBB67AE85u; // sqrt(3) - 1
inline constexpr int kRounds = 10;

constexpr Counter round(const Counter& c, const Key& k) {
    const std::uint64_t p0 = std::uint64_t{kMultiplier0} * c[0];
    const std::uint64_t p1 = std::uint64_t{kMultiplier1} * c[2];
    const auto hi = [](std::uint64_t p) { return static_cast<std::uint32_t>(p >> 32); };
    const auto lo = [](std::uint64_t p) { return static_cast<std::uint32_t>(p); };
    return {hi(p1) ^ c[1] ^ k[0], lo(p1), hi(p0) ^ c[3] ^ k[1], lo(p0)};
}

constexpr Counter generate(Counter c, Key k) {
    for (int r = 0; r < kRounds; ++r) {
        if (r > 0) {
            k[0] += kWeyl0;
            k[1] += kWeyl1;
        }
        c = round(c, k);
    }
    return c;
}

} // namespace philox

} // namespace riskengine
