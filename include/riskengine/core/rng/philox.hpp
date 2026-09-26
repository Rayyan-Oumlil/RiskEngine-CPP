#pragma once

#include <array>
#include <cstdint>

#if defined(__AVX2__)
#include <immintrin.h>
#endif

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

#if defined(__AVX2__)
// Four Philox calls at once. c[w] holds word w of four counters, one per 64-bit lane, zero-extended
// (the layout vpmuludq wants: it multiplies the low 32 bits of each 64-bit lane into a full 64-bit
// product). Only integer multiply, shift, and and xor are involved, so each lane is bit-identical
// to generate() on that lane's counter.
inline void generate_x4(__m256i (&c)[4], Key k) {
    const __m256i m0 = _mm256_set1_epi64x(kMultiplier0);
    const __m256i m1 = _mm256_set1_epi64x(kMultiplier1);
    const __m256i low32 = _mm256_set1_epi64x(0xFFFFFFFF);
    for (int r = 0; r < kRounds; ++r) {
        if (r > 0) {
            k[0] += kWeyl0;
            k[1] += kWeyl1;
        }
        const __m256i p0 = _mm256_mul_epu32(c[0], m0);
        const __m256i p1 = _mm256_mul_epu32(c[2], m1);
        const __m256i k0 = _mm256_set1_epi64x(k[0]);
        const __m256i k1 = _mm256_set1_epi64x(k[1]);
        // round(): {hi(p1) ^ c1 ^ k0, lo(p1), hi(p0) ^ c3 ^ k1, lo(p0)}
        const __m256i n0 = _mm256_xor_si256(_mm256_xor_si256(_mm256_srli_epi64(p1, 32), c[1]), k0);
        const __m256i n2 = _mm256_xor_si256(_mm256_xor_si256(_mm256_srli_epi64(p0, 32), c[3]), k1);
        c[1] = _mm256_and_si256(p1, low32);
        c[3] = _mm256_and_si256(p0, low32);
        c[0] = n0;
        c[2] = n2;
    }
}
#endif

} // namespace philox

} // namespace riskengine
