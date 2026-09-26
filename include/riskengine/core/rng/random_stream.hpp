#pragma once

#include <cstddef>
#include <cstdint>
#include <span>

#if defined(__AVX2__)
#include <immintrin.h>
#endif

#include "riskengine/core/rng/normal_icdf.hpp"
#include "riskengine/core/rng/philox.hpp"

namespace riskengine {

// Identifies a reproducible random experiment. A stochastic pricer is a pure function of
// (market, SeedKey): calling it twice with the same key after bumping the market gives exact
// common random numbers. `stream` separates independent uses under one seed; callers use streams
// below 2^30, the two top bits are reserved for internal sub-streams derived from a caller's key.
struct SeedKey {
    std::uint64_t seed;
    std::uint32_t stream = 0;
};

inline constexpr std::uint32_t kPilotStreamBit = 0x8000'0000u;    // control-variate pilot runs
inline constexpr std::uint32_t kScrambleStreamBit = 0x4000'0000u; // QMC scrambling seeds
inline constexpr std::uint32_t kReservedStreamBits = kPilotStreamBit | kScrambleStreamBit;

// The random sequence of one block of paths. Philox call j of block b uses key = seed and
// counter = (j, b, stream), and its 128 bits give uniforms 2j and 2j + 1 of the block. Blocks never
// share state, so they can run on any thread in any order.
class RandomStream {
public:
    RandomStream(SeedKey key, std::uint32_t block)
        : key_{static_cast<std::uint32_t>(key.seed), static_cast<std::uint32_t>(key.seed >> 32)},
          block_(block),
          stream_(key.stream) {}

    // Uniform on the open interval (0, 1): the midpoints of a 2^52 grid, (k + 1/2) 2^-52. The grid
    // is symmetric (1 - u is exact and on the grid) and never reaches 0 or 1.
    double uniform() {
        if (has_spare_) {
            has_spare_ = false;
            return spare_;
        }
        const philox::Counter bits = philox::generate(
            {static_cast<std::uint32_t>(index_), static_cast<std::uint32_t>(index_ >> 32), block_, stream_}, key_);
        ++index_;
        spare_ = to_uniform(bits[2], bits[3]);
        has_spare_ = true;
        return to_uniform(bits[0], bits[1]);
    }

    // Standard normal by inversion, so normal(1 - u) == -normal(u) exactly.
    double normal() { return norm_icdf(uniform()); }

    // Fills u with the next u.size() uniforms: bit-identical to calling uniform() that many times,
    // spare included. With AVX2, four Philox calls (eight uniforms) run at a time.
    void uniforms(std::span<double> u) {
        std::size_t i = 0;
        if (has_spare_ && !u.empty()) u[i++] = uniform(); // hand out the pending spare first
#if defined(__AVX2__)
        for (; u.size() - i >= 8; i += 8) {
            fill8(u.data() + i);
            index_ += 4;
        }
#endif
        for (; i < u.size(); ++i) u[i] = uniform(); // remainder; may leave a spare, as uniform() would
    }

    // Fills z with the next z.size() normals: bit-identical to calling normal() that many times
    // (same uniforms, same order), but draws and inverts them as batches, which vectorize with AVX2.
    void normals(std::span<double> z) {
        uniforms(z);
        norm_icdf(std::span<const double>(z), z);
    }

private:
#if defined(__AVX2__)
    // Uniforms 2j and 2j + 1 of Philox calls j = index_ .. index_ + 3, in order, into out[0..7].
    void fill8(double* out) const {
        const auto word = [](std::uint64_t v) { return static_cast<long long>(v & 0xFFFFFFFFu); };
        const std::uint64_t j = index_;
        __m256i c[4] = {
            _mm256_set_epi64x(word(j + 3), word(j + 2), word(j + 1), word(j)),                 // index, low
            _mm256_set_epi64x(word((j + 3) >> 32), word((j + 2) >> 32), word((j + 1) >> 32), word(j >> 32)),
            _mm256_set1_epi64x(block_),
            _mm256_set1_epi64x(stream_),
        };
        philox::generate_x4(c, key_);
        const __m256d first = to_uniform_x4(c[0], c[1]);  // to_uniform(bits[0], bits[1]) per call
        const __m256d second = to_uniform_x4(c[2], c[3]); // to_uniform(bits[2], bits[3]) per call
        // Interleave to call order: first0 second0 first1 second1 | first2 second2 first3 second3.
        const __m256d lo = _mm256_unpacklo_pd(first, second); // f0 s0 f2 s2
        const __m256d hi = _mm256_unpackhi_pd(first, second); // f1 s1 f3 s3
        _mm256_storeu_pd(out, _mm256_permute2f128_pd(lo, hi, 0x20));
        _mm256_storeu_pd(out + 4, _mm256_permute2f128_pd(lo, hi, 0x31));
    }

    // to_uniform on four lanes. k < 2^52, so double(k) is exact by the standard trick: put k in the
    // mantissa of 2^52 and subtract 2^52. (k + 0.5) and the power-of-two scaling are exact too.
    static __m256d to_uniform_x4(__m256i hi, __m256i lo) {
        const __m256i k = _mm256_srli_epi64(_mm256_or_si256(_mm256_slli_epi64(hi, 32), lo), 12);
        const __m256d two52 = _mm256_set1_pd(0x1p52);
        const __m256d kd = _mm256_sub_pd(_mm256_castsi256_pd(_mm256_or_si256(k, _mm256_castpd_si256(two52))), two52);
        return _mm256_mul_pd(_mm256_add_pd(kd, _mm256_set1_pd(0.5)), _mm256_set1_pd(0x1p-52));
    }
#endif

    static double to_uniform(std::uint32_t hi, std::uint32_t lo) {
        const std::uint64_t k = ((std::uint64_t{hi} << 32) | lo) >> 12; // top 52 bits
        return (static_cast<double>(k) + 0.5) * 0x1p-52;
    }

    philox::Key key_;
    std::uint32_t block_;
    std::uint32_t stream_;
    std::uint64_t index_ = 0;
    double spare_ = 0.0;
    bool has_spare_ = false;
};

} // namespace riskengine
