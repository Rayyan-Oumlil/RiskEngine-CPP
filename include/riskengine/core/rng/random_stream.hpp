#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "riskengine/core/simd.hpp"

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
    // spare included. With AVX2 or AVX-512F, 4 or 8 Philox calls (8 or 16 uniforms) run at a time.
    void uniforms(std::span<double> u) {
        std::size_t i = 0;
        if (has_spare_ && !u.empty()) u[i++] = uniform(); // hand out the pending spare first
#if defined(RISKENGINE_SIMD)
        using S = simd::Widest;
        for (; u.size() - i >= 2 * S::kLanes; i += 2 * S::kLanes) {
            fill_lanes<S>(u.data() + i);
            index_ += S::kLanes;
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
#if defined(RISKENGINE_SIMD)
    // Uniforms 2j and 2j + 1 of Philox calls j = index_ .. index_ + S::kLanes - 1, in order, into
    // out[0 .. 2 * S::kLanes).
    template <class S>
    void fill_lanes(double* out) const {
        constexpr int kLanes = S::kLanes;
        std::uint64_t lo[kLanes], hi[kLanes];
        for (int l = 0; l < kLanes; ++l) {
            const std::uint64_t j = index_ + static_cast<std::uint64_t>(l);
            lo[l] = j & 0xFFFFFFFFu;
            hi[l] = j >> 32;
        }
        typename S::Int c[4] = {S::from_lanes(lo), S::from_lanes(hi), S::set1(std::uint64_t{block_}),
                                S::set1(std::uint64_t{stream_})};
        philox::generate_lanes<S>(c, key_);
        // Per call: to_uniform(bits[0], bits[1]) then to_uniform(bits[2], bits[3]).
        S::store_interleaved(out, to_uniform_lanes<S>(c[0], c[1]), to_uniform_lanes<S>(c[2], c[3]));
    }

    // to_uniform on a vector. k < 2^52, so double(k) is exact by the standard trick: put k in the
    // mantissa of 2^52 and subtract 2^52. (k + 0.5) and the power-of-two scaling are exact too.
    template <class S>
    static typename S::Vec to_uniform_lanes(typename S::Int hi, typename S::Int lo) {
        const auto k = S::shr(S::bit_or(S::shl(hi, 32), lo), 12);
        const auto two52 = S::set1(0x1p52);
        const auto kd = S::sub(S::as_double(S::bit_or(k, S::as_int(two52))), two52);
        return S::mul(S::add(kd, S::set1(0.5)), S::set1(0x1p-52));
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

// Hands out one stream's normals, one or a span at a time, from a read-ahead buffer filled by
// RandomStream::normals. For consumers that need a single normal per path (a European payoff, a
// Greek estimator), which would otherwise never reach the batched, vectorized generator.
//
// Same sequence as calling normal() on the stream repeatedly: the buffer owns its stream and only
// draws further ahead on it. Nothing else may draw from that stream, which is why it is taken by
// value. Drawing ahead past the last value used wastes at most one buffer of work.
class NormalBuffer {
public:
    explicit NormalBuffer(RandomStream rng, std::size_t capacity = 512)
        : rng_(rng), buffer_(capacity), next_(capacity) {}

    double next() {
        if (next_ == buffer_.size()) refill();
        return buffer_[next_++];
    }

    void fill(std::span<double> z) {
        for (std::size_t i = 0; i < z.size();) {
            if (next_ == buffer_.size()) refill();
            const std::size_t n = std::min(z.size() - i, buffer_.size() - next_);
            std::copy_n(buffer_.begin() + static_cast<std::ptrdiff_t>(next_), n, z.begin() + static_cast<std::ptrdiff_t>(i));
            next_ += n;
            i += n;
        }
    }

private:
    void refill() {
        rng_.normals(buffer_);
        next_ = 0;
    }

    RandomStream rng_;
    std::vector<double> buffer_;
    std::size_t next_;
};

} // namespace riskengine
