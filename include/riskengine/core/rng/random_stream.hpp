#pragma once

#include <cstdint>

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

private:
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
