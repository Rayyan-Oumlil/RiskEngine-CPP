#pragma once

#include <cstdint>
#include <vector>

#include "riskengine/core/rng/philox.hpp"
#include "riskengine/core/rng/random_stream.hpp"

namespace riskengine {

inline std::uint32_t reverse_bits(std::uint32_t x) {
    x = ((x >> 1) & 0x55555555u) | ((x & 0x55555555u) << 1);
    x = ((x >> 2) & 0x33333333u) | ((x & 0x33333333u) << 2);
    x = ((x >> 4) & 0x0f0f0f0fu) | ((x & 0x0f0f0f0fu) << 4);
    x = ((x >> 8) & 0x00ff00ffu) | ((x & 0x00ff00ffu) << 8);
    return (x >> 16) | (x << 16);
}

// Nested uniform (Owen) scrambling of a 32-bit fraction, hash-based (Burley 2020, "Practical
// Hash-based Owen Scrambling", with the Laine-Karras-style hash). In the bit-reversed domain the
// hash only propagates information from low to high bits, so each output digit depends only on
// the input digits above it: a random permutation of every sub-interval, as Owen's scrambling
// requires. It preserves the net structure of the Sobol points (their accuracy) while making every
// point uniformly distributed (unbiased estimates, and an error bar from independent seeds).
inline std::uint32_t owen_scramble(std::uint32_t x, std::uint32_t seed) {
    x = reverse_bits(x);
    x ^= x * 0x3d20adeau;
    x += seed;
    x *= (seed >> 16) | 1u;
    x ^= x * 0x05526c56u;
    x ^= x * 0x53a22864u;
    return reverse_bits(x);
}

// A 32-bit fraction as a double in (0, 1): the midpoint (x + 1/2) 2^-32, never 0 or 1.
inline double to_open_unit(std::uint32_t x) { return (static_cast<double>(x) + 0.5) * 0x1p-32; }

// Independent scrambling seeds for replication r of a `dims`-dimensional point set, drawn from
// Philox on the caller's key with the reserved scramble stream: seed of dimension j is word 0 of
// Philox(key.seed, counter = (j, r, key.stream | kScrambleStreamBit)).
inline std::vector<std::uint32_t> scramble_seeds(SeedKey key, std::uint32_t replication, std::uint32_t dims) {
    const philox::Key k{static_cast<std::uint32_t>(key.seed), static_cast<std::uint32_t>(key.seed >> 32)};
    std::vector<std::uint32_t> seeds(dims);
    for (std::uint32_t j = 0; j < dims; ++j)
        seeds[j] = philox::generate({j, 0, replication, key.stream | kScrambleStreamBit}, k)[0];
    return seeds;
}

} // namespace riskengine
