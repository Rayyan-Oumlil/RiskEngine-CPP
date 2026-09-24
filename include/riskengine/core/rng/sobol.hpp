#pragma once

#include <bit>
#include <cassert>
#include <cstdint>
#include <vector>

#include "riskengine/core/rng/sobol_directions.hpp"

namespace riskengine {

// Sobol low-discrepancy sequence (Bratley & Fox 1988; Joe & Kuo 2008 direction numbers), 32-bit
// precision, up to sobol_data::kMaxDimension dimensions and 2^32 points.
//
// Points are in Gray-code order (Antonov & Saleev 1979), the order used by SciPy and most
// libraries: coordinate j of point i is the XOR of the direction numbers V_j[k] over the set bits
// k of gray(i) = i ^ (i >> 1). Point i can therefore be computed directly (a block of paths can
// start anywhere), and consecutive points differ by a single XOR (SobolCursor).
class Sobol {
public:
    static constexpr int kBits = 32;

    explicit Sobol(std::uint32_t dimensions) : dimensions_(dimensions), v_(std::size_t{dimensions} * kBits) {
        assert(dimensions >= 1 && dimensions <= sobol_data::kMaxDimension);
        for (int k = 0; k < kBits; ++k) v_[k] = 1u << (kBits - 1 - k); // van der Corput
        for (std::uint32_t j = 1; j < dimensions; ++j) {
            const sobol_data::Direction& dir = sobol_data::kDirections[j - 1];
            const unsigned s = dir.degree;
            std::uint32_t* v = &v_[std::size_t{j} * kBits];
            for (unsigned k = 0; k < s && k < kBits; ++k) v[k] = dir.m[k] << (kBits - 1 - k);
            for (unsigned k = s; k < kBits; ++k) {
                std::uint32_t x = v[k - s] ^ (v[k - s] >> s);
                for (unsigned i = 1; i < s; ++i)
                    if ((dir.coefficients >> (s - 1 - i)) & 1u) x ^= v[k - i];
                v[k] = x;
            }
        }
    }

    std::uint32_t dimensions() const { return dimensions_; }

    // Coordinate `dim` of point `index`, as the 32-bit fraction x / 2^32.
    std::uint32_t point(std::uint64_t index, std::uint32_t dim) const {
        assert(index < (std::uint64_t{1} << kBits) && dim < dimensions_);
        const auto gray = static_cast<std::uint32_t>(index ^ (index >> 1));
        std::uint32_t x = 0;
        for (int k = 0; gray >> k; ++k)
            if ((gray >> k) & 1u) x ^= v_[std::size_t{dim} * kBits + k];
        return x;
    }

    std::uint32_t direction(std::uint32_t dim, int k) const { return v_[std::size_t{dim} * kBits + k]; }

private:
    std::uint32_t dimensions_;
    std::vector<std::uint32_t> v_; // v_[dim * 32 + k]: direction number k of dimension dim
};

// Walks the points start, start + 1, ... of a Sobol sequence, one XOR per coordinate per step.
class SobolCursor {
public:
    SobolCursor(const Sobol& sobol, std::uint64_t start) : sobol_(&sobol), index_(start), x_(sobol.dimensions()) {
        for (std::uint32_t j = 0; j < sobol.dimensions(); ++j) x_[j] = sobol.point(start, j);
    }

    // The current point's coordinates (32-bit fractions); valid until advance().
    const std::vector<std::uint32_t>& point() const { return x_; }
    std::uint64_t index() const { return index_; }

    void advance() {
        const int k = std::countr_zero(index_ + 1); // gray(i + 1) = gray(i) ^ (1 << k)
        for (std::uint32_t j = 0; j < x_.size(); ++j) x_[j] ^= sobol_->direction(j, k);
        ++index_;
    }

private:
    const Sobol* sobol_;
    std::uint64_t index_;
    std::vector<std::uint32_t> x_;
};

} // namespace riskengine
