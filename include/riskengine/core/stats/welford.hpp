#pragma once

#include <cmath>
#include <cstdint>

#include "riskengine/core/estimate.hpp"

namespace riskengine {

// Streaming mean and variance (Welford 1962), with the pairwise merge of Chan, Golub and LeVeque
// (1979). Numerically stable where the naive sum of squares cancels catastrophically.
class Welford {
public:
    void add(double x) {
        ++n_;
        const double delta = x - mean_;
        mean_ += delta / static_cast<double>(n_);
        m2_ += delta * (x - mean_);
    }

    void merge(const Welford& other) {
        if (other.n_ == 0) return;
        if (n_ == 0) {
            *this = other;
            return;
        }
        const double na = static_cast<double>(n_), nb = static_cast<double>(other.n_);
        const double n = na + nb;
        const double delta = other.mean_ - mean_;
        mean_ += delta * (nb / n);
        m2_ += other.m2_ + delta * delta * (na * nb / n);
        n_ += other.n_;
    }

    std::uint64_t count() const { return n_; }
    double mean() const { return mean_; }
    // Unbiased sample variance; 0 with fewer than two observations.
    double variance() const { return n_ > 1 ? m2_ / static_cast<double>(n_ - 1) : 0.0; }
    double std_error() const { return n_ > 1 ? std::sqrt(variance() / static_cast<double>(n_)) : 0.0; }

private:
    std::uint64_t n_ = 0;
    double mean_ = 0.0;
    double m2_ = 0.0;
};

inline Estimate to_estimate(const Welford& w) {
    return Estimate{w.mean(), w.std_error(), 0.0, w.count()};
}

} // namespace riskengine
