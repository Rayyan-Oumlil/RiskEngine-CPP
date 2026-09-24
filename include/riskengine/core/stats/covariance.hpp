#pragma once

#include <cstdint>

namespace riskengine {

// Streaming covariance of (x, y) with the pairwise merge of Chan, Golub and LeVeque: the bivariate
// counterpart of Welford. Used to estimate control-variate coefficients.
class Covariance {
public:
    void add(double x, double y) {
        ++n_;
        const double inv_n = 1.0 / static_cast<double>(n_);
        const double dx = x - mean_x_;
        mean_x_ += dx * inv_n;
        mean_y_ += (y - mean_y_) * inv_n;
        m2_x_ += dx * (x - mean_x_);
        c_xy_ += dx * (y - mean_y_);
    }

    void merge(const Covariance& other) {
        if (other.n_ == 0) return;
        if (n_ == 0) {
            *this = other;
            return;
        }
        const double na = static_cast<double>(n_), nb = static_cast<double>(other.n_);
        const double n = na + nb;
        const double dx = other.mean_x_ - mean_x_, dy = other.mean_y_ - mean_y_;
        mean_x_ += dx * (nb / n);
        mean_y_ += dy * (nb / n);
        m2_x_ += other.m2_x_ + dx * dx * (na * nb / n);
        c_xy_ += other.c_xy_ + dx * dy * (na * nb / n);
        n_ += other.n_;
    }

    std::uint64_t count() const { return n_; }
    double mean_x() const { return mean_x_; }
    double mean_y() const { return mean_y_; }
    // Unbiased; 0 with fewer than two observations.
    double variance_x() const { return n_ > 1 ? m2_x_ / static_cast<double>(n_ - 1) : 0.0; }
    double covariance() const { return n_ > 1 ? c_xy_ / static_cast<double>(n_ - 1) : 0.0; }

private:
    std::uint64_t n_ = 0;
    double mean_x_ = 0.0, mean_y_ = 0.0;
    double m2_x_ = 0.0, c_xy_ = 0.0;
};

} // namespace riskengine
