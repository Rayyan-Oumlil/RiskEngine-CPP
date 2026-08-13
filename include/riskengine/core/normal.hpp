#pragma once

#include <cmath>
#include <numbers>

namespace riskengine {

// Standard normal CDF via erfc: accurate in both tails. Never compute it as 1 - norm_cdf(-x).
inline double norm_cdf(double x) {
    return 0.5 * std::erfc(-x / std::numbers::sqrt2);
}

inline double norm_pdf(double x) {
    return std::exp(-0.5 * x * x) * (0.5 * std::numbers::inv_sqrtpi * std::numbers::sqrt2);
}

} // namespace riskengine
