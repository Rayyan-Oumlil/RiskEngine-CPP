#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

#include "riskengine/core/normal.hpp"
#include "riskengine/core/rng/normal_icdf.hpp"

// Value-at-Risk and Expected Shortfall. Conventions (docs/riskengine_research.md 6.3):
// - losses are positive (loss = -P&L), and VaR_alpha is a positive loss: 99 % VaR is exceeded on
//   1 % of days;
// - the empirical VaR is the ceil(alpha N)-th smallest loss (1-based), found in O(N) with
//   nth_element; the empirical ES is the mean of the losses from that one upwards, i.e. of the
//   worst N - ceil(alpha N) + 1 losses.
namespace riskengine {

struct RiskMeasures {
    double var;
    double es;
};

// Empirical VaR and ES of a loss sample (taken by value: it is partially reordered).
inline RiskMeasures empirical_var_es(std::vector<double> losses, double alpha) {
    assert(!losses.empty() && alpha > 0.0 && alpha < 1.0);
    const std::size_t n = losses.size();
    const auto k = static_cast<std::size_t>(std::ceil(alpha * static_cast<double>(n))) - 1; // 0-based
    std::nth_element(losses.begin(), losses.begin() + static_cast<std::ptrdiff_t>(k), losses.end());
    const double var = losses[k];
    double tail = 0.0;
    for (std::size_t i = k; i < n; ++i) tail += losses[i]; // everything from k on is >= var
    return {var, tail / static_cast<double>(n - k)};
}

// Normal loss distribution with mean mu and standard deviation sd:
// VaR = mu + sd z_alpha, ES = mu + sd n(z_alpha) / (1 - alpha).
inline RiskMeasures normal_var_es(double mu, double sd, double alpha) {
    const double z = norm_icdf(alpha);
    return {mu + sd * z, mu + sd * norm_pdf(z) / (1.0 - alpha)};
}

// Cornish-Fisher expansion of the alpha-quantile of a standardized distribution with skewness
// `skew` and excess kurtosis `excess_kurtosis`: corrects the normal quantile for the first two
// departures from normality. Improves on the normal quantile for moderately skewed distributions,
// but is not a distribution (it can even be non-monotone in alpha) and degrades in the far tail.
inline double cornish_fisher_quantile(double alpha, double skew, double excess_kurtosis) {
    const double z = norm_icdf(alpha);
    return z + (z * z - 1.0) * skew / 6.0 + (z * z * z - 3.0 * z) * excess_kurtosis / 24.0 -
           (2.0 * z * z * z - 5.0 * z) * skew * skew / 36.0;
}

} // namespace riskengine
