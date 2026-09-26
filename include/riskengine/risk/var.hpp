#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "riskengine/core/normal.hpp"
#include "riskengine/core/rng/normal_icdf.hpp"
#include "riskengine/core/rng/random_stream.hpp"

// Value-at-Risk and Expected Shortfall. Conventions (docs/model_risk_report.md 7):
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

struct Interval {
    double low, high;
};

struct RiskMeasuresWithCI {
    RiskMeasures estimate;
    Interval var, es;
};

// Percentile-bootstrap confidence intervals for the empirical VaR and ES: `resamples` samples of
// size N drawn with replacement (indices from RandomStream(key, b) for resample b, so the interval
// is reproducible), and the (1 - level)/2 and (1 + level)/2 quantiles of the resampled measures.
// In the tail the quantile estimator is noisy (its variance is alpha (1 - alpha) / (N f(q)^2)),
// which is why every VaR in the report carries such an interval.
inline RiskMeasuresWithCI bootstrap_var_es(const std::vector<double>& losses, double alpha, std::uint32_t resamples,
                                           SeedKey key, double level = 0.95) {
    assert(resamples >= 10 && level > 0.0 && level < 1.0);
    const std::size_t n = losses.size();
    std::vector<double> vars, ess, sample(n);
    for (std::uint32_t b = 0; b < resamples; ++b) {
        RandomStream rng(key, b);
        for (double& x : sample)
            x = losses[std::min(n - 1, static_cast<std::size_t>(rng.uniform() * static_cast<double>(n)))];
        const RiskMeasures m = empirical_var_es(sample, alpha);
        vars.push_back(m.var);
        ess.push_back(m.es);
    }
    const auto percentile = [&](std::vector<double> v, double p) {
        const auto k = static_cast<std::size_t>(std::floor(p * static_cast<double>(v.size() - 1)));
        std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(k), v.end());
        return v[k];
    };
    const double lo = 0.5 * (1.0 - level), hi = 0.5 * (1.0 + level);
    return {empirical_var_es(losses, alpha), {percentile(vars, lo), percentile(vars, hi)},
            {percentile(ess, lo), percentile(ess, hi)}};
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
