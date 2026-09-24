#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <span>

// VaR backtests (report 7.4). An exception is a day whose realized loss exceeds that day's VaR.
namespace riskengine {

struct TestResult {
    double statistic; // likelihood-ratio statistic
    double p_value;   // probability of a statistic at least this large if the model is right
};

namespace detail {

// x log(y), with the convention 0 log 0 = 0.
inline double xlogy(double x, double y) { return x == 0.0 ? 0.0 : x * std::log(y); }

// Survival functions of chi-square with 1 and 2 degrees of freedom.
inline double chi2_1_survival(double x) { return std::erfc(std::sqrt(0.5 * std::max(0.0, x))); }
inline double chi2_2_survival(double x) { return std::exp(-0.5 * std::max(0.0, x)); }

} // namespace detail

// Kupiec (1995) proportion-of-failures test of unconditional coverage: are x exceptions in n days
// consistent with an exception probability p = 1 - VaR level? LR ~ chi-square(1) under the model.
inline TestResult kupiec_pof(std::size_t n, std::size_t x, double p) {
    assert(n > 0 && x <= n && p > 0.0 && p < 1.0);
    const double nn = static_cast<double>(n), xx = static_cast<double>(x), phat = xx / nn;
    const double lr = -2.0 * (detail::xlogy(nn - xx, 1.0 - p) + detail::xlogy(xx, p) -
                              detail::xlogy(nn - xx, 1.0 - phat) - detail::xlogy(xx, phat));
    return {lr, detail::chi2_1_survival(lr)};
}

// Christoffersen (1998) test of independence: does an exception make the next day's exception more
// likely? Compares a first-order Markov chain for the exception indicator with an i.i.d. one.
// LR ~ chi-square(1) under independence. Exceptions clustered in crises reject it.
inline TestResult christoffersen_independence(std::span<const bool> exceptions) {
    double n00 = 0, n01 = 0, n10 = 0, n11 = 0;
    for (std::size_t t = 1; t < exceptions.size(); ++t) {
        const bool prev = exceptions[t - 1], cur = exceptions[t];
        (prev ? (cur ? n11 : n10) : (cur ? n01 : n00)) += 1.0;
    }
    const double pi0 = n01 / std::max(1.0, n00 + n01), pi1 = n11 / std::max(1.0, n10 + n11);
    const double pi = (n01 + n11) / std::max(1.0, n00 + n01 + n10 + n11);
    const double lr = -2.0 * (detail::xlogy(n00 + n10, 1.0 - pi) + detail::xlogy(n01 + n11, pi) -
                              detail::xlogy(n00, 1.0 - pi0) - detail::xlogy(n01, pi0) - detail::xlogy(n10, 1.0 - pi1) -
                              detail::xlogy(n11, pi1));
    return {lr, detail::chi2_1_survival(lr)};
}

// Conditional coverage: POF + independence, chi-square(2) under the model.
inline TestResult christoffersen_conditional_coverage(std::span<const bool> exceptions, double p) {
    std::size_t x = 0;
    for (bool e : exceptions) x += e ? 1 : 0;
    const double lr = kupiec_pof(exceptions.size(), x, p).statistic + christoffersen_independence(exceptions).statistic;
    return {lr, detail::chi2_2_survival(lr)};
}

// Basel traffic light for a 99 % VaR over 250 days: green up to 4 exceptions, yellow 5 to 9, red 10+.
enum class BaselZone { Green, Yellow, Red };

inline BaselZone basel_zone(std::size_t exceptions_in_250_days) {
    if (exceptions_in_250_days <= 4) return BaselZone::Green;
    if (exceptions_in_250_days <= 9) return BaselZone::Yellow;
    return BaselZone::Red;
}

inline const char* to_string(BaselZone z) {
    return z == BaselZone::Green ? "green" : (z == BaselZone::Yellow ? "yellow" : "red");
}

} // namespace riskengine
