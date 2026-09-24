#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

#include "riskengine/methods/analytic/black_scholes.hpp"

namespace riskengine {

enum class ImpliedVolStatus {
    Ok,
    ZeroTimeValue,  // price equals the discounted forward intrinsic value: vol = 0
    BelowIntrinsic, // price violates the lower no-arbitrage bound
    AboveMaximum,   // price violates the upper no-arbitrage bound
    NotConverged,
};

struct ImpliedVolResult {
    double vol;
    ImpliedVolStatus status;
    int iterations = 0;

    bool ok() const { return status == ImpliedVolStatus::Ok || status == ImpliedVolStatus::ZeroTimeValue; }
};

namespace detail {

// Brent's method (Brent 1973, "zero") on [a, b] with f(a) < 0 < f(b), both finite.
// Terminates on an x-tolerance of a few ulps, so the result is as precise as f allows.
template <class F>
std::pair<double, int> brent_root(F&& f, double a, double b, double fa, double fb, int max_iter) {
    constexpr double eps = std::numeric_limits<double>::epsilon();
    double c = a, fc = fa, d = b - a, e = d;
    for (int iter = 1; iter <= max_iter; ++iter) {
        if ((fb > 0.0) == (fc > 0.0)) {
            c = a;
            fc = fa;
            d = e = b - a;
        }
        if (std::abs(fc) < std::abs(fb)) {
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }
        const double tol = 2.0 * eps * std::abs(b);
        const double m = 0.5 * (c - b);
        if (std::abs(m) <= tol || fb == 0.0) return {b, iter};

        if (std::abs(e) >= tol && std::abs(fa) > std::abs(fb)) {
            // Inverse quadratic interpolation, or secant when only two points are distinct.
            double p, q;
            const double s = fb / fa;
            if (a == c) {
                p = 2.0 * m * s;
                q = 1.0 - s;
            } else {
                const double qa = fa / fc, r = fb / fc;
                p = s * (2.0 * m * qa * (qa - r) - (b - a) * (r - 1.0));
                q = (qa - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0.0) q = -q; else p = -p;
            if (2.0 * p < std::min(3.0 * m * q - std::abs(tol * q), std::abs(e * q))) {
                e = d;
                d = p / q;
            } else {
                d = m;
                e = m;
            }
        } else {
            d = m;
            e = m;
        }
        a = b;
        fa = fb;
        b += std::abs(d) > tol ? d : (m > 0.0 ? tol : -tol);
        fb = f(b);
    }
    return {b, -1};
}

} // namespace detail

// Implied volatility of a European option under Black-Scholes-Merton.
//
// Method: put-call parity maps the quote to the out-of-the-money option (whose price is pure
// time value), then Brent's method solves
// log(price_otm(sigma)) = log(target) on a bracket. Working in log price keeps the problem
// well conditioned in the deep wings, where the price spans hundreds of orders of magnitude.
// Newton is deliberately not used: it diverges where vega -> 0.
//
// Conditioning: an in-the-money quote carries its time value next to a much larger intrinsic
// value, so the time value is only known to ~eps * price in absolute terms. The recovered vol is
// then only as accurate as eps * price / vega allows; out-of-the-money quotes do not lose this.
//
// Jaeckel's "Let's Be Rational" (2015) remains the long-term target (see
// docs/riskengine_research.md 3.1); this bracketed solver is the robust fallback.
inline ImpliedVolResult implied_vol(const VanillaOption& o, double price, Spot spot, Rate rate, Rate div) {
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    const double s = spot.value, k = o.strike.value, r = rate.value, q = div.value, t = o.maturity.value;
    const double df_r = std::exp(-r * t), df_q = std::exp(-q * t);
    const double fwd_intrinsic = s * df_q - k * df_r; // C - P by parity

    // Out-of-the-money option in the forward sense, and its price implied by the quote.
    const OptionType otm_type = fwd_intrinsic > 0.0 ? OptionType::Put : OptionType::Call;
    double target = price;
    if (o.type == OptionType::Call && fwd_intrinsic > 0.0) target = price - fwd_intrinsic;
    if (o.type == OptionType::Put && fwd_intrinsic < 0.0) target = price + fwd_intrinsic;

    // Converting an in-the-money quote subtracts the intrinsic value, which costs a few ulps of
    // the terms involved: a time value inside that noise is indistinguishable from zero.
    const double rounding = 4.0 * std::numeric_limits<double>::epsilon() * (s * df_q + k * df_r);
    const double upper = otm_type == OptionType::Call ? s * df_q : k * df_r;
    if (!(target >= -rounding)) return {nan, ImpliedVolStatus::BelowIntrinsic};
    if (target <= 0.0) return {0.0, ImpliedVolStatus::ZeroTimeValue};
    if (target >= upper) return {nan, ImpliedVolStatus::AboveMaximum};
    if (t <= 0.0) return {nan, ImpliedVolStatus::AboveMaximum}; // expired option has no time value

    const double log_target = std::log(target);
    auto g = [&](double sigma) {
        return std::log(detail::bs_price(otm_type, s, k, r, q, sigma, t)) - log_target;
    };

    // Bracket: g is increasing in sigma, g(0+) = -inf (or very negative), g(inf) = log(upper/target) > 0.
    double hi = 0.2, g_hi = g(hi);
    for (int i = 0; g_hi <= 0.0; ++i) {
        if (i == 64) return {nan, ImpliedVolStatus::NotConverged};
        hi *= 2.0;
        g_hi = g(hi);
    }
    double lo = hi, g_lo = g_hi;
    for (int i = 0; g_lo > 0.0; ++i) {
        if (i == 2000) return {nan, ImpliedVolStatus::NotConverged};
        hi = lo;
        g_hi = g_lo;
        lo *= 0.5;
        g_lo = g(lo);
    }
    // The price may have underflowed to 0 at lo (g = -inf): bisect until both ends are finite.
    for (int i = 0; !std::isfinite(g_lo); ++i) {
        if (i == 200) return {nan, ImpliedVolStatus::NotConverged};
        const double mid = 0.5 * (lo + hi);
        const double g_mid = g(mid);
        if (g_mid > 0.0) {
            hi = mid;
            g_hi = g_mid;
        } else {
            lo = mid;
            g_lo = g_mid;
        }
    }
    if (g_hi == 0.0) return {hi, ImpliedVolStatus::Ok};
    if (g_lo == 0.0) return {lo, ImpliedVolStatus::Ok};

    const auto [vol, iterations] = detail::brent_root(g, lo, hi, g_lo, g_hi, 200);
    if (iterations < 0) return {vol, ImpliedVolStatus::NotConverged, 200};
    return {vol, ImpliedVolStatus::Ok, iterations};
}

} // namespace riskengine
