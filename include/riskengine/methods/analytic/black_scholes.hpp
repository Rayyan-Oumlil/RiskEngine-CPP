#pragma once

#include <cassert>
#include <cmath>

#include "riskengine/core/estimate.hpp"
#include "riskengine/core/market.hpp"
#include "riskengine/core/normal.hpp"
#include "riskengine/greeks/greeks.hpp"
#include "riskengine/payoffs/vanilla.hpp"

namespace riskengine {

namespace detail {

// Everything the price and the Greeks share. N(x) and N(-x) are both computed through erfc
// rather than as 1 - N(x), so that deep in-the-money and out-of-the-money values keep full
// relative precision.
struct BsTerms {
    double df_r;          // e^{-rT}
    double df_q;          // e^{-qT}
    double n_d1, n_md1;   // N(d1), N(-d1)
    double n_d2, n_md2;   // N(d2), N(-d2)
    double pdf_d1;        // n(d1)
    double pdf_d2;        // n(d2)
    double sigma_sqrt_t;
    bool degenerate;      // sigma * sqrt(T) == 0: no diffusion left, price is discounted forward intrinsic
};

inline BsTerms bs_terms(double s, double k, double r, double q, double sigma, double t) {
    assert(s > 0.0 && k > 0.0 && sigma >= 0.0 && t >= 0.0);
    BsTerms b{};
    b.df_r = std::exp(-r * t);
    b.df_q = std::exp(-q * t);
    b.sigma_sqrt_t = sigma * std::sqrt(t);
    b.degenerate = !(b.sigma_sqrt_t > 0.0); // also catches underflow of sigma * sqrt(T)

    if (b.degenerate) {
        // Limit sigma*sqrt(T) -> 0: N(d1), N(d2) become a step in the forward moneyness.
        // At the money forward the limit of N(d1) is exactly 1/2.
        const double fwd = s * b.df_q;
        const double strike = k * b.df_r;
        const double step = fwd > strike ? 1.0 : (fwd < strike ? 0.0 : 0.5);
        b.n_d1 = b.n_d2 = step;
        b.n_md1 = b.n_md2 = 1.0 - step;
        b.pdf_d1 = b.pdf_d2 = 0.0;
        return b;
    }

    const double d1 = (std::log(s / k) + (r - q) * t) / b.sigma_sqrt_t + 0.5 * b.sigma_sqrt_t;
    const double d2 = d1 - b.sigma_sqrt_t;
    b.n_d1 = norm_cdf(d1);
    b.n_md1 = norm_cdf(-d1);
    b.n_d2 = norm_cdf(d2);
    b.n_md2 = norm_cdf(-d2);
    b.pdf_d1 = norm_pdf(d1);
    b.pdf_d2 = norm_pdf(d2);
    return b;
}

inline double bs_price(OptionType type, double s, double k, double r, double q, double sigma, double t) {
    const BsTerms b = bs_terms(s, k, r, q, sigma, t);
    return type == OptionType::Call ? s * b.df_q * b.n_d1 - k * b.df_r * b.n_d2
                                    : k * b.df_r * b.n_md2 - s * b.df_q * b.n_md1;
}

} // namespace detail

// Closed-form Black-Scholes-Merton price of a European option with a continuous dividend yield.
// Preconditions: spot > 0, strike > 0, vol >= 0, maturity >= 0. A zero vol or zero maturity is
// handled explicitly (discounted forward intrinsic value), never through a division by zero.
inline double black_scholes_price(const VanillaOption& o, const MarketState& m) {
    return detail::bs_price(o.type, m.spot.value, o.strike.value, m.rate.value, m.div.value,
                            m.vol.value, o.maturity.value);
}

// Closed-form Greeks, in the raw units of docs/conventions.md. In the degenerate case
// (sigma * sqrt(T) == 0) gamma and vega are 0 and delta, theta, rho are those of the
// discounted forward intrinsic value.
inline Greeks black_scholes_greeks(const VanillaOption& o, const MarketState& m) {
    const double s = m.spot.value, k = o.strike.value, r = m.rate.value, q = m.div.value;
    const double sigma = m.vol.value, t = o.maturity.value;
    const detail::BsTerms b = detail::bs_terms(s, k, r, q, sigma, t);

    const double gamma = b.degenerate ? 0.0 : b.df_q * b.pdf_d1 / (s * b.sigma_sqrt_t);
    const double vega = s * b.df_q * b.pdf_d1 * std::sqrt(t);
    const double time_decay = b.degenerate ? 0.0 : -s * b.df_q * b.pdf_d1 * sigma / (2.0 * std::sqrt(t));

    if (o.type == OptionType::Call) {
        return Greeks{
            .delta = b.df_q * b.n_d1,
            .gamma = gamma,
            .vega = vega,
            .theta = time_decay - r * k * b.df_r * b.n_d2 + q * s * b.df_q * b.n_d1,
            .rho = k * t * b.df_r * b.n_d2,
        };
    }
    return Greeks{
        .delta = -b.df_q * b.n_md1,
        .gamma = gamma,
        .vega = vega,
        .theta = time_decay + r * k * b.df_r * b.n_md2 - q * s * b.df_q * b.n_md1,
        .rho = -k * t * b.df_r * b.n_md2,
    };
}

// Pricer over the closed form. Satisfies the Pricer concept.
struct BlackScholes {
    VanillaOption option;

    Estimate price(const MarketState& m) const { return Estimate{black_scholes_price(option, m)}; }
    Greeks greeks(const MarketState& m) const { return black_scholes_greeks(option, m); }
};

} // namespace riskengine
