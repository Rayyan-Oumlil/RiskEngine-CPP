#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>

#include "riskengine/methods/analytic/black_scholes.hpp"

namespace riskengine {

// Discretely monitored geometric-average Asian option, fixings at t_i = i T / n (i = 1..n), payoff
// max(0, phi (G - K)) with G = (prod S_{t_i})^{1/n}. Under GBM, ln G is normal with
//     mean     ln S + (r - q - sigma^2/2) T (n + 1) / (2n)
//     variance sigma^2 T (n + 1)(2n + 1) / (6 n^2),
// so the price is Black-76 on the forward F = E[G] with total variance sigma_G^2. It is written as
// a Black-Scholes price with spot F, dividend yield r and vol sigma_G / sqrt(T), which reduces to
// the plain Black-Scholes price for n = 1. Used as the control variate for arithmetic Asians.
inline double geometric_asian_price(OptionType type, Strike strike, Maturity maturity, std::uint32_t fixings,
                                    const MarketState& m) {
    assert(fixings >= 1);
    const double n = static_cast<double>(fixings);
    const double t = maturity.value, sigma = m.vol.value, r = m.rate.value;
    const double var_g = sigma * sigma * t * (n + 1.0) * (2.0 * n + 1.0) / (6.0 * n * n);
    const double mean_g = std::log(m.spot.value) + (r - m.div.value - 0.5 * sigma * sigma) * t * (n + 1.0) / (2.0 * n);
    const double forward = std::exp(mean_g + 0.5 * var_g);
    const double vol = t > 0.0 ? std::sqrt(var_g / t) : 0.0;
    return detail::bs_price(type, forward, strike.value, r, r, vol, t);
}

} // namespace riskengine
