#pragma once

#include "riskengine/methods/analytic/black_scholes.hpp"

namespace riskengine {

// Cash-or-nothing digital paying 1 at maturity if the option finishes in the money:
// call = e^{-rT} N(d2), put = e^{-rT} N(-d2). Same preconditions and degenerate branch as
// black_scholes_price: at sigma * sqrt(T) == 0 the payout is a step in forward moneyness.
inline double digital_price(const VanillaOption& o, const MarketState& m) {
    const detail::BsTerms b = detail::bs_terms(m.spot.value, o.strike.value, m.rate.value, m.div.value,
                                               m.vol.value, o.maturity.value);
    return b.df_r * (o.type == OptionType::Call ? b.n_d2 : b.n_md2);
}

// dV/dS = +-e^{-rT} n(d2) / (S sigma sqrt(T)). In the degenerate case the price is a step, whose
// slope is 0 everywhere except exactly at the strike; 0 is returned.
inline double digital_delta(const VanillaOption& o, const MarketState& m) {
    const double s = m.spot.value;
    const detail::BsTerms b = detail::bs_terms(s, o.strike.value, m.rate.value, m.div.value, m.vol.value,
                                               o.maturity.value);
    if (b.degenerate) return 0.0;
    const double delta = b.df_r * b.pdf_d2 / (s * b.sigma_sqrt_t);
    return o.type == OptionType::Call ? delta : -delta;
}

} // namespace riskengine
