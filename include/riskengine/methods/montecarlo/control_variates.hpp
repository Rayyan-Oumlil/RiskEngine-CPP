#pragma once

#include <cmath>
#include <concepts>
#include <cstdint>
#include <span>

#include "riskengine/core/market.hpp"
#include "riskengine/methods/analytic/geometric_asian.hpp"
#include "riskengine/payoffs/asian.hpp"

namespace riskengine {

// A control variate X is evaluated on the same fixings as the payoff Y and must have a known
// discounted expectation mu = E[e^{-rT} X]. The engine then averages Y - beta (X - mu), with the
// variance-minimizing beta = Cov(Y, X) / Var(X) estimated on an independent pilot run.
template <class C>
concept ControlVariate = requires(const C& c, std::span<const double> fixings, const MarketState& m, Maturity t,
                                  std::uint32_t steps) {
    { c(fixings) } -> std::convertible_to<double>;
    { c.expectation(m, t, steps) } -> std::convertible_to<double>;
};

struct NoControl {};

// X = S_T, with E[e^{-rT} S_T] = S e^{-qT} for any model in which the discounted price (with
// dividends reinvested) is a martingale. Strong for in-the-money calls, weak far out of the money.
struct TerminalSpotControl {
    double operator()(std::span<const double> fixings) const { return fixings.back(); }

    double expectation(const MarketState& m, Maturity t, std::uint32_t) const {
        return m.spot.value * std::exp(-m.div.value * t.value);
    }
};

// X = geometric-average Asian payoff on the same fixings, whose price is known in closed form under
// GBM. Correlation with the arithmetic Asian is typically above 0.99.
struct GeometricAsianControl {
    double strike;
    OptionType type;

    double operator()(std::span<const double> fixings) const { return GeometricAsianPayoff{strike, type}(fixings); }

    double expectation(const MarketState& m, Maturity t, std::uint32_t steps) const {
        return geometric_asian_price(type, Strike{strike}, t, steps, m);
    }
};

static_assert(ControlVariate<TerminalSpotControl> && ControlVariate<GeometricAsianControl>);

} // namespace riskengine
