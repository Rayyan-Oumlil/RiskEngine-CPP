#pragma once

#include <algorithm>
#include <vector>

#include "riskengine/greeks/greeks.hpp"
#include "riskengine/methods/analytic/black_scholes.hpp"

namespace riskengine {

// A portfolio of European options on one underlying plus a position in the underlying itself,
// valued in closed form. Maturities are measured from today; value(m, elapsed) ages every option
// by `elapsed` years (an expired option is worth its intrinsic value).
struct Portfolio {
    struct Holding {
        double quantity;
        VanillaOption option;
    };

    std::vector<Holding> options;
    double stock = 0.0; // units of the underlying

    double value(const MarketState& m, double elapsed = 0.0) const {
        double v = stock * m.spot.value;
        for (const Holding& h : options) {
            VanillaOption aged = h.option;
            aged.maturity.value = std::max(0.0, aged.maturity.value - elapsed);
            v += h.quantity * black_scholes_price(aged, m);
        }
        return v;
    }

    // Raw sensitivities of the whole book (docs/conventions.md units); the stock adds to delta only.
    Greeks greeks(const MarketState& m) const {
        Greeks total{stock, 0.0, 0.0, 0.0, 0.0};
        for (const Holding& h : options) {
            const Greeks g = black_scholes_greeks(h.option, m);
            total.delta += h.quantity * g.delta;
            total.gamma += h.quantity * g.gamma;
            total.vega += h.quantity * g.vega;
            total.theta += h.quantity * g.theta;
            total.rho += h.quantity * g.rho;
        }
        return total;
    }
};

// Short one straddle (call + put struck at K), with the underlying position that makes the book
// delta-neutral today: the position of the report's central demonstration (report 7.1).
inline Portfolio delta_hedged_short_straddle(Strike k, Maturity t, const MarketState& m) {
    Portfolio p{{{-1.0, {k, t, OptionType::Call}}, {-1.0, {k, t, OptionType::Put}}}, 0.0};
    p.stock = -p.greeks(m).delta;
    return p;
}

} // namespace riskengine
