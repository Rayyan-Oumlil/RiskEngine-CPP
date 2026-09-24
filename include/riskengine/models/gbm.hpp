#pragma once

#include <cmath>
#include <cstddef>
#include <span>

#include "riskengine/core/market.hpp"

namespace riskengine {

// S_T = S e^{(r - q - sigma^2/2) T + sigma sqrt(T) z}, generic in the scalar type of the market
// inputs so that Dual inputs give exact pathwise derivatives (delta, vega, rho) of S_T.
template <class T>
T gbm_terminal_spot(const T& spot, const T& rate, const T& div, const T& vol, double t, double z) {
    using std::exp;
    return spot * exp((rate - div - 0.5 * vol * vol) * t + vol * std::sqrt(t) * z);
}

// Geometric Brownian motion under the risk-neutral measure, dS = (r - q) S dt + sigma S dW.
// The log-space step is exact for any dt, so GBM paths carry no discretization bias.
class GBM {
public:
    using State = double; // spot
    static constexpr std::size_t factors = 1;

    explicit GBM(const MarketState& m)
        : spot_(m.spot.value),
          drift_(m.rate.value - m.div.value - 0.5 * m.vol.value * m.vol.value),
          vol_(m.vol.value) {}

    State initial_state() const { return spot_; }

    State step(State s, double dt, std::span<const double> z) const {
        return s * std::exp(drift_ * dt + vol_ * std::sqrt(dt) * z[0]);
    }

    double spot(State s) const { return s; }

private:
    double spot_;
    double drift_; // r - q - sigma^2 / 2
    double vol_;
};

} // namespace riskengine
