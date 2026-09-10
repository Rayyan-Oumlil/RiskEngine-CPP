#pragma once

// Shared by the Greeks experiments (report 6): replicate one estimator and summarize its error.

#include <cmath>
#include <cstdint>
#include <string>

#include "riskengine/methods/analytic/digital.hpp"
#include "riskengine/methods/montecarlo/greeks.hpp"
#include "riskengine/payoffs/digital.hpp"
#include "riskengine/payoffs/vanilla.hpp"

namespace riskengine::greeks_study {

// Error of a single estimate, measured over independent replications (SeedKey{seed, r}).
struct ErrorStats {
    double mean;
    double sd;   // spread of one estimate
    double bias; // mean - exact
    double rmse; // sqrt(bias^2 + sd^2), the error a user of one estimate faces
};

template <class Payoff>
ErrorStats replicate(Greek greek, GreekMethod method, const Payoff& payoff, Maturity t, const MarketState& m,
                     const GreekConfig& config, double exact, std::uint32_t replications, std::uint64_t seed) {
    Welford w;
    for (std::uint32_t r = 0; r < replications; ++r)
        w.add(mc_greek(greek, method, payoff, t, m, SeedKey{seed, r}, config).value);
    const double sd = std::sqrt(w.variance()), bias = w.mean() - exact;
    return {w.mean(), sd, bias, std::sqrt(bias * bias + sd * sd)};
}

// Closed-form delta or gamma of a call or a digital call.
inline double exact_greek(bool digital, Greek greek, double strike, Maturity t, const MarketState& m) {
    const VanillaOption o{Strike{strike}, t, OptionType::Call};
    if (digital) return greek == Greek::Delta ? digital_delta(o, m) : digital_gamma(o, m);
    const Greeks g = black_scholes_greeks(o, m);
    return greek == Greek::Delta ? g.delta : g.gamma;
}

inline const char* to_string(Greek g) { return g == Greek::Delta ? "delta" : "gamma"; }

} // namespace riskengine::greeks_study
