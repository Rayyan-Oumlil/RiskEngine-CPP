#pragma once

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>

#include "riskengine/concepts.hpp"
#include "riskengine/core/simulation.hpp"
#include "riskengine/greeks/dual.hpp"
#include "riskengine/models/gbm.hpp"

namespace riskengine {

// Monte Carlo estimators of delta and gamma for a terminal payoff under GBM (report 6).
//
// Every estimator is an average of i.i.d. per-path quantities, so each Estimate carries an honest
// standard error, and each is a pure function of (market, key, config): the numbers of block b come
// from RandomStream(key, b), as in the pricing engine.
//
//   FdIndependent    central difference of prices simulated with independent normals for S(1+h),
//                    S(1-h) (and S for gamma): variance O(1/(N h^2)) for delta, O(1/(N h^4)) for gamma.
//   FdCrn            the same bumps on the same normal (common random numbers): variance O(1/N) for a
//                    Lipschitz payoff, but O(1/(N h)) when the payoff jumps (a digital delta) or kinks
//                    under a second difference (a call gamma).
//   Pathwise         d/dS of the discounted payoff along the path: e^{-rT} f'(S_T) S_T / S. Unbiased
//                    when f is Lipschitz; for a digital f' = 0 a.e., so it converges to 0.
//   PathwiseDual     the same derivative computed by forward automatic differentiation (Dual), a
//                    cross-check of the hand-written formula.
//   LikelihoodRatio  differentiates the density instead of the payoff: e^{-rT} f(S_T) times the score
//                    Z / (S sigma sqrt(T)) for delta, (Z^2 - 1 - Z sigma sqrt(T)) / (S^2 sigma^2 T) for
//                    gamma. Unbiased for any payoff, but its variance grows as T -> 0.
//   Mixed            gamma as the likelihood ratio applied to the pathwise delta:
//                    e^{-rT} f'(S_T) S_T (Z / (sigma sqrt(T)) - 1) / S^2 (Glasserman 2003, 7.3).
enum class GreekMethod { FdIndependent, FdCrn, Pathwise, PathwiseDual, LikelihoodRatio, Mixed };

inline const char* to_string(GreekMethod m) {
    switch (m) {
        case GreekMethod::FdIndependent: return "fd_independent";
        case GreekMethod::FdCrn: return "fd_crn";
        case GreekMethod::Pathwise: return "pathwise";
        case GreekMethod::PathwiseDual: return "pathwise_dual";
        case GreekMethod::LikelihoodRatio: return "likelihood_ratio";
        case GreekMethod::Mixed: return "mixed";
    }
    return "unknown";
}

enum class Greek { Delta, Gamma };

struct GreekConfig {
    std::uint64_t paths;
    double bump = 0.01;        // finite differences only: relative bump h, S -> S (1 +- h)
    std::uint32_t blocks = 64; // as in BlockPlan: part of the result's identity
    unsigned threads = 1;      // never changes the result
};

// Pathwise and PathwiseDual estimate delta only (a kinked payoff has no pathwise second derivative);
// Mixed estimates gamma only. mc_greek rejects any other combination with an exception, in every
// build, rather than return a Greek under the wrong name.
inline bool supports(GreekMethod method, Greek greek) {
    if (method == GreekMethod::Pathwise || method == GreekMethod::PathwiseDual) return greek == Greek::Delta;
    if (method == GreekMethod::Mixed) return greek == Greek::Gamma;
    return true;
}

// Terminal payoffs usable by every estimator: scalar-generic (for Dual) with an a.e. derivative.
template <class F>
concept DifferentiableTerminalPayoff = TerminalPayoff<F> && requires(const F& f, double s, Dual d) {
    { f.derivative(s) } -> std::convertible_to<double>;
    { f(d) } -> std::convertible_to<Dual>;
};

template <DifferentiableTerminalPayoff Payoff>
Estimate mc_greek(Greek greek, GreekMethod method, const Payoff& payoff, Maturity maturity, const MarketState& m,
                  SeedKey key, const GreekConfig& config) {
    if (!supports(method, greek))
        throw std::invalid_argument(std::string(to_string(method)) + " does not estimate this Greek");
    const double s = m.spot.value, r = m.rate.value, q = m.div.value, sigma = m.vol.value, t = maturity.value;
    const double discount = std::exp(-r * t);
    const double sst = sigma * std::sqrt(t);
    const auto terminal = [&](double spot, double z) { return gbm_terminal_spot(spot, r, q, sigma, t, z); };
    // Realized bumped spots: divide by the step actually represented in floating point.
    const double up = s * (1.0 + config.bump), down = s * (1.0 - config.bump);
    const double half_step = 0.5 * (up - down);

    const Welford w = reduce_blocks(BlockPlan{config.paths, config.blocks, config.threads},
                                    [&](std::uint32_t b, std::uint64_t n) {
        RandomStream rng(key, b);
        Welford acc;
        for (std::uint64_t i = 0; i < n; ++i) {
            double sample = 0.0;
            switch (method) {
                case GreekMethod::FdIndependent: {
                    const double f_up = payoff(terminal(up, rng.normal()));
                    const double f_down = payoff(terminal(down, rng.normal()));
                    if (greek == Greek::Delta) {
                        sample = (f_up - f_down) / (up - down);
                    } else {
                        const double f_mid = payoff(terminal(s, rng.normal()));
                        sample = (f_up - 2.0 * f_mid + f_down) / (half_step * half_step);
                    }
                    break;
                }
                case GreekMethod::FdCrn: {
                    const double z = rng.normal();
                    const double f_up = payoff(terminal(up, z)), f_down = payoff(terminal(down, z));
                    sample = greek == Greek::Delta
                                 ? (f_up - f_down) / (up - down)
                                 : (f_up - 2.0 * payoff(terminal(s, z)) + f_down) / (half_step * half_step);
                    break;
                }
                case GreekMethod::Pathwise: {
                    const double st = terminal(s, rng.normal());
                    sample = payoff.derivative(st) * st / s;
                    break;
                }
                case GreekMethod::PathwiseDual: {
                    const Dual st = gbm_terminal_spot(Dual::variable(s), Dual(r), Dual(q), Dual(sigma), t, rng.normal());
                    sample = static_cast<Dual>(payoff(st)).d;
                    break;
                }
                case GreekMethod::LikelihoodRatio: {
                    const double z = rng.normal();
                    const double f = payoff(terminal(s, z));
                    sample = greek == Greek::Delta ? f * z / (s * sst)
                                                   : f * (z * z - 1.0 - z * sst) / (s * s * sst * sst);
                    break;
                }
                case GreekMethod::Mixed: {
                    const double z = rng.normal();
                    const double st = terminal(s, z);
                    sample = payoff.derivative(st) * st * (z / sst - 1.0) / (s * s);
                    break;
                }
            }
            acc.add(discount * sample);
        }
        return acc;
    });
    return to_estimate(w);
}

} // namespace riskengine
