// Quickstart: one European call priced by every method in the library, plus a Monte Carlo Greek
// and the same call under Heston and Merton. Build with the rest of the project and run:
//   build/examples/quickstart
#include <cstdio>

#include "riskengine/methods/analytic/black_scholes.hpp"
#include "riskengine/methods/analytic/implied_vol.hpp"
#include "riskengine/methods/montecarlo/engine.hpp"
#include "riskengine/methods/montecarlo/greeks.hpp"
#include "riskengine/methods/tree/binomial.hpp"
#include "riskengine/models/gbm.hpp"
#include "riskengine/models/heston.hpp"
#include "riskengine/models/merton.hpp"
#include "riskengine/payoffs/vanilla.hpp"

using namespace riskengine;

int main() {
    const MarketState market{Spot{100}, Rate{0.05}, Rate{0.0}, Vol{0.20}};
    const VanillaOption call{Strike{100}, Maturity{1.0}, OptionType::Call};

    // Closed form and Greeks (raw units, see docs/conventions.md).
    const double bs = black_scholes_price(call, market);
    const Greeks g = black_scholes_greeks(call, market);
    std::printf("Black-Scholes        %.6f   delta %.4f  gamma %.4f  vega/pt %.4f\n", bs, g.delta, g.gamma,
                vega_per_vol_point(g));

    // Implied vol round trip.
    const ImpliedVolResult iv = implied_vol(call, bs, market.spot, market.rate, market.div);
    std::printf("Implied vol          %.12f (%s)\n", iv.vol, to_string(iv.status));

    // Trees: CRR oscillates with n, Leisen-Reimer converges at second order.
    std::printf("CRR, n = 1000        %.6f   error %+.2e\n", binomial_price(TreeMethod::Crr, call, Exercise::European, market, 1000),
                binomial_price(TreeMethod::Crr, call, Exercise::European, market, 1000) - bs);
    std::printf("Leisen-Reimer, 1001  %.6f   error %+.2e\n",
                binomial_price(TreeMethod::LeisenReimer, call, Exercise::European, market, 1001),
                binomial_price(TreeMethod::LeisenReimer, call, Exercise::European, market, 1001) - bs);
    const VanillaOption put{Strike{100}, Maturity{1.0}, OptionType::Put};
    std::printf("American put (BBS-R) %.6f   European put %.6f\n",
                binomial_price(TreeMethod::BbsRichardson, put, Exercise::American, market, 2000),
                black_scholes_price(put, market));

    // Monte Carlo: every estimate carries its standard error; the result does not depend on threads.
    const MonteCarlo<GBM, VanillaPayoff> mc(VanillaPayoff{100.0, OptionType::Call}, Maturity{1.0},
                                            MonteCarloConfig{.paths = 1u << 20, .threads = 4, .antithetic = true});
    const Estimate e = mc.price(market, SeedKey{2026});
    std::printf("Monte Carlo (anti.)  %.6f ± %.6f  (%.1f standard errors from BS)\n", e.value, e.std_error,
                (e.value - bs) / e.std_error);

    // A Monte Carlo Greek: pathwise delta.
    const Estimate delta = mc_greek(Greek::Delta, GreekMethod::Pathwise, VanillaPayoff{100.0, OptionType::Call},
                                    Maturity{1.0}, market, SeedKey{7}, GreekConfig{.paths = 1u << 18, .threads = 4});
    std::printf("Pathwise delta       %.4f ± %.4f  (exact %.4f)\n", delta.value, delta.std_error, g.delta);

    // The same call under Heston and Merton (closed forms).
    const HestonParams heston{0.04, 2.0, 0.04, 0.3, -0.7};
    const MertonParams merton{0.25, -0.2, 0.1};
    std::printf("Heston (CF)          %.6f\n", heston_price(OptionType::Call, 100, 100, 0.05, 0.0, 1.0, heston));
    std::printf("Merton (series)      %.6f\n", merton_price(OptionType::Call, 100, 100, 0.05, 0.0, 0.2, 1.0, merton));
    return 0;
}
