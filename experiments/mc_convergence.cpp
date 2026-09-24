// Report 5.1: Monte Carlo convergence. Prices the canonical at-the-money call with N = 2^10 ... 2^22
// paths (an independent seed per N) and records the estimate, its standard error and its actual
// error against Black-Scholes. The standard error must fall as N^{-1/2}; the actual error must stay
// within a few standard errors of zero, not beat the rate.

#include <cmath>
#include <cstdint>

#include "harness/experiment.hpp"
#include "riskengine/methods/analytic/black_scholes.hpp"
#include "riskengine/methods/montecarlo/engine.hpp"
#include "riskengine/models/gbm.hpp"

using namespace riskengine;

int main(int argc, char** argv) {
    return harness::run("mc_convergence", argc, argv, [](harness::Experiment& exp) {
        const MarketState market{Spot{100}, Rate{0.05}, Rate{0.0}, Vol{0.2}};
        const VanillaOption call{Strike{100}, Maturity{1.0}, OptionType::Call};
        const double exact = black_scholes_price(call, market);
        constexpr int first_power = 10, last_power = 22;
        exp.param("payoff", "European call, K = 100, T = 1");
        exp.param("market", "S = 100, r = 0.05, q = 0, sigma = 0.2");
        exp.param("exact", exact);
        exp.param("seed", "SeedKey{2026, p} for N = 2^p");
        exp.param("blocks", 64);

        auto csv = exp.csv({"paths", "estimate", "std_error", "abs_error", "z_score"});
        for (int p = first_power; p <= last_power; ++p) {
            const std::uint64_t paths = std::uint64_t{1} << p;
            const MonteCarlo<GBM, VanillaPayoff> pricer({100.0, OptionType::Call}, Maturity{1.0},
                                                       MonteCarloConfig{.paths = paths, .threads = 4});
            const Estimate e = pricer.price(market, SeedKey{2026, static_cast<std::uint32_t>(p)});
            csv.row(paths, e.value, e.std_error, std::abs(e.value - exact), (e.value - exact) / e.std_error);
        }
    });
}
