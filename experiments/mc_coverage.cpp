// Report 5.1: are the standard errors honest? 1,000 independent pricings (one seed each) of an
// at-the-money call and an out-of-the-money digital, 10,000 paths each. If the standard error is
// right, the z-scores (estimate - exact) / SE are standard normal and ~95 % of the 95 % confidence
// intervals contain the exact price (binomial standard deviation 0.69 % over 1,000 trials).

#include <cmath>
#include <cstdint>

#include "harness/experiment.hpp"
#include "riskengine/methods/analytic/black_scholes.hpp"
#include "riskengine/methods/analytic/digital.hpp"
#include "riskengine/methods/montecarlo/engine.hpp"
#include "riskengine/models/gbm.hpp"
#include "riskengine/payoffs/digital.hpp"

using namespace riskengine;

int main(int argc, char** argv) {
    return harness::run("mc_coverage", argc, argv, [](harness::Experiment& exp) {
        const MarketState market{Spot{100}, Rate{0.05}, Rate{0.0}, Vol{0.2}};
        constexpr int replications = 1000;
        const MonteCarloConfig cfg{.paths = 10'000, .threads = 4};
        const double call_exact = black_scholes_price(VanillaOption{Strike{100}, Maturity{1.0}, OptionType::Call}, market);
        const double digital_exact = digital_price(VanillaOption{Strike{130}, Maturity{1.0}, OptionType::Call}, market);
        exp.param("market", "S = 100, r = 0.05, q = 0, sigma = 0.2, T = 1");
        exp.param("replications", replications);
        exp.param("paths_per_replication", cfg.paths);
        exp.param("seed", "SeedKey{i, 0} for the call, SeedKey{i, 1} for the digital, i = 0..999");
        exp.param("call_exact", call_exact);
        exp.param("digital_exact", digital_exact);

        const MonteCarlo<GBM, VanillaPayoff> call({100.0, OptionType::Call}, Maturity{1.0}, cfg);
        const MonteCarlo<GBM, DigitalPayoff> digital({130.0, OptionType::Call}, Maturity{1.0}, cfg);
        auto csv = exp.csv({"payoff", "replication", "estimate", "std_error", "z_score", "covered_95"});
        for (int i = 0; i < replications; ++i) {
            const auto seed = static_cast<std::uint64_t>(i);
            const Estimate c = call.price(market, SeedKey{seed, 0});
            const Estimate d = digital.price(market, SeedKey{seed, 1});
            const double zc = (c.value - call_exact) / c.std_error, zd = (d.value - digital_exact) / d.std_error;
            csv.row("atm_call", i, c.value, c.std_error, zc, std::abs(zc) <= 1.959963984540054 ? 1 : 0);
            csv.row("otm_digital_k130", i, d.value, d.std_error, zd, std::abs(zd) <= 1.959963984540054 ? 1 : 0);
        }
    });
}
