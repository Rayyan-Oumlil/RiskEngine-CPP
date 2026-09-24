// Report 5.3: randomized quasi-Monte Carlo against pseudo-random Monte Carlo.
//
// For each problem and N = 2^6 ... 2^16, the error of a single N-point estimate is measured as the
// standard deviation over 64 independent replications (independent seeds for pseudo-random,
// independent Owen scramblings for RQMC), which is exactly the spread a user of one estimate sees.
// Pseudo-random Monte Carlo gives N^{-1/2}; scrambled Sobol does better the smoother and the
// lower-dimensional the problem: a call (1-D, a kink), a digital (1-D, a jump), a 12-fixing
// arithmetic Asian (12-D, a kink) and a digital on that average (12-D, a jump that is not aligned
// with any coordinate), the path payoffs with and without a Brownian bridge.

#include <cmath>
#include <cstdint>
#include <string>
#include <utility>

#include "harness/experiment.hpp"
#include "riskengine/methods/montecarlo/engine.hpp"
#include "riskengine/models/gbm.hpp"
#include "riskengine/payoffs/asian.hpp"
#include "riskengine/payoffs/digital.hpp"

using namespace riskengine;

namespace {

const MarketState kMarket{Spot{100}, Rate{0.05}, Rate{0.0}, Vol{0.2}};
constexpr Maturity kMaturity{1.0};
constexpr std::uint32_t kReplications = 64;

MonteCarloConfig config(std::uint64_t paths, std::uint32_t steps, Sampling sampling, bool bridge) {
    return MonteCarloConfig{.paths = paths, .steps = steps, .threads = 4, .sampling = sampling,
                            .replications = kReplications, .brownian_bridge = bridge};
}

// Mean and single-estimate standard deviation over kReplications replications.
template <class Payoff>
std::pair<double, double> measure(const Payoff& payoff, std::uint64_t paths, std::uint32_t steps, Sampling sampling,
                                  bool bridge) {
    const MonteCarlo<GBM, Payoff> pricer(payoff, kMaturity, config(paths, steps, sampling, bridge));
    if (sampling == Sampling::RandomizedQmc) {
        const Estimate e = pricer.price(kMarket, SeedKey{2026});
        return {e.value, e.std_error * std::sqrt(static_cast<double>(kReplications))};
    }
    Welford replicates;
    for (std::uint32_t r = 0; r < kReplications; ++r) replicates.add(pricer.price(kMarket, SeedKey{2026, r}).value);
    return {replicates.mean(), std::sqrt(replicates.variance())};
}

} // namespace

int main(int argc, char** argv) {
    return harness::run("qmc_convergence", argc, argv, [](harness::Experiment& exp) {
        exp.param("market", "S = 100, r = 0.05, q = 0, sigma = 0.2, T = 1");
        exp.param("replications", kReplications);
        exp.param("seeds", "pseudo-random: SeedKey{2026, r}; RQMC: SeedKey{2026}, scramblings r = 0..63");
        exp.param("problems", "ATM call; digital call K = 130; arithmetic Asian ATM, 12 fixings; "
                              "digital on the 12-fixing arithmetic average, K = 100");

        auto csv = exp.csv({"problem", "method", "paths", "mean", "sd_single_estimate"});
        auto sweep = [&](const std::string& problem, const std::string& method, auto payoff, std::uint32_t steps,
                         Sampling sampling, bool bridge) {
            for (int p = 6; p <= 16; ++p) {
                const std::uint64_t n = std::uint64_t{1} << p;
                const auto [mean, sd] = measure(payoff, n, steps, sampling, bridge);
                csv.row(problem, method, n, mean, sd);
            }
        };
        const VanillaPayoff call{100.0, OptionType::Call};
        const DigitalPayoff digital{130.0, OptionType::Call};
        const ArithmeticAsianPayoff asian{100.0, OptionType::Call};
        const ArithmeticAsianDigitalPayoff asian_digital{100.0, OptionType::Call};
        for (auto [method, sampling] : {std::pair{"pseudo_random", Sampling::PseudoRandom},
                                        std::pair{"rqmc", Sampling::RandomizedQmc}}) {
            sweep("atm_call", method, call, 1, sampling, false);
            sweep("digital_k130", method, digital, 1, sampling, false);
            sweep("asian_12", method, asian, 12, sampling, false);
            sweep("asian_digital_12", method, asian_digital, 12, sampling, false);
        }
        sweep("asian_12", "rqmc_bridge", asian, 12, Sampling::RandomizedQmc, true);
        sweep("asian_digital_12", "rqmc_bridge", asian_digital, 12, Sampling::RandomizedQmc, true);
    });
}
