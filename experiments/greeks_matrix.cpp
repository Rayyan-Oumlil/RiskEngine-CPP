// Report 6.4: recommendation matrix. Payoff (call, digital) x regime (ATM, deep OTM, one-week) x
// Greek (delta, gamma) x estimator, at N = 2^16 paths per estimate. Finite differences use a fixed
// relative bump of 1 %, a common desk convention, rather than the oracle bump of greeks_vs_n.
// Error is the RMSE of a single estimate over 32 replications; efficiency is 1 / (MSE x time), per
// estimate, so an estimator twice as accurate at the same cost is 4 times as efficient. Times are
// single-threaded, the minimum of three runs; the statistical columns are deterministic.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>

#include "greeks_common.hpp"
#include "harness/experiment.hpp"

using namespace riskengine;
using namespace riskengine::greeks_study;

namespace {

template <class Payoff>
double seconds_per_estimate(Greek greek, GreekMethod method, const Payoff& payoff, Maturity t, const MarketState& m,
                            GreekConfig cfg) {
    cfg.threads = 1;
    double best = std::numeric_limits<double>::infinity();
    for (int r = 0; r < 3; ++r) {
        const auto start = std::chrono::steady_clock::now();
        (void)mc_greek(greek, method, payoff, t, m, SeedKey{1}, cfg);
        const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
        best = std::min(best, elapsed.count());
    }
    return best;
}

} // namespace

int main(int argc, char** argv) {
    return harness::run("greeks_matrix", argc, argv, [](harness::Experiment& exp) {
        const MarketState market{Spot{100}, Rate{0.05}, Rate{0.0}, Vol{0.2}};
        constexpr std::uint64_t paths = 1 << 16;
        constexpr std::uint32_t replications = 32;
        constexpr double bump = 0.01;
        exp.param("market", "S = 100, r = 0.05, q = 0, sigma = 0.2");
        exp.param("regimes", "atm: K = 100, T = 1; otm: K = 150, T = 1; short: K = 100, T = 1/52");
        exp.param("paths", paths);
        exp.param("replications", replications);
        exp.param("fd_bump", bump);
        exp.param("seeds", "SeedKey{2026, r}, r = 0..31; timing SeedKey{1}, single thread, min of 3");

        auto csv = exp.csv({"payoff", "regime", "greek", "method", "exact", "bias", "sd", "rmse", "relative_rmse",
                            "seconds_per_estimate", "efficiency"});
        struct Regime {
            const char* name;
            double strike;
            double maturity;
        };
        for (const Regime& regime : {Regime{"atm", 100, 1.0}, Regime{"otm", 150, 1.0}, Regime{"short", 100, 1.0 / 52}})
            for (bool is_digital : {false, true})
                for (Greek greek : {Greek::Delta, Greek::Gamma}) {
                    const Maturity t{regime.maturity};
                    const VanillaPayoff call{regime.strike, OptionType::Call};
                    const DigitalPayoff digital{regime.strike, OptionType::Call};
                    const double exact = exact_greek(is_digital, greek, regime.strike, t, market);
                    for (GreekMethod method : {GreekMethod::FdIndependent, GreekMethod::FdCrn, GreekMethod::Pathwise,
                                               GreekMethod::LikelihoodRatio, GreekMethod::Mixed}) {
                        if (!supports(method, greek)) continue;
                        const GreekConfig cfg{.paths = paths, .bump = bump, .threads = 4};
                        const ErrorStats e =
                            is_digital ? replicate(greek, method, digital, t, market, cfg, exact, replications, 2026)
                                       : replicate(greek, method, call, t, market, cfg, exact, replications, 2026);
                        const double seconds = is_digital ? seconds_per_estimate(greek, method, digital, t, market, cfg)
                                                          : seconds_per_estimate(greek, method, call, t, market, cfg);
                        csv.row(is_digital ? "digital" : "call", regime.name, to_string(greek), to_string(method), exact,
                                e.bias, e.sd, e.rmse, e.rmse / std::abs(exact), seconds, 1.0 / (e.rmse * e.rmse * seconds));
                    }
                }
    });
}
