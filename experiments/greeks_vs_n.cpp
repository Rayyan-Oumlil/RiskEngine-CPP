// Report 6.1-6.3: convergence rate of each Greek estimator. For N = 2^8 ... 2^18 paths, the RMSE of
// a single estimate over 32 replications; finite differences use, at each N, the best bump of the
// grid 10^(k/4) (an oracle choice: the rate it gives is the best that estimator can achieve). Theory
// for the RMSE slopes: independent seeds -1/3 (delta) and -1/4 (gamma); CRN -1/2 for the call delta,
// -2/5 for the digital delta and the call gamma, -2/7 for the digital gamma; unbiased estimators
// -1/2, except where they converge to the wrong value (pathwise digital delta, mixed digital gamma).

#include <cmath>
#include <cstdint>
#include <limits>

#include "greeks_common.hpp"
#include "harness/experiment.hpp"

using namespace riskengine;
using namespace riskengine::greeks_study;

int main(int argc, char** argv) {
    return harness::run("greeks_vs_n", argc, argv, [](harness::Experiment& exp) {
        const MarketState market{Spot{100}, Rate{0.05}, Rate{0.0}, Vol{0.2}};
        const Maturity t{1.0};
        constexpr std::uint32_t replications = 32;
        exp.param("market", "S = K = 100, r = 0.05, q = 0, sigma = 0.2, T = 1");
        exp.param("replications", replications);
        exp.param("seeds", "SeedKey{2026, r}, r = 0..31");
        exp.param("h_grid", "best of 10^(k/4), k = -16..-2, per N (finite differences only)");

        auto csv = exp.csv({"payoff", "greek", "method", "paths", "best_h", "bias", "sd", "rmse", "exact"});
        const VanillaPayoff call{100, OptionType::Call};
        const DigitalPayoff digital{100, OptionType::Call};
        for (bool is_digital : {false, true})
            for (Greek greek : {Greek::Delta, Greek::Gamma}) {
                const double exact = exact_greek(is_digital, greek, 100, t, market);
                const char* payoff_name = is_digital ? "digital" : "call";
                for (GreekMethod method : {GreekMethod::FdIndependent, GreekMethod::FdCrn, GreekMethod::Pathwise,
                                           GreekMethod::LikelihoodRatio, GreekMethod::Mixed}) {
                    if (!supports(method, greek)) continue;
                    const bool fd = method == GreekMethod::FdIndependent || method == GreekMethod::FdCrn;
                    for (int p = 8; p <= 18; ++p) {
                        const std::uint64_t paths = std::uint64_t{1} << p;
                        ErrorStats best{0, 0, 0, std::numeric_limits<double>::infinity()};
                        double best_h = std::numeric_limits<double>::quiet_NaN();
                        for (int k = fd ? -16 : 0; k <= (fd ? -2 : 0); ++k) {
                            const double h = std::pow(10.0, k / 4.0);
                            const GreekConfig cfg{.paths = paths, .bump = h, .threads = 4};
                            const ErrorStats e =
                                is_digital ? replicate(greek, method, digital, t, market, cfg, exact, replications, 2026)
                                           : replicate(greek, method, call, t, market, cfg, exact, replications, 2026);
                            if (e.rmse < best.rmse) {
                                best = e;
                                if (fd) best_h = h;
                            }
                        }
                        csv.row(payoff_name, to_string(greek), to_string(method), paths, best_h, best.bias, best.sd,
                                best.rmse, exact);
                    }
                }
            }
    });
}
