// Report 6.1-6.2: error of finite-difference Greeks against the relative bump h, with independent
// seeds and with common random numbers, for an ATM call and an ATM digital (delta and gamma), at a
// fixed N = 2^16 paths. Bias grows as h^2 to the right; with independent seeds the variance grows as
// 1/h^2 (delta) or 1/h^4 (gamma) to the left; with CRN it does not grow at all for the call delta,
// but still grows as 1/h for the digital delta and the call gamma. The unbiased estimators
// (pathwise, likelihood ratio, mixed) are recorded at the same N, as h-independent reference levels.

#include <cmath>
#include <cstdint>
#include <limits>

#include "greeks_common.hpp"
#include "harness/experiment.hpp"

using namespace riskengine;
using namespace riskengine::greeks_study;

int main(int argc, char** argv) {
    return harness::run("greeks_vs_h", argc, argv, [](harness::Experiment& exp) {
        const MarketState market{Spot{100}, Rate{0.05}, Rate{0.0}, Vol{0.2}};
        const Maturity t{1.0};
        constexpr std::uint64_t paths = 1 << 16;
        constexpr std::uint32_t replications = 32;
        exp.param("market", "S = K = 100, r = 0.05, q = 0, sigma = 0.2, T = 1");
        exp.param("paths", paths);
        exp.param("replications", replications);
        exp.param("seeds", "SeedKey{2026, r}, r = 0..31, the same for every row");
        exp.param("h_grid", "10^(k/4), k = -16..-2");

        auto csv = exp.csv({"payoff", "greek", "method", "h", "mean", "sd", "bias", "rmse", "exact"});
        const VanillaPayoff call{100, OptionType::Call};
        const DigitalPayoff digital{100, OptionType::Call};
        for (bool is_digital : {false, true})
            for (Greek greek : {Greek::Delta, Greek::Gamma}) {
                const double exact = exact_greek(is_digital, greek, 100, t, market);
                const char* payoff_name = is_digital ? "digital" : "call";
                auto row = [&](GreekMethod method, double h) {
                    const GreekConfig cfg{.paths = paths, .bump = std::isnan(h) ? 0.01 : h, .threads = 4};
                    const ErrorStats e = is_digital
                                             ? replicate(greek, method, digital, t, market, cfg, exact, replications, 2026)
                                             : replicate(greek, method, call, t, market, cfg, exact, replications, 2026);
                    csv.row(payoff_name, to_string(greek), to_string(method), h, e.mean, e.sd, e.bias, e.rmse, exact);
                };
                for (int k = -16; k <= -2; ++k)
                    for (GreekMethod m : {GreekMethod::FdIndependent, GreekMethod::FdCrn}) row(m, std::pow(10.0, k / 4.0));
                const double no_bump = std::numeric_limits<double>::quiet_NaN();
                for (GreekMethod m : {GreekMethod::Pathwise, GreekMethod::LikelihoodRatio, GreekMethod::Mixed})
                    if (supports(m, greek)) row(m, no_bump);
            }
    });
}
