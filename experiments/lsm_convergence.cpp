// Report: Longstaff-Schwartz convergence. Prices the canonical American put with N = 2^12 ... 2^18
// paths (an independent seed per N, 50 exercise dates) and records the estimate, its standard error,
// its regression bias estimate, and its actual error against the tree reference (3.2, n = 20,000).
// The standard error must fall as N^{-1/2}, same as the European Monte Carlo engine; the bias
// estimate should not grow with N, since it comes from the fixed {1, S, S^2} basis, not from path
// count.

#include <cmath>
#include <cstdint>

#include "harness/experiment.hpp"
#include "riskengine/methods/montecarlo/longstaff_schwartz.hpp"
#include "riskengine/methods/tree/binomial.hpp"
#include "riskengine/models/gbm.hpp"

using namespace riskengine;

int main(int argc, char** argv) {
    return harness::run("lsm_convergence", argc, argv, [](harness::Experiment& exp) {
        const MarketState market{Spot{100}, Rate{0.05}, Rate{0.0}, Vol{0.2}};
        const VanillaOption put{Strike{100}, Maturity{1.0}, OptionType::Put};
        const double tree_reference = binomial_price(TreeMethod::Crr, put, Exercise::American, market, 20000);
        constexpr int first_power = 12, last_power = 18;
        constexpr std::uint32_t steps = 50;
        exp.param("payoff", "American put, K = 100, T = 1");
        exp.param("market", "S = 100, r = 0.05, q = 0, sigma = 0.2");
        exp.param("tree_reference", tree_reference);
        exp.param("tree_reference_source", "CRR, n = 20,000 (docs/model_risk_report.md, Appendix C)");
        exp.param("exercise_dates", steps);
        exp.param("basis", "{1, S, S^2} on in-the-money paths");
        exp.param("seed", "SeedKey{2026, p} for N = 2^p");

        auto csv = exp.csv({"paths", "estimate", "std_error", "bias_estimate", "abs_error", "z_score"});
        for (int p = first_power; p <= last_power; ++p) {
            const std::uint64_t paths = std::uint64_t{1} << p;
            const LongstaffSchwartz<GBM> lsm(put, LongstaffSchwartzConfig{.paths = paths, .steps = steps});
            const Estimate e = lsm.price(market, SeedKey{2026, static_cast<std::uint32_t>(p)});
            csv.row(paths, e.value, e.std_error, e.discretization, std::abs(e.value - tree_reference),
                    (e.value - tree_reference) / e.std_error);
        }
    });
}
