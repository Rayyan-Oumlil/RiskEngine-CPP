// Report 5.2: efficiency of variance-reduction techniques, including where they fail.
//
// Efficiency = 1 / (variance per path x CPU time per path) (Glasserman): a technique that halves
// the variance but doubles the cost is worth nothing. "Variance per path" is the estimator's
// variance times the number of simulated paths (an antithetic pair counts as two paths), so every
// technique is charged for the paths it simulates. Times are single-threaded, the minimum of three
// runs, and include the pilot run of a control variate. The statistical columns are deterministic;
// the timing columns depend on the machine recorded in the metadata.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "harness/experiment.hpp"
#include "riskengine/methods/analytic/black_scholes.hpp"
#include "riskengine/methods/montecarlo/engine.hpp"
#include "riskengine/models/gbm.hpp"
#include "riskengine/payoffs/asian.hpp"

using namespace riskengine;

namespace {

const MarketState kMarket{Spot{100}, Rate{0.05}, Rate{0.0}, Vol{0.2}};
constexpr Maturity kMaturity{1.0};
constexpr std::uint64_t kPaths = 1'000'000;
constexpr int kRepeats = 3;

MonteCarloConfig config(bool antithetic, std::uint32_t steps = 1) {
    return MonteCarloConfig{.paths = kPaths, .steps = steps, .threads = 1, .antithetic = antithetic};
}

struct Row {
    Estimate estimate;
    double seconds;
};

template <class Pricer>
Row measure(const Pricer& pricer) {
    Row row{pricer.price(kMarket, SeedKey{2026}), std::numeric_limits<double>::infinity()};
    for (int r = 0; r < kRepeats; ++r) {
        const auto start = std::chrono::steady_clock::now();
        const Estimate e = pricer.price(kMarket, SeedKey{2026});
        const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
        row.seconds = std::min(row.seconds, elapsed.count());
        if (e.value != row.estimate.value) throw std::runtime_error("pricer is not deterministic");
    }
    return row;
}

double call_price(double k) {
    return black_scholes_price(VanillaOption{Strike{k}, kMaturity, OptionType::Call}, kMarket);
}

} // namespace

int main(int argc, char** argv) {
    return harness::run("mc_efficiency", argc, argv, [](harness::Experiment& exp) {
        exp.param("market", "S = 100, r = 0.05, q = 0, sigma = 0.2, T = 1");
        exp.param("paths", kPaths);
        exp.param("threads", 1);
        exp.param("timing", "min of 3 runs, steady_clock, pilot run included");
        exp.param("seed", "SeedKey{2026} for every row");
        exp.param("pilot_paths", MonteCarloConfig{.paths = 0}.pilot_paths);

        auto csv = exp.csv({"case", "technique", "estimate", "std_error", "exact", "variance_per_path",
                            "seconds_per_path", "efficiency", "efficiency_vs_plain"});
        double plain_efficiency = 0.0;
        auto record = [&](const std::string& name, const std::string& technique, const Row& row, double exact) {
            const double var_per_path = row.estimate.std_error * row.estimate.std_error * static_cast<double>(kPaths);
            const double sec_per_path = row.seconds / static_cast<double>(kPaths);
            const double efficiency = 1.0 / (var_per_path * sec_per_path);
            if (technique == "plain") plain_efficiency = efficiency;
            csv.row(name, technique, row.estimate.value, row.estimate.std_error, exact, var_per_path, sec_per_path,
                    efficiency, efficiency / plain_efficiency);
        };
        using Call = MonteCarlo<GBM, VanillaPayoff>;
        using CallSpot = MonteCarlo<GBM, VanillaPayoff, TerminalSpotControl>;

        for (double k : {100.0, 70.0, 160.0}) {
            const std::string name = k == 100.0 ? "atm_call" : (k < 100.0 ? "itm_call_k70" : "otm_call_k160");
            const VanillaPayoff payoff{k, OptionType::Call};
            record(name, "plain", measure(Call(payoff, kMaturity, config(false))), call_price(k));
            record(name, "antithetic", measure(Call(payoff, kMaturity, config(true))), call_price(k));
            record(name, "spot_control", measure(CallSpot(payoff, kMaturity, config(false))), call_price(k));
            record(name, "antithetic+spot_control", measure(CallSpot(payoff, kMaturity, config(true))), call_price(k));
        }

        const double straddle_exact =
            call_price(100.0) + black_scholes_price(VanillaOption{Strike{100}, kMaturity, OptionType::Put}, kMarket);
        using Straddle = MonteCarlo<GBM, StraddlePayoff>;
        record("atm_straddle", "plain", measure(Straddle({100.0}, kMaturity, config(false))), straddle_exact);
        record("atm_straddle", "antithetic", measure(Straddle({100.0}, kMaturity, config(true))), straddle_exact);

        constexpr std::uint32_t fixings = 12;
        const ArithmeticAsianPayoff asian{100.0, OptionType::Call};
        const GeometricAsianControl geometric{100.0, OptionType::Call};
        const double no_closed_form = std::numeric_limits<double>::quiet_NaN();
        using Asian = MonteCarlo<GBM, ArithmeticAsianPayoff>;
        using AsianGeo = MonteCarlo<GBM, ArithmeticAsianPayoff, GeometricAsianControl>;
        record("asian_atm_12_fixings", "plain", measure(Asian(asian, kMaturity, config(false, fixings))), no_closed_form);
        record("asian_atm_12_fixings", "antithetic", measure(Asian(asian, kMaturity, config(true, fixings))), no_closed_form);
        record("asian_atm_12_fixings", "geometric_control",
               measure(AsianGeo(asian, kMaturity, config(false, fixings), geometric)), no_closed_form);
        record("asian_atm_12_fixings", "antithetic+geometric_control",
               measure(AsianGeo(asian, kMaturity, config(true, fixings), geometric)), no_closed_form);
    });
}
