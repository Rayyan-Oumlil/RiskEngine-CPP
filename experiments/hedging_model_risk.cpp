// Report 8.2: hedging in the wrong world. A one-year ATM call (K = 100) is sold at its Black-Scholes
// price at 20 % vol and delta-hedged with the Black-Scholes delta at 20 %, rebalanced 16, 63, 252 or
// 1,008 times a year (each four times the last, and each dividing the simulation grid), while the world follows Black-Scholes, Heston or Merton. The three are
// calibrated to the same one-year ATM price (as in model_risk_exotics), so the premium is fair in
// every world and any systematic loss comes from the hedge, not from the price.
//
// Each world simulates 20,000 paths on a grid of 1,008 steps (Heston by QE); the four hedging
// frequencies rebalance on sub-grids of the same paths, so their P&L differ only by the hedge. The
// P&L is the discounted value at maturity of premium + hedge - payoff, per option.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "harness/experiment.hpp"
#include "riskengine/core/rng/random_stream.hpp"
#include "riskengine/methods/analytic/black_scholes.hpp"
#include "riskengine/models/gbm.hpp"
#include "riskengine/models/heston.hpp"
#include "riskengine/models/merton.hpp"
#include "riskengine/risk/var.hpp"

using namespace riskengine;

namespace {

constexpr double kSpot = 100, kStrike = 100, kRate = 0.05, kT = 1.0, kVol = 0.2;
constexpr std::uint32_t kGrid = 1008, kPaths = 20'000;
constexpr std::uint32_t kFrequencies[] = {16, 63, 252, 1008};
// Rebalancing dates must fall on the grid, or the last interval would be shorter than the others.
static_assert(kGrid % 16 == 0 && kGrid % 63 == 0 && kGrid % 252 == 0 && kGrid % 1008 == 0);

double bs_delta(double s, double tau) {
    const VanillaOption o{Strike{kStrike}, Maturity{tau}, OptionType::Call};
    return black_scholes_greeks(o, MarketState{Spot{s}, Rate{kRate}, Rate{0.0}, Vol{kVol}}).delta;
}

// Hedging P&L of every path, for every frequency, in one world. pnl[f][path].
template <class Model>
std::vector<std::vector<double>> hedge(const Model& model, std::uint32_t world) {
    const double dt = kT / kGrid, premium = black_scholes_price(VanillaOption{Strike{kStrike}, Maturity{kT}, OptionType::Call},
                                                                MarketState{Spot{kSpot}, Rate{kRate}, Rate{0.0}, Vol{kVol}});
    std::vector<std::vector<double>> pnl(std::size(kFrequencies), std::vector<double>(kPaths));
    std::vector<double> spots(kGrid + 1), z(Model::factors);
    for (std::uint32_t p = 0; p < kPaths; ++p) {
        RandomStream rng(SeedKey{2026, world}, p);
        auto state = model.initial_state();
        spots[0] = model.spot(state);
        for (std::uint32_t i = 1; i <= kGrid; ++i) {
            for (double& v : z) v = rng.normal();
            state = model.step(state, dt, z);
            spots[i] = model.spot(state);
        }
        for (std::size_t f = 0; f < std::size(kFrequencies); ++f) {
            const std::uint32_t every = kGrid / kFrequencies[f];
            double delta = bs_delta(spots[0], kT);
            double cash = premium - delta * spots[0];
            std::uint32_t last = 0;
            for (std::uint32_t i = every; i < kGrid; i += every) {
                cash *= std::exp(kRate * dt * (i - last)); // interest on the steps actually elapsed
                const double next = bs_delta(spots[i], kT - i * dt);
                cash -= (next - delta) * spots[i];
                delta = next;
                last = i;
            }
            cash *= std::exp(kRate * dt * (kGrid - last));
            const double value = cash + delta * spots[kGrid] - std::max(spots[kGrid] - kStrike, 0.0);
            pnl[f][p] = std::exp(-kRate * kT) * value;
        }
    }
    return pnl;
}

} // namespace

int main(int argc, char** argv) {
    return harness::run("hedging_model_risk", argc, argv, [](harness::Experiment& exp) {
        // The calibration of model_risk_exotics: the same one-year ATM price in all three worlds.
        const HestonParams heston{0.046610715587438928, 1.5, 0.046610715587438928, 0.6, -0.7};
        const MertonParams merton{0.25, -0.2, 0.1};
        const double merton_vol = 0.17129687129725085;
        exp.param("option", "short 1 ATM call, K = 100, T = 1, sold at its Black-Scholes price at 20 % vol");
        exp.param("hedge", "Black-Scholes delta at 20 % vol, self-financing, cash at r = 5 %");
        exp.param("worlds", "Black-Scholes 20 %; Heston kappa 1.5, xi 0.6, rho -0.7, v0 = theta = 0.0466107; "
                            "Merton lambda 0.25, jumps N(-0.2, 0.1^2), diffusion 17.13 % (model_risk_exotics)");
        exp.param("paths", kPaths);
        exp.param("grid_steps", kGrid);
        exp.param("seed", "RandomStream(SeedKey{2026, world}, path), world = 0, 1, 2");

        auto csv = exp.csv({"world", "rebalances_per_year", "mean_pnl", "sd_pnl", "mean_std_error", "q01", "q99",
                            "es99_loss"});
        auto report = [&](const char* world, const std::vector<std::vector<double>>& pnl) {
            for (std::size_t f = 0; f < std::size(kFrequencies); ++f) {
                std::vector<double> x = pnl[f];
                double mean = 0.0;
                for (double v : x) mean += v / x.size();
                double var = 0.0;
                for (double v : x) var += (v - mean) * (v - mean) / (x.size() - 1);
                std::vector<double> losses(x.size());
                std::transform(x.begin(), x.end(), losses.begin(), [](double v) { return -v; });
                const RiskMeasures loss = empirical_var_es(losses, 0.99);
                std::sort(x.begin(), x.end());
                csv.row(world, kFrequencies[f], mean, std::sqrt(var), std::sqrt(var / x.size()),
                        x[x.size() / 100], x[x.size() - 1 - x.size() / 100], loss.es);
            }
        };
        report("black_scholes", hedge(GBM(MarketState{Spot{kSpot}, Rate{kRate}, Rate{0.0}, Vol{kVol}}), 0));
        report("heston", hedge(HestonQE(MarketState{Spot{kSpot}, Rate{kRate}, Rate{0.0}, Vol{0.0}}, heston), 1));
        report("merton", hedge(Merton(MarketState{Spot{kSpot}, Rate{kRate}, Rate{0.0}, Vol{merton_vol}}, merton), 2));
    });
}
