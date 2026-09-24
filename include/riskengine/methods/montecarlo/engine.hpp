#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

#include "riskengine/concepts.hpp"
#include "riskengine/core/simulation.hpp"
#include "riskengine/core/stats/covariance.hpp"
#include "riskengine/methods/montecarlo/control_variates.hpp"

namespace riskengine {

struct MonteCarloConfig {
    std::uint64_t paths;               // simulated paths (antithetic: pairs = paths / 2)
    std::uint32_t steps = 1;           // time steps, equally spaced; also the fixing dates of a path payoff
    std::uint32_t blocks = 64;         // fixed work split: part of the experiment's identity (see BlockPlan)
    unsigned threads = 1;              // never changes the result
    bool antithetic = false;           // pair each path with its reflection z -> -z
    std::uint64_t pilot_paths = 10'000; // control variates only: paths used to estimate beta
};

// Streams with the high bit set are reserved for pilot runs, so a pilot can never reuse the
// numbers of a caller's own stream.
inline constexpr std::uint32_t kPilotStreamBit = 0x8000'0000u;

// Generic Monte Carlo pricer: any PathModel, any terminal or path payoff, optionally with
// antithetic variates and a control variate.
//
// Path i of block b draws its normals from RandomStream(key, b), step-major then factor-major.
// The payoff is evaluated on the terminal spot (TerminalPayoff) or on the spot at the end of each
// step (PathPayoff), then discounted at the risk-free rate. With antithetic variates one sample is
// the average of a path and its reflection, so the standard error accounts for their correlation.
// With a control variate X of known mean mu, the samples are Y - beta (X - mu), where beta is
// estimated on an independent pilot run (estimating it on the same paths would bias the price).
// Because the pilot is independent, the estimate is unbiased for any beta, and the reported
// standard error is the one conditional on the pilot's beta; the sampling error of beta itself
// only affects the variance at second order.
// The price is a pure function of (market, key, paths, steps, blocks, pilot_paths).
template <PathModel Model, class Payoff, class Control = NoControl>
    requires(TerminalPayoff<Payoff> || PathPayoff<Payoff>) &&
            (std::same_as<Control, NoControl> || ControlVariate<Control>)
class MonteCarlo {
    static constexpr bool kHasControl = !std::is_same_v<Control, NoControl>;

public:
    MonteCarlo(Payoff payoff, Maturity maturity, MonteCarloConfig config, Control control = {})
        : payoff_(payoff), maturity_(maturity), config_(config), control_(control) {
        assert(config_.steps >= 1);
        assert(!config_.antithetic || (config_.paths % 2 == 0 && config_.pilot_paths % 2 == 0));
    }

    Estimate price(const MarketState& m, SeedKey key) const {
        assert((key.stream & kPilotStreamBit) == 0);
        const Model model(m);
        const double discount = std::exp(-m.rate.value * maturity_.value);

        if constexpr (!kHasControl) {
            return to_estimate(run<Welford>(model, key, config_.paths, [&](Welford& acc, double y, double) {
                acc.add(discount * y);
            }));
        } else {
            const Covariance pilot =
                run<Covariance>(model, SeedKey{key.seed, key.stream | kPilotStreamBit}, config_.pilot_paths,
                                [](Covariance& acc, double y, double x) { acc.add(x, y); });
            const double var_x = pilot.variance_x();
            const double beta = var_x > 0.0 ? pilot.covariance() / var_x : 0.0;
            const double mu = control_.expectation(m, maturity_, config_.steps);
            return to_estimate(run<Welford>(model, key, config_.paths, [&](Welford& acc, double y, double x) {
                acc.add(discount * y - beta * (discount * x - mu));
            }));
        }
    }

    const MonteCarloConfig& config() const { return config_; }

private:
    // Simulates `paths` paths (pairs when antithetic) and folds each sample's undiscounted payoff y
    // and control value x (0 without a control) into an accumulator of type Acc.
    template <class Acc, class Fold>
    Acc run(const Model& model, SeedKey key, std::uint64_t paths, Fold fold) const {
        const double dt = maturity_.value / config_.steps;
        const std::size_t draws = static_cast<std::size_t>(config_.steps) * Model::factors;
        const std::uint64_t samples = config_.antithetic ? paths / 2 : paths;

        return reduce_blocks(BlockPlan{samples, config_.blocks, config_.threads}, [&](std::uint32_t b, std::uint64_t n) {
            RandomStream rng(key, b);
            std::vector<double> z(draws), fixings(config_.steps);
            Acc acc;
            for (std::uint64_t i = 0; i < n; ++i) {
                for (double& v : z) v = rng.normal();
                auto [y, x] = evaluate(model, dt, z, fixings);
                if (config_.antithetic) {
                    for (double& v : z) v = -v;
                    const auto [y2, x2] = evaluate(model, dt, z, fixings);
                    y = 0.5 * (y + y2);
                    x = 0.5 * (x + x2);
                }
                fold(acc, y, x);
            }
            return acc;
        });
    }

    // Payoff and control value on one path.
    std::pair<double, double> evaluate(const Model& model, double dt, std::span<const double> z,
                                       std::vector<double>& fixings) const {
        typename Model::State s = model.initial_state();
        for (std::uint32_t step = 0; step < config_.steps; ++step) {
            s = model.step(s, dt, z.subspan(step * Model::factors, Model::factors));
            fixings[step] = model.spot(s);
        }
        const std::span<const double> path(fixings);
        double y;
        if constexpr (TerminalPayoff<Payoff>) {
            y = payoff_(path.back());
        } else {
            y = payoff_(path);
        }
        double x = 0.0;
        if constexpr (kHasControl) x = control_(path);
        return {y, x};
    }

    Payoff payoff_;
    Maturity maturity_;
    MonteCarloConfig config_;
    Control control_;
};

} // namespace riskengine
