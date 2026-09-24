#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>
#include <span>
#include <vector>

#include "riskengine/concepts.hpp"
#include "riskengine/core/simulation.hpp"

namespace riskengine {

struct MonteCarloConfig {
    std::uint64_t paths;         // simulated paths (antithetic: pairs = paths / 2)
    std::uint32_t steps = 1;     // time steps, equally spaced; also the fixing dates of a path payoff
    std::uint32_t blocks = 64;   // fixed work split: part of the experiment's identity (see BlockPlan)
    unsigned threads = 1;        // never changes the result
    bool antithetic = false;     // pair each path with its reflection z -> -z
};

// Generic Monte Carlo pricer: any PathModel, any terminal or path payoff.
//
// Path i of block b draws its normals from RandomStream(key, b), step-major then factor-major.
// The payoff is evaluated on the terminal spot (TerminalPayoff) or on the spot at the end of each
// step (PathPayoff), then discounted at the risk-free rate. With antithetic variates, one sample is
// the average of a path and its reflection, so the standard error accounts for their correlation.
// The price is a pure function of (market, key, config.paths, config.steps, config.blocks).
template <PathModel Model, class Payoff>
    requires TerminalPayoff<Payoff> || PathPayoff<Payoff>
class MonteCarlo {
public:
    MonteCarlo(Payoff payoff, Maturity maturity, MonteCarloConfig config)
        : payoff_(payoff), maturity_(maturity), config_(config) {
        assert(config_.steps >= 1);
        assert(!config_.antithetic || config_.paths % 2 == 0);
    }

    Estimate price(const MarketState& m, SeedKey key) const {
        const Model model(m);
        const double discount = std::exp(-m.rate.value * maturity_.value);
        const double dt = maturity_.value / config_.steps;
        const std::size_t draws = static_cast<std::size_t>(config_.steps) * Model::factors;
        const std::uint64_t samples = config_.antithetic ? config_.paths / 2 : config_.paths;

        const Welford w = reduce_blocks(
            BlockPlan{samples, config_.blocks, config_.threads}, [&](std::uint32_t b, std::uint64_t n) {
                RandomStream rng(key, b);
                std::vector<double> z(draws), z_reflected(draws), fixings(config_.steps);
                Welford acc;
                for (std::uint64_t i = 0; i < n; ++i) {
                    for (double& x : z) x = rng.normal();
                    double y = payoff_on_path(model, dt, z, fixings);
                    if (config_.antithetic) {
                        for (std::size_t j = 0; j < draws; ++j) z_reflected[j] = -z[j];
                        y = 0.5 * (y + payoff_on_path(model, dt, z_reflected, fixings));
                    }
                    acc.add(discount * y);
                }
                return acc;
            });
        return to_estimate(w);
    }

    const MonteCarloConfig& config() const { return config_; }

private:
    double payoff_on_path(const Model& model, double dt, std::span<const double> z,
                          std::vector<double>& fixings) const {
        typename Model::State s = model.initial_state();
        for (std::uint32_t step = 0; step < config_.steps; ++step) {
            s = model.step(s, dt, z.subspan(step * Model::factors, Model::factors));
            fixings[step] = model.spot(s);
        }
        if constexpr (TerminalPayoff<Payoff>) {
            return payoff_(fixings.back());
        } else {
            return payoff_(std::span<const double>(fixings));
        }
    }

    Payoff payoff_;
    Maturity maturity_;
    MonteCarloConfig config_;
};

} // namespace riskengine
