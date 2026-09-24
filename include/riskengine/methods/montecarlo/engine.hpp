#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

#include "riskengine/concepts.hpp"
#include "riskengine/core/rng/owen_scramble.hpp"
#include "riskengine/core/rng/sobol.hpp"
#include "riskengine/core/simulation.hpp"
#include "riskengine/core/stats/covariance.hpp"
#include "riskengine/methods/montecarlo/brownian_bridge.hpp"
#include "riskengine/methods/montecarlo/control_variates.hpp"

namespace riskengine {

enum class Sampling {
    PseudoRandom,  // Philox normals
    RandomizedQmc, // Owen-scrambled Sobol points, independent scramblings for the error bar
};

struct MonteCarloConfig {
    std::uint64_t paths;                // simulated paths (antithetic: pairs = paths / 2); QMC: per replication
    std::uint32_t steps = 1;            // time steps, equally spaced; also the fixing dates of a path payoff
    std::uint32_t blocks = 64;          // fixed work split: part of the experiment's identity (see BlockPlan)
    unsigned threads = 1;               // never changes the result
    bool antithetic = false;            // pair each path with its reflection z -> -z
    std::uint64_t pilot_paths = 10'000; // control variates only: pseudo-random paths used to estimate beta
    Sampling sampling = Sampling::PseudoRandom;
    std::uint32_t replications = 16;    // randomized QMC only: independent scramblings
    bool brownian_bridge = false;       // build paths by Brownian bridge (matters for QMC, steps > 1)
};

// Generic Monte Carlo pricer: any PathModel, any terminal or path payoff, optionally with
// antithetic variates, a control variate, randomized quasi-Monte Carlo and a Brownian bridge.
//
// Each path consumes steps x factors standard normals, step-major then factor-major (with a
// Brownian bridge, the normals are first mapped to increments factor by factor, in bridge order).
// - Pseudo-random: path i of block b draws its normals from RandomStream(key, b).
// - Randomized QMC: sample i uses Sobol point i, one coordinate per normal, Owen-scrambled with the
//   seeds of replication r (scramble_seeds(key, r, dims)) and inverted by norm_icdf. The estimate is
//   the mean of `replications` independent replications and its standard error is their standard
//   deviation over sqrt(replications): scrambling makes each replication an unbiased estimate, which
//   is what gives quasi-random numbers an error bar.
// The payoff is evaluated on the terminal spot (TerminalPayoff) or on the spot at the end of each
// step (PathPayoff), then discounted at the risk-free rate. With antithetic variates one sample is
// the average of a path and its reflection, so the standard error accounts for their correlation.
// With a control variate X of known mean mu, the samples are Y - beta (X - mu), where beta is
// estimated on an independent pseudo-random pilot run (estimating it on the same paths would bias
// the price). Because the pilot is independent, the estimate is unbiased for any beta, and the
// reported standard error is the one conditional on the pilot's beta; the sampling error of beta
// itself only affects the variance at second order.
// The price is a pure function of (market, key, config) minus `threads`.
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
        assert(config_.sampling == Sampling::PseudoRandom ||
               (config_.replications >= 2 && draws() <= sobol_data::kMaxDimension));
    }

    Estimate price(const MarketState& m, SeedKey key) const {
        assert((key.stream & kReservedStreamBits) == 0);
        const Model model(m);
        const double discount = std::exp(-m.rate.value * maturity_.value);

        double beta = 0.0, mu = 0.0;
        if constexpr (kHasControl) {
            const Covariance pilot =
                run<Covariance>(model, PseudoRandomNormals{SeedKey{key.seed, key.stream | kPilotStreamBit}},
                                config_.pilot_paths, [](Covariance& acc, double y, double x) { acc.add(x, y); });
            const double var_x = pilot.variance_x();
            beta = var_x > 0.0 ? pilot.covariance() / var_x : 0.0;
            mu = control_.expectation(m, maturity_, config_.steps);
        }
        const auto fold = [&](Welford& acc, double y, double x) {
            if constexpr (kHasControl) {
                acc.add(discount * y - beta * (discount * x - mu));
            } else {
                (void)x;
                acc.add(discount * y);
            }
        };

        if (config_.sampling == Sampling::PseudoRandom)
            return to_estimate(run<Welford>(model, PseudoRandomNormals{key}, config_.paths, fold));

        const Sobol sobol(draws());
        Welford replicates;
        for (std::uint32_t r = 0; r < config_.replications; ++r) {
            const std::vector<std::uint32_t> seeds = scramble_seeds(key, r, draws());
            replicates.add(run<Welford>(model, QuasiRandomNormals{&sobol, &seeds}, config_.paths, fold).mean());
        }
        Estimate e = to_estimate(replicates);
        e.samples = replicates.count() * samples_per_run(config_.paths);
        return e;
    }

    const MonteCarloConfig& config() const { return config_; }

private:
    // Sources of standard normals. Each block of a run gets its own cursor.
    struct PseudoRandomNormals {
        SeedKey key;

        struct Cursor {
            RandomStream rng;
            void fill(std::span<double> z) {
                for (double& v : z) v = rng.normal();
            }
        };
        Cursor cursor(std::uint32_t block, std::uint64_t) const { return {RandomStream(key, block)}; }
    };

    struct QuasiRandomNormals {
        const Sobol* sobol;
        const std::vector<std::uint32_t>* seeds;

        struct Cursor {
            SobolCursor points;
            const std::vector<std::uint32_t>* seeds;
            void fill(std::span<double> z) {
                const std::vector<std::uint32_t>& x = points.point();
                for (std::size_t j = 0; j < z.size(); ++j)
                    z[j] = norm_icdf(to_open_unit(owen_scramble(x[j], (*seeds)[j])));
                points.advance();
            }
        };
        Cursor cursor(std::uint32_t, std::uint64_t first_sample) const {
            return {SobolCursor(*sobol, first_sample), seeds};
        }
    };

    std::uint32_t draws() const { return config_.steps * static_cast<std::uint32_t>(Model::factors); }
    std::uint64_t samples_per_run(std::uint64_t paths) const { return config_.antithetic ? paths / 2 : paths; }

    // Simulates `paths` paths (pairs when antithetic) with normals from `source`, and folds each
    // sample's undiscounted payoff y and control value x (0 without a control) into an Acc.
    template <class Acc, class Source, class Fold>
    Acc run(const Model& model, const Source& source, std::uint64_t paths, Fold fold) const {
        const double dt = maturity_.value / config_.steps;
        const BlockPlan plan{samples_per_run(paths), config_.blocks, config_.threads};
        const bool bridge = config_.brownian_bridge && config_.steps > 1;
        const BrownianBridge bridge_map(config_.steps);

        return reduce_blocks(plan, [&](std::uint32_t b, std::uint64_t n) {
            auto normals = source.cursor(b, plan.block_start(b));
            std::vector<double> z(draws()), fixings(config_.steps), in(config_.steps), out(config_.steps);
            Acc acc;
            for (std::uint64_t i = 0; i < n; ++i) {
                normals.fill(z);
                if (bridge) to_increments(bridge_map, z, in, out);
                auto [y, x] = evaluate(model, dt, z, fixings);
                if (config_.antithetic) {
                    for (double& v : z) v = -v; // the bridge is linear: reflecting after it is the same
                    const auto [y2, x2] = evaluate(model, dt, z, fixings);
                    y = 0.5 * (y + y2);
                    x = 0.5 * (x + x2);
                }
                fold(acc, y, x);
            }
            return acc;
        });
    }

    // Brownian bridge, factor by factor: normals z[i * F + f] in bridge order become the increments
    // z[k * F + f] of step k.
    static void to_increments(const BrownianBridge& bridge, std::span<double> z, std::vector<double>& in,
                              std::vector<double>& out) {
        constexpr std::size_t factors = Model::factors;
        for (std::size_t f = 0; f < factors; ++f) {
            for (std::size_t i = 0; i < in.size(); ++i) in[i] = z[i * factors + f];
            bridge.transform(in, out);
            for (std::size_t k = 0; k < out.size(); ++k) z[k * factors + f] = out[k];
        }
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
