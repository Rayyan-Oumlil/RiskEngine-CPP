#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "riskengine/concepts.hpp"
#include "riskengine/core/estimate.hpp"
#include "riskengine/core/market.hpp"
#include "riskengine/core/rng/random_stream.hpp"
#include "riskengine/core/stats/welford.hpp"
#include "riskengine/payoffs/vanilla.hpp"

namespace riskengine {

struct LongstaffSchwartzConfig {
    std::uint64_t paths;      // simulated paths; half are held out for the exercise-boundary fit
                               // and half for an unbiased re-pricing pass (Longstaff & Schwartz 2001)
    std::uint32_t steps = 50; // equally spaced exercise dates, including maturity
};

namespace detail {

// Fits y ~ b0 + b1 x + b2 x^2 by ordinary least squares (3x3 normal equations, solved directly:
// the basis is fixed and tiny, so a closed-form Cramer solve is simpler and just as accurate as an
// iterative or QR solver here). Returns {0, 0, 0} if fewer than 3 points or the design is singular
// (e.g. every in-the-money path has the same spot), so the continuation value falls back to 0 and
// that path is never exercised early on that step - conservative, not a NaN.
struct QuadraticFit {
    double b0 = 0.0, b1 = 0.0, b2 = 0.0;

    double operator()(double x) const { return b0 + b1 * x + b2 * x * x; }
};

inline QuadraticFit fit_quadratic(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() < 3) return {};

    double s0 = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0; // sum x^0 .. x^4
    double t0 = 0, t1 = 0, t2 = 0;                 // sum y, x*y, x^2*y
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double xi = x[i], yi = y[i];
        const double xi2 = xi * xi;
        s0 += 1.0;
        s1 += xi;
        s2 += xi2;
        s3 += xi2 * xi;
        s4 += xi2 * xi2;
        t0 += yi;
        t1 += xi * yi;
        t2 += xi2 * yi;
    }

    // Solve [s0 s1 s2; s1 s2 s3; s2 s3 s4] [b0 b1 b2]^T = [t0 t1 t2]^T by Cramer's rule.
    const double det = s0 * (s2 * s4 - s3 * s3) - s1 * (s1 * s4 - s3 * s2) + s2 * (s1 * s3 - s2 * s2);
    if (!(std::abs(det) > 0.0)) return {}; // singular (or all-NaN): no fit, treat as never-exercise

    const double det0 = t0 * (s2 * s4 - s3 * s3) - s1 * (t1 * s4 - s3 * t2) + s2 * (t1 * s3 - s2 * t2);
    const double det1 = s0 * (t1 * s4 - s3 * t2) - t0 * (s1 * s4 - s3 * s2) + s2 * (s1 * t2 - t1 * s2);
    const double det2 = s0 * (s2 * t2 - t1 * s3) - s1 * (s1 * t2 - t1 * s2) + t0 * (s1 * s3 - s2 * s2);

    return QuadraticFit{det0 / det, det1 / det, det2 / det};
}

} // namespace detail

// American option pricing by least-squares Monte Carlo (Longstaff & Schwartz 2001, "Valuing
// American Options by Simulation: A Simple Least-Squares Approach", RFS 14(1)).
//
// Unlike MonteCarlo<Model, Payoff> (methods/montecarlo/engine.hpp), which folds each path into a
// running accumulator and never stores it, this method needs every path's full trajectory at once:
// the continuation value at each exercise date is estimated by regressing, across all paths
// in-the-money at that date, the payoff each path actually received later against its own spot at
// that date. This is fundamentally a backward pass over stored paths, not a streaming fold, so it
// lives in its own class with O(paths x steps) memory instead of the O(1) engine above.
//
// Basis: {1, S, S^2} on the paths in the money at each date, the standard "regress on the spot and
// its square, ignore out-of-the-money paths" baseline used in most textbook and practitioner LSM
// implementations (the original paper suggests Laguerre polynomials; the two choices agree closely
// near the money and 1/S/S^2 is simpler to reason about and to check by hand).
//
// Estimate::discretization on the returned price is the well-known low bias of the regression
// estimator (Longstaff & Schwartz 2001, section 4): the fitted continuation value is only an
// approximation, so the exercise decision it drives is never better than optimal, which pushes the
// price down. It is reported as (in-sample regression estimate) - (out-of-sample price on the fitted
// boundary), not folded into std_error, so a caller can tell statistical noise from this bias.
template <PathModel Model>
class LongstaffSchwartz {
public:
    LongstaffSchwartz(VanillaOption option, LongstaffSchwartzConfig config)
        requires(!ParametrizedModel<Model>)
        : option_(option), config_(config) {
        assert(config_.steps >= 1);
        assert(Model::factors == 1); // the regression basis is a function of the single spot
    }

    Estimate price(const MarketState& m, SeedKey key) const {
        const Model model(m);
        const double dt = option_.maturity.value / config_.steps;
        const double discount_step = std::exp(-m.rate.value * dt);
        const VanillaPayoff payoff{option_.strike.value, option_.type};

        // Two independent halves, per Longstaff & Schwartz 2001 section 3: fitting and pricing the
        // exercise boundary on the same paths overstates the value (the regression can "see" a
        // path's own future when deciding whether it should have exercised). std_error and
        // discretization below come from the second, out-of-sample half.
        const std::uint64_t half = config_.paths / 2;
        const std::vector<std::vector<double>> fit_paths = simulate_paths(model, dt, key, half, 0);
        const std::vector<std::vector<double>> price_paths =
            simulate_paths(model, dt, key, config_.paths - half, half);

        const std::vector<detail::QuadraticFit> boundary = fit_exercise_boundary(fit_paths, payoff, discount_step);

        const Estimate in_sample = price_along_boundary(fit_paths, boundary, payoff, discount_step);
        const Estimate out_sample = price_along_boundary(price_paths, boundary, payoff, discount_step);

        Estimate result = out_sample;
        result.discretization = in_sample.value - out_sample.value; // regression's look-ahead bias
        result.samples = fit_paths.size() + price_paths.size();
        return result;
    }

private:
    // path[i][k] is the spot of path i at exercise date k (k = 0 .. steps - 1); date steps - 1 is
    // maturity. Rows generated from block `block_offset + i`, so the fit and pricing halves are
    // statistically independent draws under one key.
    std::vector<std::vector<double>> simulate_paths(const Model& model, double dt, SeedKey key,
                                                     std::uint64_t n, std::uint64_t block_offset) const {
        std::vector<std::vector<double>> paths(n, std::vector<double>(config_.steps));
        for (std::uint64_t i = 0; i < n; ++i) {
            RandomStream rng(key, static_cast<std::uint32_t>(block_offset + i));
            typename Model::State s = model.initial_state();
            for (std::uint32_t k = 0; k < config_.steps; ++k) {
                const double z = rng.normal();
                s = model.step(s, dt, std::span<const double>(&z, 1));
                paths[i][k] = model.spot(s);
            }
        }
        return paths;
    }

    // Backward induction: at each date but the last, regress the discounted cash flow each
    // in-the-money path actually received afterward on its spot, then mark early exercise wherever
    // the immediate payoff beats the fitted continuation value. Returns one fit per exercise date
    // before maturity, indexed by date (a default-constructed fit, "never exercise", at a date with
    // fewer than 3 in-the-money paths to regress on).
    std::vector<detail::QuadraticFit> fit_exercise_boundary(const std::vector<std::vector<double>>& paths,
                                                             const VanillaPayoff& payoff,
                                                             double discount_step) const {
        const std::size_t n = paths.size();
        std::vector<double> cash_flow(n);
        std::vector<std::uint32_t> exercised_at(n, config_.steps - 1); // default: only at maturity
        for (std::size_t i = 0; i < n; ++i) cash_flow[i] = payoff(paths[i].back());

        std::vector<detail::QuadraticFit> fits(config_.steps - 1);

        for (std::uint32_t k = config_.steps - 1; k-- > 0;) {
            std::vector<double> itm_spot, itm_continuation;
            for (std::size_t i = 0; i < n; ++i) {
                const double immediate = payoff(paths[i][k]);
                if (immediate <= 0.0) continue; // out of the money: never optimal to exercise here
                const double steps_ahead = static_cast<double>(exercised_at[i] - k);
                const double discounted_future = cash_flow[i] * std::pow(discount_step, steps_ahead);
                itm_spot.push_back(paths[i][k]);
                itm_continuation.push_back(discounted_future);
            }

            const detail::QuadraticFit fit = detail::fit_quadratic(itm_spot, itm_continuation);
            fits[k] = fit;

            for (std::size_t i = 0; i < n; ++i) {
                const double immediate = payoff(paths[i][k]);
                if (immediate <= 0.0) continue;
                if (immediate > fit(paths[i][k])) {
                    cash_flow[i] = immediate;
                    exercised_at[i] = k;
                }
            }
        }

        return fits;
    }

    // Prices a set of paths along an already-fitted exercise boundary: walk forward, exercise the
    // first date the immediate payoff beats the fitted continuation value (or hold to maturity),
    // discount that one cash flow back to time 0. This is the out-of-sample step of Longstaff &
    // Schwartz 2001 section 3, and what makes the reported price and std_error unbiased.
    Estimate price_along_boundary(const std::vector<std::vector<double>>& paths,
                                  const std::vector<detail::QuadraticFit>& boundary, const VanillaPayoff& payoff,
                                  double discount_step) const {
        Welford acc;
        for (const std::vector<double>& path : paths) {
            double pv = 0.0;
            bool exercised = false;
            for (std::uint32_t k = 0; k + 1 < config_.steps; ++k) {
                const double immediate = payoff(path[k]);
                if (immediate > 0.0 && k < boundary.size() && boundary[k](path[k]) < immediate) {
                    pv = immediate * std::pow(discount_step, static_cast<double>(k + 1));
                    exercised = true;
                    break;
                }
            }
            if (!exercised) {
                const double maturity_payoff = payoff(path.back());
                pv = maturity_payoff * std::pow(discount_step, static_cast<double>(config_.steps));
            }
            acc.add(pv);
        }
        return to_estimate(acc);
    }

    VanillaOption option_;
    LongstaffSchwartzConfig config_;
};

} // namespace riskengine
