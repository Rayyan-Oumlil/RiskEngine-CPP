#pragma once

#include <concepts>
#include <cstddef>
#include <span>

#include "riskengine/core/estimate.hpp"
#include "riskengine/core/market.hpp"
#include "riskengine/core/rng/random_stream.hpp"

namespace riskengine {

// A pricer owns its instrument and prices it against a market state.
template <class P>
concept Pricer = requires(const P& p, const MarketState& m) {
    { p.price(m) } -> std::same_as<Estimate>;
};

// A stochastic pricer is a pure function of (market, seed key): no mutable RNG state, so calling
// it again with the same key after bumping the market gives exact common random numbers.
template <class P>
concept StochasticPricer = requires(const P& p, const MarketState& m, SeedKey k) {
    { p.price(m, k) } -> std::same_as<Estimate>;
};

// A model with parameters beyond the market state (stochastic volatility, jumps) declares them as
// M::Params and is built from (market, params); GBM is built from the market alone.
template <class M>
concept ParametrizedModel = requires { typename M::Params; } &&
                            std::constructible_from<M, const MarketState&, const typename M::Params&>;

struct NoModelParams {};

template <class M>
struct model_params {
    using type = NoModelParams;
};
template <ParametrizedModel M>
struct model_params<M> {
    using type = typename M::Params;
};
template <class M>
using model_params_t = typename model_params<M>::type;

// Stochastic dynamics, decoupled from the solvers. `factors` normals drive one step; spot()
// exposes the traded price from the (possibly multi-dimensional) state.
template <class M>
concept PathModel = (std::constructible_from<M, const MarketState&> || ParametrizedModel<M>) &&
                    requires(const M& m, typename M::State s, double dt, std::span<const double> z) {
                        { M::factors } -> std::convertible_to<std::size_t>;
                        { m.initial_state() } -> std::same_as<typename M::State>;
                        { m.step(s, dt, z) } -> std::same_as<typename M::State>;
                        { m.spot(s) } -> std::convertible_to<double>;
                    };

// Payoffs on the terminal spot, and on the fixings S(t_1), ..., S(t_n) of a path.
template <class F>
concept TerminalPayoff = requires(const F& f, double s) {
    { f(s) } -> std::convertible_to<double>;
};

template <class F>
concept PathPayoff = requires(const F& f, std::span<const double> fixings) {
    { f(fixings) } -> std::convertible_to<double>;
};

} // namespace riskengine
