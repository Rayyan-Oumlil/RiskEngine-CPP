#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "riskengine/core/estimate.hpp"
#include "riskengine/core/market.hpp"
#include "riskengine/methods/analytic/black_scholes.hpp"
#include "riskengine/payoffs/vanilla.hpp"

namespace riskengine {

// Binomial trees for European and American vanilla options (report 4).
//
//   Crr            Cox-Ross-Rubinstein: u = e^{sigma sqrt(dt)}, d = 1/u. Converges as 1/n, but the error
//                  changes sign with the position of the strike between the terminal nodes, so it
//                  oscillates from one n to the next (for an at-the-money strike, even against odd n).
//   CrrAveraged    (V_n + V_{n+1}) / 2: cancels the even/odd part of that oscillation.
//   LeisenReimer   Leisen and Reimer (1996): u, d and p chosen by Peizer-Pratt inversion so that the
//                  tree is centred on the strike; odd n only, error O(1/n^2) for European options.
//   Bbs            Broadie and Detemple (1996): CRR in which the last step is replaced by the
//                  Black-Scholes price over one step (and, for an American option, the maximum of
//                  that price and the exercise value). Smooths the payoff kink: monotone 1/n error.
//   BbsRichardson  2 BBS(n) - BBS(n/2): Richardson extrapolation of the 1/n term; even n only. It is
//                  second order only where the BBS error is a smooth c/n: at the money its constant c
//                  differs between even and odd n, so n must be a multiple of 4 (report 4.2).
//
// The backward induction runs in place on one vector (O(n) memory), with the node spots taken from
// precomputed powers of u and d. Every method rejects parameters for which the up probability
// leaves (0, 1): such a tree prices with negative probabilities, i.e. admits arbitrage.
enum class TreeMethod { Crr, CrrAveraged, LeisenReimer, Bbs, BbsRichardson };
enum class Exercise { European, American };

inline const char* to_string(TreeMethod m) {
    switch (m) {
        case TreeMethod::Crr: return "crr";
        case TreeMethod::CrrAveraged: return "crr_averaged";
        case TreeMethod::LeisenReimer: return "leisen_reimer";
        case TreeMethod::Bbs: return "bbs";
        case TreeMethod::BbsRichardson: return "bbs_richardson";
    }
    return "unknown";
}

namespace detail {

struct Lattice {
    double u, d, p, discount; // up and down factors, up probability, one-step discount factor
};

// Peizer-Pratt method 2 inversion of the binomial distribution (Leisen and Reimer 1996).
inline double peizer_pratt(double z, unsigned n) {
    const double m = n + 1.0 / 3.0 + 0.1 / (n + 1.0);
    const double root = std::sqrt(1.0 - std::exp(-(z / m) * (z / m) * (n + 1.0 / 6.0)));
    return z >= 0.0 ? 0.5 + 0.5 * root : 0.5 - 0.5 * root;
}

inline Lattice make_lattice(TreeMethod method, double s, double k, double r, double q, double sigma, double t,
                            unsigned n) {
    const double dt = t / n;
    const double growth = std::exp((r - q) * dt);
    Lattice l{};
    l.discount = std::exp(-r * dt);
    if (method == TreeMethod::LeisenReimer) {
        const double sst = sigma * std::sqrt(t);
        const double d1 = (std::log(s / k) + (r - q) * t) / sst + 0.5 * sst;
        l.p = peizer_pratt(d1 - sst, n);
        l.u = growth * peizer_pratt(d1, n) / l.p;
        l.d = (growth - l.p * l.u) / (1.0 - l.p);
    } else {
        l.u = std::exp(sigma * std::sqrt(dt));
        l.d = 1.0 / l.u;
        l.p = (growth - l.d) / (l.u - l.d);
    }
    if (!(l.p > 0.0 && l.p < 1.0 && l.d > 0.0))
        throw std::invalid_argument("binomial tree: up probability " + std::to_string(l.p) +
                                    " outside (0, 1); use more steps (dt < sigma^2 / (r - q)^2)");
    return l;
}

// Backward induction through an n-step tree rooted at spot s; returns the value at the root. `bbs`
// replaces the last step by Black-Scholes over one step. If `nodes` is given, it receives the three
// values at step 2 (lowest spot first) and `spots` their spots: what the Greeks are read from.
inline double backward_induction(const VanillaOption& o, Exercise ex, const Lattice& l, double s, double r,
                                 double q, double sigma, double dt, unsigned n, bool bbs,
                                 std::vector<double>* nodes = nullptr, std::vector<double>* spots = nullptr) {
    const double k = o.strike.value;
    const bool call = o.type == OptionType::Call;
    const bool american = ex == Exercise::American;
    // pow_u[j] = u^j, pow_d[j] = d^j, each from one exp so rounding does not accumulate.
    std::vector<double> pow_u(n + 1), pow_d(n + 1);
    const double log_u = std::log(l.u), log_d = std::log(l.d);
    for (unsigned j = 0; j <= n; ++j) {
        pow_u[j] = std::exp(j * log_u);
        pow_d[j] = std::exp(j * log_d);
    }
    const auto spot = [&](unsigned i, unsigned j) { return s * pow_u[j] * pow_d[i - j]; };
    const auto intrinsic = [&](double x) { return std::max(call ? x - k : k - x, 0.0); };

    std::vector<double> v(n + 1);
    unsigned top = n; // the step whose values v holds
    if (bbs) {
        top = n - 1;
        for (unsigned j = 0; j <= top; ++j) {
            const double x = spot(top, j);
            const double continuation = bs_price(o.type, x, k, r, q, sigma, dt);
            v[j] = american ? std::max(continuation, intrinsic(x)) : continuation;
        }
    } else {
        for (unsigned j = 0; j <= n; ++j) v[j] = intrinsic(spot(n, j));
    }
    const double pu = l.discount * l.p, pd = l.discount * (1.0 - l.p);
    // Far out of the money, node values decay geometrically through the subnormal range, where
    // arithmetic is up to ~100 times slower on x86: a European call at n = 10,000 took 9 times
    // longer than its O(n^2) cost predicts (report 9). Such a value cannot move a price above
    // 1e-300, so it is flushed to zero.
    constexpr double kTiny = std::numeric_limits<double>::min();
    for (unsigned i = top; i-- > 0;) {
        for (unsigned j = 0; j <= i; ++j) {
            v[j] = pu * v[j + 1] + pd * v[j];
            if (v[j] < kTiny) v[j] = 0.0;
            if (american) v[j] = std::max(v[j], intrinsic(spot(i, j)));
        }
        if (i == 2 && nodes) {
            nodes->assign(v.begin(), v.begin() + 3);
            if (spots) *spots = {spot(2, 0), spot(2, 1), spot(2, 2)};
        }
    }
    return v[0];
}

inline void check_tree_inputs(const MarketState& m, const VanillaOption& o, unsigned n) {
    if (!(m.spot.value > 0.0 && o.strike.value > 0.0 && m.vol.value > 0.0 && o.maturity.value > 0.0))
        throw std::invalid_argument("binomial tree: needs spot, strike, vol and maturity > 0");
    if (n < 1) throw std::invalid_argument("binomial tree: needs at least one step");
}

// One plain tree (no averaging, no extrapolation).
inline double tree_price(TreeMethod base, const VanillaOption& o, Exercise ex, const MarketState& m, unsigned n) {
    const double s = m.spot.value, r = m.rate.value, q = m.div.value, sigma = m.vol.value, t = o.maturity.value;
    const Lattice l = make_lattice(base, s, o.strike.value, r, q, sigma, t, n);
    return backward_induction(o, ex, l, s, r, q, sigma, t / n, n, base == TreeMethod::Bbs);
}

} // namespace detail

// Price of a European or American vanilla option on a binomial tree with `steps` steps.
// LeisenReimer needs an odd number of steps and BbsRichardson an even one; both throw otherwise.
inline double binomial_price(TreeMethod method, const VanillaOption& o, Exercise ex, const MarketState& m,
                             unsigned steps) {
    detail::check_tree_inputs(m, o, steps);
    switch (method) {
        case TreeMethod::Crr:
        case TreeMethod::Bbs:
            return detail::tree_price(method, o, ex, m, steps);
        case TreeMethod::CrrAveraged:
            return 0.5 * (detail::tree_price(TreeMethod::Crr, o, ex, m, steps) +
                          detail::tree_price(TreeMethod::Crr, o, ex, m, steps + 1));
        case TreeMethod::LeisenReimer:
            if (steps % 2 == 0) throw std::invalid_argument("Leisen-Reimer needs an odd number of steps");
            return detail::tree_price(method, o, ex, m, steps);
        case TreeMethod::BbsRichardson:
            if (steps % 2 != 0) throw std::invalid_argument("BBS-Richardson needs an even number of steps");
            return 2.0 * detail::tree_price(TreeMethod::Bbs, o, ex, m, steps) -
                   detail::tree_price(TreeMethod::Bbs, o, ex, m, steps / 2);
    }
    throw std::invalid_argument("unknown tree method");
}

struct TreeGreeks {
    double price, delta, gamma, theta; // raw units of docs/conventions.md; theta per year
};

// Price, delta, gamma and theta read off the nodes of an extended CRR (or BBS) tree: the tree is
// rooted two steps before today (Pelsser and Vorst 1994), so that its step-2 nodes are
// S d^2, S, S u^2 today and the finite differences between them are centred on S. Theta compares
// today's middle node with the root, same spot, two steps earlier. The price is the n-step price.
inline TreeGreeks binomial_greeks(TreeMethod method, const VanillaOption& o, Exercise ex, const MarketState& m,
                                  unsigned steps) {
    detail::check_tree_inputs(m, o, steps);
    if (method != TreeMethod::Crr && method != TreeMethod::Bbs)
        throw std::invalid_argument(std::string("binomial_greeks: ") + to_string(method) +
                                    " is not supported (use crr or bbs)");
    const double s = m.spot.value, r = m.rate.value, q = m.div.value, sigma = m.vol.value, t = o.maturity.value;
    const double dt = t / steps;
    const unsigned n = steps + 2;
    const detail::Lattice l = detail::make_lattice(TreeMethod::Crr, s, o.strike.value, r, q, sigma, n * dt, n);
    std::vector<double> v, x;
    const double root =
        detail::backward_induction(o, ex, l, s, r, q, sigma, dt, n, method == TreeMethod::Bbs, &v, &x);
    const double delta = (v[2] - v[0]) / (x[2] - x[0]);
    const double gamma = ((v[2] - v[1]) / (x[2] - x[1]) - (v[1] - v[0]) / (x[1] - x[0])) / (0.5 * (x[2] - x[0]));
    return TreeGreeks{.price = v[1], .delta = delta, .gamma = gamma, .theta = (v[1] - root) / (2.0 * dt)};
}

// Pricer over a binomial tree. Satisfies the Pricer concept; the Estimate carries no standard
// error (the tree is deterministic) and its discretization error is the subject of report 4.
struct BinomialTree {
    VanillaOption option;
    Exercise exercise = Exercise::European;
    TreeMethod method = TreeMethod::Crr;
    unsigned steps = 1000;

    Estimate price(const MarketState& m) const {
        return Estimate{binomial_price(method, option, exercise, m, steps)};
    }
};

} // namespace riskengine
