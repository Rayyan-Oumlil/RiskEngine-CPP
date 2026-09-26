#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <stdexcept>

#include "riskengine/concepts.hpp"
#include "riskengine/methods/analytic/black_scholes.hpp"
#include "riskengine/methods/tree/binomial.hpp"

using namespace riskengine;

static_assert(Pricer<BinomialTree>);

namespace {

const MarketState kMarket{Spot{100}, Rate{0.05}, Rate{0.0}, Vol{0.2}};
const VanillaOption kCall{Strike{100}, Maturity{1}, OptionType::Call};
const VanillaOption kPut{Strike{100}, Maturity{1}, OptionType::Put};

double error(TreeMethod method, const VanillaOption& o, unsigned n, const MarketState& m = kMarket) {
    return binomial_price(method, o, Exercise::European, m, n) - black_scholes_price(o, m);
}

} // namespace

TEST_CASE("Invalid trees are rejected, never priced", "[tree]") {
    // dt = 1 > sigma^2 / (r - q)^2 = 0.0004: the CRR up probability exceeds 1.
    const MarketState steep{Spot{100}, Rate{0.5}, Rate{0.0}, Vol{0.01}};
    CHECK_THROWS_AS(binomial_price(TreeMethod::Crr, kCall, Exercise::European, steep, 1), std::invalid_argument);
    CHECK_NOTHROW(binomial_price(TreeMethod::Crr, kCall, Exercise::European, steep, 3000));
    CHECK_THROWS_AS(binomial_price(TreeMethod::LeisenReimer, kCall, Exercise::European, kMarket, 100),
                    std::invalid_argument);
    CHECK_THROWS_AS(binomial_price(TreeMethod::BbsRichardson, kCall, Exercise::European, kMarket, 101),
                    std::invalid_argument);
    CHECK_THROWS_AS(binomial_price(TreeMethod::Crr, kCall, Exercise::European, kMarket, 0), std::invalid_argument);
    const MarketState flat{Spot{100}, Rate{0.05}, Rate{0.0}, Vol{0.0}};
    CHECK_THROWS_AS(binomial_price(TreeMethod::Crr, kCall, Exercise::European, flat, 10), std::invalid_argument);
    CHECK_THROWS_AS(binomial_greeks(TreeMethod::LeisenReimer, kCall, Exercise::European, kMarket, 101),
                    std::invalid_argument);
}

TEST_CASE("Reference values at n = 20,000 (docs/model_risk_report.md, Appendix C)", "[tree][american]") {
    // CRR, recomputed independently (docs/model_risk_report.md, Appendix C).
    CHECK(std::abs(binomial_price(TreeMethod::Crr, kPut, Exercise::European, kMarket, 20000) - 5.573426) < 5e-7);
    CHECK(std::abs(binomial_price(TreeMethod::Crr, kPut, Exercise::American, kMarket, 20000) - 6.090333) < 5e-7);
}

TEST_CASE("European trees converge to Black-Scholes at their theoretical order", "[tree]") {
    // CRR: O(1/n), with the sign of the error alternating between even and odd n at the money.
    CHECK(error(TreeMethod::Crr, kCall, 100) < 0.0);
    CHECK(error(TreeMethod::Crr, kCall, 101) > 0.0);
    const double crr_ratio = error(TreeMethod::Crr, kCall, 1000) / error(TreeMethod::Crr, kCall, 100);
    CHECK(crr_ratio > 0.09);
    CHECK(crr_ratio < 0.11);
    // Leisen-Reimer: O(1/n^2), so ten times the steps divides the error by about a hundred.
    const double lr_ratio = error(TreeMethod::LeisenReimer, kCall, 1001) / error(TreeMethod::LeisenReimer, kCall, 101);
    CHECK(lr_ratio > 0.008);
    CHECK(lr_ratio < 0.012);
    CHECK(std::abs(error(TreeMethod::LeisenReimer, kCall, 1001)) < 1e-6);
    // BBS is monotone in n; Richardson removes its 1/n term.
    CHECK(error(TreeMethod::Bbs, kCall, 1000) > 0.0);
    CHECK(error(TreeMethod::Bbs, kCall, 1001) > 0.0);
    CHECK(std::abs(error(TreeMethod::BbsRichardson, kCall, 1000)) < 1e-6);
    // Averaging consecutive CRR trees cancels most of the even/odd oscillation.
    CHECK(std::abs(error(TreeMethod::CrrAveraged, kCall, 1000)) < 0.1 * std::abs(error(TreeMethod::Crr, kCall, 1000)));
}

TEST_CASE("European trees satisfy put-call parity", "[tree]") {
    const MarketState m{Spot{100}, Rate{0.05}, Rate{0.02}, Vol{0.3}};
    const VanillaOption call{Strike{110}, Maturity{0.5}, OptionType::Call};
    const VanillaOption put{Strike{110}, Maturity{0.5}, OptionType::Put};
    const double forward = 100 * std::exp(-0.02 * 0.5) - 110 * std::exp(-0.05 * 0.5);
    for (TreeMethod method : {TreeMethod::Crr, TreeMethod::LeisenReimer, TreeMethod::Bbs}) {
        INFO(to_string(method));
        const double c = binomial_price(method, call, Exercise::European, m, 301);
        const double p = binomial_price(method, put, Exercise::European, m, 301);
        CHECK(std::abs(c - p - forward) < 1e-11);
    }
}

TEST_CASE("American invariants", "[tree][american]") {
    for (TreeMethod method : {TreeMethod::Crr, TreeMethod::LeisenReimer, TreeMethod::Bbs}) {
        INFO(to_string(method));
        // Without dividends an American call is never exercised early (Merton): the tree gives the
        // European value exactly, since the continuation value always exceeds the exercise value.
        CHECK(binomial_price(method, kCall, Exercise::American, kMarket, 501) ==
              binomial_price(method, kCall, Exercise::European, kMarket, 501));
        // An American put is worth more than the European one, and at least its exercise value.
        CHECK(binomial_price(method, kPut, Exercise::American, kMarket, 501) >
              binomial_price(method, kPut, Exercise::European, kMarket, 501) + 0.5);
        const VanillaOption deep{Strike{150}, Maturity{1}, OptionType::Put};
        CHECK(binomial_price(method, deep, Exercise::American, kMarket, 501) >= 50.0);
        // With a dividend yield, early exercise of the call has value.
        const MarketState div{Spot{100}, Rate{0.05}, Rate{0.08}, Vol{0.2}};
        CHECK(binomial_price(method, kCall, Exercise::American, div, 501) >
              binomial_price(method, kCall, Exercise::European, div, 501) + 1e-3);
    }
    // The American put methods agree with each other to their discretization error.
    const double lr = binomial_price(TreeMethod::LeisenReimer, kPut, Exercise::American, kMarket, 4001);
    const double bbsr = binomial_price(TreeMethod::BbsRichardson, kPut, Exercise::American, kMarket, 4000);
    CHECK(std::abs(lr - bbsr) < 2e-4);
}

TEST_CASE("Greeks from the extended tree match Black-Scholes", "[tree][greeks]") {
    for (const VanillaOption& o : {kCall, kPut}) {
        const Greeks bs = black_scholes_greeks(o, kMarket);
        for (TreeMethod method : {TreeMethod::Crr, TreeMethod::Bbs}) {
            INFO(to_string(method));
            const TreeGreeks g = binomial_greeks(method, o, Exercise::European, kMarket, 2000);
            // The price is the one of the plain n-step tree.
            CHECK(std::abs(g.price - binomial_price(method, o, Exercise::European, kMarket, 2000)) < 1e-11);
            CHECK(std::abs(g.delta - bs.delta) < 1e-4);
            CHECK(std::abs(g.gamma / bs.gamma - 1.0) < 1e-3);
            CHECK(std::abs(g.theta / bs.theta - 1.0) < 1e-3);
        }
    }
    // American put: a larger delta in absolute value and a larger gamma than the European one.
    const TreeGreeks am = binomial_greeks(TreeMethod::Crr, kPut, Exercise::American, kMarket, 2000);
    const Greeks eu = black_scholes_greeks(kPut, kMarket);
    CHECK(am.delta < eu.delta);
    CHECK(am.gamma > eu.gamma);
}

TEST_CASE("The tree pricer is a Pricer", "[tree]") {
    const BinomialTree tree{kPut, Exercise::American, TreeMethod::LeisenReimer, 1001};
    const Estimate e = tree.price(kMarket);
    CHECK(e.value == binomial_price(TreeMethod::LeisenReimer, kPut, Exercise::American, kMarket, 1001));
    CHECK(e.std_error == 0.0);
}
