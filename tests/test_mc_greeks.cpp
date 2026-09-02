#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <stdexcept>
#include <utility>

#include "riskengine/methods/analytic/digital.hpp"
#include "riskengine/methods/montecarlo/greeks.hpp"
#include "riskengine/payoffs/digital.hpp"
#include "riskengine/payoffs/vanilla.hpp"

using namespace riskengine;

static_assert(DifferentiableTerminalPayoff<VanillaPayoff> && DifferentiableTerminalPayoff<DigitalPayoff>);

namespace {

const MarketState kMarket{Spot{100}, Rate{0.05}, Rate{0.0}, Vol{0.2}};
const Maturity kT{1.0};
const VanillaOption kCall{Strike{100}, kT, OptionType::Call};
const VanillaPayoff kCallPayoff{100, OptionType::Call};
const DigitalPayoff kDigitalPayoff{100, OptionType::Call};

void check_within_4se(const Estimate& e, double exact) {
    INFO("MC " << e.value << " +- " << e.std_error << ", exact " << exact);
    CHECK(std::abs(e.value - exact) < 4.0 * e.std_error);
}

MarketState at_spot(double s) { return MarketState{Spot{s}, kMarket.rate, kMarket.div, kMarket.vol}; }

} // namespace

TEST_CASE("Unbiased estimators agree with the closed-form Greeks within 4 SE", "[mc_greeks]") {
    const GreekConfig cfg{.paths = 200'000, .threads = 4};
    const Greeks bs = black_scholes_greeks(kCall, kMarket);
    for (GreekMethod m : {GreekMethod::Pathwise, GreekMethod::PathwiseDual, GreekMethod::LikelihoodRatio}) {
        INFO(to_string(m));
        check_within_4se(mc_greek(Greek::Delta, m, kCallPayoff, kT, kMarket, SeedKey{1}, cfg), bs.delta);
    }
    for (GreekMethod m : {GreekMethod::LikelihoodRatio, GreekMethod::Mixed}) {
        INFO(to_string(m));
        check_within_4se(mc_greek(Greek::Gamma, m, kCallPayoff, kT, kMarket, SeedKey{2}, cfg), bs.gamma);
    }
    check_within_4se(mc_greek(Greek::Delta, GreekMethod::LikelihoodRatio, kDigitalPayoff, kT, kMarket, SeedKey{3}, cfg),
                     digital_delta(kCall, kMarket));
    check_within_4se(mc_greek(Greek::Gamma, GreekMethod::LikelihoodRatio, kDigitalPayoff, kT, kMarket, SeedKey{4}, cfg),
                     digital_gamma(kCall, kMarket));
}

TEST_CASE("Finite differences are unbiased for the finite-difference quotient of the prices", "[mc_greeks]") {
    // E[FD estimator] is exactly the same central difference applied to the true prices, at the
    // same realized bumps: that is its target, whatever h.
    for (double h : {0.2, 0.05, 0.01}) {
        const GreekConfig cfg{.paths = 200'000, .bump = h, .threads = 4};
        const double up = 100 * (1 + h), down = 100 * (1 - h), half = 0.5 * (up - down);
        for (bool digital : {false, true}) {
            INFO("h=" << h << (digital ? " digital" : " call"));
            const auto price = [&](double s) {
                return digital ? digital_price(kCall, at_spot(s)) : black_scholes_price(kCall, at_spot(s));
            };
            const double delta_quotient = (price(up) - price(down)) / (up - down);
            const double gamma_quotient = (price(up) - 2 * price(100) + price(down)) / (half * half);
            for (GreekMethod m : {GreekMethod::FdIndependent, GreekMethod::FdCrn}) {
                INFO(to_string(m));
                const auto delta = digital ? mc_greek(Greek::Delta, m, kDigitalPayoff, kT, kMarket, SeedKey{5}, cfg)
                                           : mc_greek(Greek::Delta, m, kCallPayoff, kT, kMarket, SeedKey{5}, cfg);
                const auto gamma = digital ? mc_greek(Greek::Gamma, m, kDigitalPayoff, kT, kMarket, SeedKey{6}, cfg)
                                           : mc_greek(Greek::Gamma, m, kCallPayoff, kT, kMarket, SeedKey{6}, cfg);
                check_within_4se(delta, delta_quotient);
                check_within_4se(gamma, gamma_quotient);
            }
        }
    }
}

TEST_CASE("Pathwise on a digital converges, with zero variance, to the wrong answer", "[mc_greeks]") {
    const GreekConfig cfg{.paths = 100'000, .threads = 4};
    for (GreekMethod m : {GreekMethod::Pathwise, GreekMethod::PathwiseDual}) {
        const Estimate e = mc_greek(Greek::Delta, m, kDigitalPayoff, kT, kMarket, SeedKey{7}, cfg);
        CHECK(e.value == 0.0);
        CHECK(e.std_error == 0.0);
    }
    const Estimate mixed = mc_greek(Greek::Gamma, GreekMethod::Mixed, kDigitalPayoff, kT, kMarket, SeedKey{7}, cfg);
    CHECK(mixed.value == 0.0);
    CHECK(digital_delta(kCall, kMarket) > 0.018); // while the true delta is 0.0188
}

TEST_CASE("Automatic differentiation reproduces the hand-written pathwise delta", "[mc_greeks][dual]") {
    const GreekConfig cfg{.paths = 50'000, .threads = 2};
    for (double k : {80.0, 100.0, 130.0}) {
        const VanillaPayoff call{k, OptionType::Call};
        const Estimate hand = mc_greek(Greek::Delta, GreekMethod::Pathwise, call, kT, kMarket, SeedKey{8}, cfg);
        const Estimate dual = mc_greek(Greek::Delta, GreekMethod::PathwiseDual, call, kT, kMarket, SeedKey{8}, cfg);
        CHECK(std::abs(hand.value - dual.value) <= 1e-14 * hand.value);
        CHECK(std::abs(hand.std_error - dual.std_error) <= 1e-12 * hand.std_error);
    }
}

TEST_CASE("Common random numbers: what they fix and what they do not", "[mc_greeks]") {
    const auto se = [&](Greek g, GreekMethod m, const auto& payoff, double h) {
        return mc_greek(g, m, payoff, kT, kMarket, SeedKey{9}, GreekConfig{.paths = 100'000, .bump = h, .threads = 4})
            .std_error;
    };
    // Call delta: CRN removes the 1/h blow-up of independent seeds (at h = 1 % the SE falls ~18x);
    // the SE then no longer depends on h.
    CHECK(se(Greek::Delta, GreekMethod::FdCrn, kCallPayoff, 0.01) <
          se(Greek::Delta, GreekMethod::FdIndependent, kCallPayoff, 0.01) / 10.0);
    CHECK(se(Greek::Delta, GreekMethod::FdCrn, kCallPayoff, 0.001) < 1.5 * se(Greek::Delta, GreekMethod::FdCrn, kCallPayoff, 0.1));
    // Digital delta and call gamma: even with CRN the SE grows as h^{-1/2} (a factor ~10 over h / 100).
    CHECK(se(Greek::Delta, GreekMethod::FdCrn, kDigitalPayoff, 0.001) > 5.0 * se(Greek::Delta, GreekMethod::FdCrn, kDigitalPayoff, 0.1));
    CHECK(se(Greek::Gamma, GreekMethod::FdCrn, kCallPayoff, 0.001) > 5.0 * se(Greek::Gamma, GreekMethod::FdCrn, kCallPayoff, 0.1));
}

TEST_CASE("Greek estimates are thread-invariant", "[mc_greeks][determinism]") {
    for (GreekMethod m : {GreekMethod::FdIndependent, GreekMethod::FdCrn, GreekMethod::LikelihoodRatio}) {
        const Estimate a = mc_greek(Greek::Gamma, m, kCallPayoff, kT, kMarket, SeedKey{10}, GreekConfig{.paths = 20'000, .threads = 1});
        const Estimate b = mc_greek(Greek::Gamma, m, kCallPayoff, kT, kMarket, SeedKey{10}, GreekConfig{.paths = 20'000, .threads = 8});
        CHECK(a.value == b.value);
        CHECK(a.std_error == b.std_error);
    }
    CHECK(!supports(GreekMethod::Pathwise, Greek::Gamma));
    CHECK(!supports(GreekMethod::Mixed, Greek::Delta));
}

TEST_CASE("Unsupported estimator / Greek pairs are rejected in every build", "[mc_greeks]") {
    const GreekConfig cfg{.paths = 100};
    for (auto [greek, method] : {std::pair{Greek::Gamma, GreekMethod::Pathwise}, std::pair{Greek::Gamma, GreekMethod::PathwiseDual},
                                 std::pair{Greek::Delta, GreekMethod::Mixed}})
        CHECK_THROWS_AS(mc_greek(greek, method, kCallPayoff, kT, kMarket, SeedKey{11}, cfg), std::invalid_argument);
}
