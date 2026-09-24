#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "riskengine/greeks/dual.hpp"
#include "riskengine/models/gbm.hpp"
#include "riskengine/payoffs/digital.hpp"
#include "riskengine/payoffs/vanilla.hpp"

using namespace riskengine;

TEST_CASE("Dual numbers differentiate elementary expressions exactly", "[dual]") {
    for (double x : {0.3, 1.0, 2.5}) {
        INFO("x=" << x);
        // f(x) = x^2 / (1 + x) + exp(x) log(x) - sqrt(x)
        const Dual d = Dual::variable(x);
        const Dual f = d * d / (1.0 + d) + exp(d) * log(d) - sqrt(d);
        const double value = x * x / (1 + x) + std::exp(x) * std::log(x) - std::sqrt(x);
        const double slope = (x * x + 2 * x) / ((1 + x) * (1 + x)) + std::exp(x) * (std::log(x) + 1 / x) -
                             0.5 / std::sqrt(x);
        CHECK(std::abs(f.v - value) <= 1e-15 * std::abs(value) + 1e-15);
        CHECK(std::abs(f.d - slope) <= 1e-14 * std::abs(slope));
    }
    CHECK((-Dual::variable(2.0)).d == -1.0);
    CHECK((Dual(3.0) - Dual::variable(1.0)).d == -1.0);
}

TEST_CASE("Pathwise derivatives of the GBM terminal spot", "[dual]") {
    const double s = 100, r = 0.05, q = 0.02, sigma = 0.25, t = 0.75;
    for (double z : {-1.3, 0.0, 0.8}) {
        const double st = gbm_terminal_spot(s, r, q, sigma, t, z);
        const Dual by_spot = gbm_terminal_spot(Dual::variable(s), Dual(r), Dual(q), Dual(sigma), t, z);
        const Dual by_vol = gbm_terminal_spot(Dual(s), Dual(r), Dual(q), Dual::variable(sigma), t, z);
        CHECK(by_spot.v == st);
        CHECK(std::abs(by_spot.d - st / s) <= 1e-15 * st / s);                                    // dS_T/dS
        CHECK(std::abs(by_vol.d - st * (std::sqrt(t) * z - sigma * t)) <= 1e-13 * st);            // dS_T/dsigma
    }
}

TEST_CASE("Payoff derivative() agrees with automatic differentiation", "[dual]") {
    const VanillaPayoff call{100, OptionType::Call}, put{100, OptionType::Put};
    const StraddlePayoff straddle{100};
    const DigitalPayoff digital{100, OptionType::Call};
    for (double s : {80.0, 99.9, 100.1, 130.0}) {
        INFO("S=" << s);
        const Dual x = Dual::variable(s);
        CHECK(call(x).v == call(s));
        CHECK(call(x).d == call.derivative(s));
        CHECK(put(x).d == put.derivative(s));
        CHECK(straddle(x).d == straddle.derivative(s));
        CHECK(digital(x).v == digital(s));
        CHECK(digital(x).d == 0.0); // the jump is invisible pathwise
        CHECK(digital.derivative(s) == 0.0);
    }
}
