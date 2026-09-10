#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

#include "riskengine/methods/analytic/digital.hpp"
#include "riskengine/methods/analytic/geometric_asian.hpp"

using namespace riskengine;
using Catch::Matchers::WithinRel;

TEST_CASE("Digital matches high-precision reference values to 1e-9", "[digital][reference]") {
    // tools/bs_reference.py: mpmath at 50 digits, delta and gamma as numerical derivatives of the price.
    struct Ref {
        double k, r, q, sigma, t;
        OptionType type;
        double price, delta, gamma;
    };
    const Ref refs[] = {
        {100, 0.05, 0.0, 0.2, 1.0, OptionType::Call, 0.5323248154537634, 0.018762017345846894, -0.00032833530355232064},
        {100, 0.05, 0.0, 0.2, 1.0, OptionType::Put, 0.41890460904695061, -0.018762017345846894, 0.00032833530355232064},
        {95, 0.03, 0.02, 0.25, 0.5, OptionType::Call, 0.58217685761318141, 0.021651006619506537, -0.00049827310971151245},
        {95, 0.03, 0.02, 0.25, 0.5, OptionType::Put, 0.40293508198988125, -0.021651006619506537, 0.00049827310971151245},
    };
    for (const auto& ref : refs) {
        const VanillaOption o{Strike{ref.k}, Maturity{ref.t}, ref.type};
        const MarketState m{Spot{100}, Rate{ref.r}, Rate{ref.q}, Vol{ref.sigma}};
        CHECK_THAT(digital_price(o, m), WithinRel(ref.price, 1e-9));
        CHECK_THAT(digital_delta(o, m), WithinRel(ref.delta, 1e-9));
        CHECK_THAT(digital_gamma(o, m), WithinRel(ref.gamma, 1e-9));
    }
}

TEST_CASE("Digital call + put pays the discounted unit, delta matches finite differences", "[digital]") {
    const double r = 0.03, q = 0.01;
    for (double k : {50.0, 90.0, 100.0, 125.0, 200.0})
        for (double t : {0.1, 1.0, 5.0}) {
            const MarketState m{Spot{100}, Rate{r}, Rate{q}, Vol{0.3}};
            const VanillaOption call{Strike{k}, Maturity{t}, OptionType::Call};
            const VanillaOption put{Strike{k}, Maturity{t}, OptionType::Put};
            CHECK(std::abs(digital_price(call, m) + digital_price(put, m) - std::exp(-r * t)) <= 1e-15);

            const double h = 1e-4 * 100;
            const auto at = [&](double s) { return digital_price(call, MarketState{Spot{s}, Rate{r}, Rate{q}, Vol{0.3}}); };
            const double fd = (at(100 + h) - at(100 - h)) / (2 * h);
            CHECK(std::abs(digital_delta(call, m) - fd) <= 1e-6 * std::abs(fd) + 1e-10);
            CHECK(digital_delta(put, m) == -digital_delta(call, m));
        }
    const VanillaOption expired{Strike{90}, Maturity{0}, OptionType::Call};
    const MarketState m{Spot{100}, Rate{r}, Rate{q}, Vol{0.3}};
    CHECK(digital_price(expired, m) == 1.0);
    CHECK(digital_delta(expired, m) == 0.0);
    CHECK(digital_gamma(expired, m) == 0.0);
}

TEST_CASE("Geometric Asian with one fixing is the European option", "[asian]") {
    const MarketState m{Spot{100}, Rate{0.05}, Rate{0.02}, Vol{0.25}};
    for (double k : {80.0, 100.0, 120.0})
        for (OptionType type : {OptionType::Call, OptionType::Put}) {
            const double european = black_scholes_price(VanillaOption{Strike{k}, Maturity{1.0}, type}, m);
            CHECK_THAT(geometric_asian_price(type, Strike{k}, Maturity{1.0}, 1, m), WithinRel(european, 1e-13));
        }
}

TEST_CASE("Geometric Asian converges to the continuous-monitoring formula", "[asian]") {
    // Continuous geometric average (Kemna-Vorst): ln G ~ N(ln S + (r - q - sigma^2/2) T/2, sigma^2 T/3).
    // Written out separately here, so the discrete formula is checked against an independent limit.
    const double s = 100, r = 0.05, q = 0.02, sigma = 0.25, t = 1.0, k = 100;
    const double var = sigma * sigma * t / 3.0;
    const double mean = std::log(s) + (r - q - 0.5 * sigma * sigma) * t / 2.0;
    const double d1 = (mean - std::log(k) + var) / std::sqrt(var), d2 = d1 - std::sqrt(var);
    const double continuous = std::exp(-r * t) * (std::exp(mean + 0.5 * var) * norm_cdf(d1) - k * norm_cdf(d2));

    const MarketState m{Spot{s}, Rate{r}, Rate{q}, Vol{sigma}};
    const double discrete = geometric_asian_price(OptionType::Call, Strike{k}, Maturity{t}, 1'000'000, m);
    CHECK_THAT(discrete, WithinRel(continuous, 1e-5));
    // Monitoring more often averages over more of the path: less variance, a cheaper ATM call.
    CHECK(geometric_asian_price(OptionType::Call, Strike{k}, Maturity{t}, 12, m) <
          geometric_asian_price(OptionType::Call, Strike{k}, Maturity{t}, 4, m));
}
