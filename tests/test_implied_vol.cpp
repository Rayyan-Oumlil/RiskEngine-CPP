#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "riskengine/methods/analytic/implied_vol.hpp"

using namespace riskengine;
using Catch::Matchers::WithinRel;

TEST_CASE("Implied vol recovers the canonical vol", "[iv]") {
    const VanillaOption call{Strike{100}, Maturity{1}, OptionType::Call};
    const auto res = implied_vol(call, 10.450583572185567, Spot{100}, Rate{0.05}, Rate{0});
    REQUIRE(res.status == ImpliedVolStatus::Ok);
    CHECK_THAT(res.vol, WithinRel(0.2, 1e-13));
}

TEST_CASE("Implied vol round trip over moneyness 0.5-2 and 1 day to 5 years", "[iv][roundtrip]") {
    // Phase 1 gate: price -> vol -> price on a grid covering the wings, for both the
    // out-of-the-money quote (tests the solver) and the in-the-money quote (tests conditioning).
    constexpr double eps = std::numeric_limits<double>::epsilon();
    const double s = 100, r = 0.03, q = 0.01;
    const std::vector<double> strikes = {50, 60, 75, 90, 100, 110, 125, 150, 200};
    const std::vector<double> maturities = {1.0 / 365, 7.0 / 365, 30.0 / 365, 0.25, 1, 2, 5};
    const std::vector<double> vols = {0.05, 0.1, 0.2, 0.5, 1.0};

    int otm_solved = 0, otm_underflow = 0, otm_total = 0;
    double worst_otm_error = 0.0;
    for (double k : strikes)
        for (double t : maturities)
            for (double sigma : vols) {
                const MarketState m{Spot{s}, Rate{r}, Rate{q}, Vol{sigma}};
                const bool call_is_otm = s * std::exp(-q * t) <= k * std::exp(-r * t);
                for (OptionType type : {OptionType::Call, OptionType::Put}) {
                    const bool otm = (type == OptionType::Call) == call_is_otm;
                    INFO("K=" << k << " T=" << t << " sigma=" << sigma
                              << (type == OptionType::Put ? " put" : " call") << (otm ? " (OTM)" : " (ITM)"));
                    const VanillaOption o{Strike{k}, Maturity{t}, type};
                    const double price = black_scholes_price(o, m);
                    const auto res = implied_vol(o, price, Spot{s}, Rate{r}, Rate{q});
                    REQUIRE(res.ok());
                    otm_total += otm ? 1 : 0;

                    const double vega = black_scholes_greeks(o, m).vega;
                    if (res.status == ImpliedVolStatus::ZeroTimeValue) {
                        // Allowed only when the time value is below what the quote can carry:
                        // underflow for an OTM quote, the rounding of the intrinsic for an ITM one.
                        if (otm) CHECK(price == 0.0);
                        else CHECK(vega * sigma <= 16 * eps * (s + k));
                        otm_underflow += otm ? 1 : 0;
                        continue;
                    }
                    const double vol_error = std::abs(res.vol - sigma) / sigma;
                    if (otm) {
                        CHECK(vol_error <= 1e-9);
                        worst_otm_error = std::max(worst_otm_error, vol_error);
                        ++otm_solved;
                    } else {
                        // Time value known to ~eps * (price + intrinsic terms): vol to that / vega.
                        CHECK(vol_error <= 1e-9 + 8 * eps * (price + s + k) / (sigma * vega));
                    }
                    const double repriced =
                        black_scholes_price(o, MarketState{Spot{s}, Rate{r}, Rate{q}, Vol{res.vol}});
                    CHECK(std::abs(repriced - price) <= 1e-13 * (s + k));
                }
            }
    // Every OTM quote with a representable price is solved; the rest underflowed to exactly 0
    // (e.g. 1 day, 5 % vol, K/S = 2), where there is no vol to recover.
    UNSCOPED_INFO("OTM solved " << otm_solved << "/" << otm_total << " (" << otm_underflow
                                << " underflowed to 0), worst OTM rel vol error " << worst_otm_error);
    CHECK(otm_solved + otm_underflow == otm_total);
    CHECK(otm_solved >= 280);
}

TEST_CASE("Implied vol reports no-arbitrage violations instead of returning garbage", "[iv][bounds]") {
    const double s = 100, r = 0.05, q = 0.0, t = 1.0;
    const VanillaOption call{Strike{80}, Maturity{t}, OptionType::Call};
    const VanillaOption put{Strike{120}, Maturity{t}, OptionType::Put};
    const double call_intrinsic = s - 80 * std::exp(-r * t);
    const double put_intrinsic = 120 * std::exp(-r * t) - s;

    CHECK(implied_vol(call, call_intrinsic - 0.01, Spot{s}, Rate{r}, Rate{q}).status == ImpliedVolStatus::BelowIntrinsic);
    CHECK(implied_vol(put, put_intrinsic - 0.01, Spot{s}, Rate{r}, Rate{q}).status == ImpliedVolStatus::BelowIntrinsic);
    CHECK(implied_vol(call, s + 1.0, Spot{s}, Rate{r}, Rate{q}).status == ImpliedVolStatus::AboveMaximum);
    CHECK(implied_vol(put, 120.0, Spot{s}, Rate{r}, Rate{q}).status == ImpliedVolStatus::AboveMaximum);
    CHECK(implied_vol(call, -1.0, Spot{s}, Rate{r}, Rate{q}).status == ImpliedVolStatus::BelowIntrinsic);

    const VanillaOption otm_call{Strike{120}, Maturity{t}, OptionType::Call};
    const auto zero = implied_vol(otm_call, 0.0, Spot{s}, Rate{r}, Rate{q});
    CHECK(zero.status == ImpliedVolStatus::ZeroTimeValue);
    CHECK(zero.vol == 0.0);
}

TEST_CASE("Implied vol handles very high vols without Newton-style divergence", "[iv]") {
    const double s = 100, r = 0.02, q = 0.0;
    for (double sigma : {2.0, 3.0, 5.0}) {
        const VanillaOption o{Strike{100}, Maturity{1}, OptionType::Call};
        const double price = black_scholes_price(o, MarketState{Spot{s}, Rate{r}, Rate{q}, Vol{sigma}});
        const auto res = implied_vol(o, price, Spot{s}, Rate{r}, Rate{q});
        REQUIRE(res.status == ImpliedVolStatus::Ok);
        CHECK_THAT(res.vol, WithinRel(sigma, 1e-8));
    }
}
