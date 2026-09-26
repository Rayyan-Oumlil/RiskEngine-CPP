#include <catch2/catch_test_macros.hpp>

#include <bit>
#include <cmath>
#include <cstdint>

#include "riskengine/methods/analytic/black_scholes.hpp"
#include "riskengine/methods/montecarlo/longstaff_schwartz.hpp"
#include "riskengine/methods/tree/binomial.hpp"
#include "riskengine/models/gbm.hpp"

using namespace riskengine;

namespace {

const MarketState kMarket{Spot{100}, Rate{0.05}, Rate{0.0}, Vol{0.2}};
const VanillaOption kPut{Strike{100}, Maturity{1}, OptionType::Put};

} // namespace

TEST_CASE("LSM American put agrees with the tree within a few standard errors", "[lsm][montecarlo]") {
    const LongstaffSchwartz<GBM> lsm(kPut, LongstaffSchwartzConfig{.paths = 100'000, .steps = 50});
    const Estimate e = lsm.price(kMarket, SeedKey{2026});

    // Reference: docs/model_risk_report.md Appendix C, CRR at n = 20,000.
    constexpr double kTreeReference = 6.090333;
    INFO("LSM price " << e.value << " +/- " << e.std_error << " (bias estimate " << e.discretization << ")");
    CHECK(std::abs(e.value - kTreeReference) < 4.0 * e.std_error);
}

TEST_CASE("LSM price is between the European price and the intrinsic value", "[lsm][montecarlo]") {
    const LongstaffSchwartz<GBM> lsm(kPut, LongstaffSchwartzConfig{.paths = 50'000, .steps = 25});
    const Estimate e = lsm.price(kMarket, SeedKey{7});
    const double european = black_scholes_price(kPut, kMarket);
    const double intrinsic = std::max(0.0, kPut.strike.value - kMarket.spot.value);

    // Early exercise only adds value: American >= European. It is also worth at least what
    // exercising immediately pays.
    CHECK(e.value > european - 4.0 * e.std_error);
    CHECK(e.value >= intrinsic - 4.0 * e.std_error);
}

TEST_CASE("LSM is a pure function of (market, key)", "[lsm][montecarlo][determinism]") {
    const LongstaffSchwartz<GBM> lsm(kPut, LongstaffSchwartzConfig{.paths = 10'000, .steps = 20});
    const Estimate a = lsm.price(kMarket, SeedKey{42});
    const Estimate b = lsm.price(kMarket, SeedKey{42});
    CHECK(a.value == b.value);
    CHECK(a.std_error == b.std_error);
}

TEST_CASE("American call without dividends never exercises early (Merton)", "[lsm][montecarlo]") {
    // Merton's theorem: with q = 0, the American call equals the European call exactly, because
    // early exercise only forfeits the remaining time value. LSM cannot prove this the way the
    // tree can (no dividend means the regression never finds an early exercise worth taking, but
    // it is still a statistical estimate) - the invariant to check is that it does not
    // systematically find spurious early-exercise value, not that it matches to machine precision.
    const VanillaOption call{Strike{100}, Maturity{1}, OptionType::Call};
    const LongstaffSchwartz<GBM> lsm(call, LongstaffSchwartzConfig{.paths = 100'000, .steps = 50});
    const Estimate e = lsm.price(kMarket, SeedKey{11});
    const double european = black_scholes_price(call, kMarket);
    INFO("LSM call " << e.value << " +/- " << e.std_error << ", European " << european);
    CHECK(std::abs(e.value - european) < 4.0 * e.std_error);
}

TEST_CASE("Deep in and out of the money put", "[lsm][montecarlo][edge]") {
    // Deep ITM: exercise now is close to optimal, so LSM should land near the intrinsic value.
    const MarketState deep_itm{Spot{40}, Rate{0.05}, Rate{0.0}, Vol{0.2}};
    const LongstaffSchwartz<GBM> lsm_itm(kPut, LongstaffSchwartzConfig{.paths = 50'000, .steps = 25});
    const Estimate itm = lsm_itm.price(deep_itm, SeedKey{3});
    CHECK(itm.value > 55.0); // intrinsic is 60; deep ITM American puts are worth close to it
    CHECK(itm.value <= 60.0 + 4.0 * itm.std_error);

    // Deep OTM: the option is nearly worthless, but must never be negative or NaN.
    const MarketState deep_otm{Spot{200}, Rate{0.05}, Rate{0.0}, Vol{0.2}};
    const LongstaffSchwartz<GBM> lsm_otm(kPut, LongstaffSchwartzConfig{.paths = 50'000, .steps = 25});
    const Estimate otm = lsm_otm.price(deep_otm, SeedKey{4});
    CHECK(otm.value >= 0.0);
    CHECK(otm.value < 1.0);
    CHECK_FALSE(std::isnan(otm.value));
}

TEST_CASE("Single exercise date reduces to a European option", "[lsm][montecarlo][edge]") {
    // With one exercise date (maturity only), there is no early-exercise decision to fit, so LSM
    // is just a plain Monte Carlo European price.
    const LongstaffSchwartz<GBM> lsm(kPut, LongstaffSchwartzConfig{.paths = 50'000, .steps = 1});
    const Estimate e = lsm.price(kMarket, SeedKey{5});
    const double european = black_scholes_price(kPut, kMarket);
    CHECK(std::abs(e.value - european) < 4.0 * e.std_error);
}

TEST_CASE("LSM is bit-identical for any thread count", "[lsm][montecarlo][determinism]") {
    // The price, its standard error and the bias estimate must be the same bits, not just close:
    // the regression sums and the final average always run on one thread, in path order.
    const auto bits = [](double x) { return std::bit_cast<std::uint64_t>(x); };
    const auto price = [](std::uint64_t paths, std::uint32_t steps, unsigned threads) {
        return LongstaffSchwartz<GBM>(kPut, {.paths = paths, .steps = steps, .threads = threads})
            .price(kMarket, SeedKey{2026});
    };
    for (const auto& [paths, steps] : {std::pair<std::uint64_t, std::uint32_t>{20'000, 50}, {20'001, 7}, {5, 3}, {9'000, 1}}) {
        const Estimate one = price(paths, steps, 1);
        for (unsigned threads : {2u, 3u, 7u, 16u, 64u}) { // odd splits, more than cores, more than paths
            const Estimate many = price(paths, steps, threads);
            CHECK(bits(many.value) == bits(one.value));
            CHECK(bits(many.std_error) == bits(one.std_error));
            CHECK(bits(many.discretization) == bits(one.discretization));
        }
    }
}
