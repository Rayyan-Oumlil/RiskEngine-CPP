#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <cmath>

#include "riskengine/methods/analytic/digital.hpp"
#include "riskengine/methods/analytic/geometric_asian.hpp"
#include "riskengine/methods/montecarlo/engine.hpp"
#include "riskengine/models/gbm.hpp"
#include "riskengine/payoffs/asian.hpp"
#include "riskengine/payoffs/digital.hpp"

using namespace riskengine;

static_assert(PathModel<GBM>);
static_assert(StochasticPricer<MonteCarlo<GBM, VanillaPayoff>>);
static_assert(StochasticPricer<MonteCarlo<GBM, ArithmeticAsianPayoff>>);
static_assert(TerminalPayoff<DigitalPayoff> && !PathPayoff<DigitalPayoff>);
static_assert(PathPayoff<GeometricAsianPayoff> && !TerminalPayoff<GeometricAsianPayoff>);

namespace {

const MarketState kMarket{Spot{100}, Rate{0.05}, Rate{0.02}, Vol{0.25}};
constexpr double kMaturity = 1.0;

// The Phase 4 acceptance test: |MC - exact| < 4 standard errors, with a fixed seed so the test is
// deterministic (at 4 SE an unseeded test would fail ~6e-5 of the time).
void check_within_4se(const Estimate& e, double exact) {
    INFO("MC " << e.value << " +- " << e.std_error << ", exact " << exact);
    CHECK(e.std_error > 0.0);
    CHECK(std::abs(e.value - exact) < 4.0 * e.std_error);
}

} // namespace

TEST_CASE("European payoffs agree with the closed forms within 4 SE", "[mc]") {
    const auto steps = GENERATE(1u, 12u); // GBM steps are exact: no bias from time stepping
    for (double k : {80.0, 100.0, 120.0})
        for (OptionType type : {OptionType::Call, OptionType::Put}) {
            INFO("K=" << k << " steps=" << steps << (type == OptionType::Call ? " call" : " put"));
            const VanillaOption o{Strike{k}, Maturity{kMaturity}, type};
            const MonteCarloConfig cfg{.paths = 200'000, .steps = steps, .threads = 4};
            const MonteCarlo<GBM, VanillaPayoff> vanilla({k, type}, Maturity{kMaturity}, cfg);
            check_within_4se(vanilla.price(kMarket, SeedKey{1}), black_scholes_price(o, kMarket));
            const MonteCarlo<GBM, DigitalPayoff> digital({k, type}, Maturity{kMaturity}, cfg);
            check_within_4se(digital.price(kMarket, SeedKey{2}), digital_price(o, kMarket));
        }
}

TEST_CASE("Straddle and geometric Asian agree with the closed forms within 4 SE", "[mc]") {
    const MonteCarloConfig cfg{.paths = 200'000, .steps = 12, .threads = 4};
    const MonteCarlo<GBM, StraddlePayoff> straddle({100.0}, Maturity{kMaturity}, cfg);
    const double straddle_exact = black_scholes_price(VanillaOption{Strike{100}, Maturity{kMaturity}, OptionType::Call}, kMarket) +
                                  black_scholes_price(VanillaOption{Strike{100}, Maturity{kMaturity}, OptionType::Put}, kMarket);
    check_within_4se(straddle.price(kMarket, SeedKey{3}), straddle_exact);

    for (OptionType type : {OptionType::Call, OptionType::Put}) {
        const MonteCarlo<GBM, GeometricAsianPayoff> asian({100.0, type}, Maturity{kMaturity}, cfg);
        check_within_4se(asian.price(kMarket, SeedKey{4}),
                         geometric_asian_price(type, Strike{100}, Maturity{kMaturity}, cfg.steps, kMarket));
    }
}

TEST_CASE("Antithetic variates: unbiased, and tighter on a monotone payoff", "[mc][antithetic]") {
    const VanillaOption o{Strike{100}, Maturity{kMaturity}, OptionType::Call};
    const MonteCarloConfig plain{.paths = 200'000, .threads = 4};
    MonteCarloConfig anti = plain;
    anti.antithetic = true;
    const Estimate p = MonteCarlo<GBM, VanillaPayoff>({100.0, OptionType::Call}, Maturity{kMaturity}, plain).price(kMarket, SeedKey{5});
    const Estimate a = MonteCarlo<GBM, VanillaPayoff>({100.0, OptionType::Call}, Maturity{kMaturity}, anti).price(kMarket, SeedKey{5});
    check_within_4se(a, black_scholes_price(o, kMarket));
    CHECK(a.samples == plain.paths / 2);
    // Same number of simulated paths: the pairs are negatively correlated for a monotone payoff.
    CHECK(a.std_error < 0.8 * p.std_error);
}

TEST_CASE("The price is a pure function of (market, key, paths, steps, blocks)", "[mc][determinism]") {
    const MonteCarloConfig one{.paths = 50'000, .steps = 4, .threads = 1};
    MonteCarloConfig many = one;
    many.threads = 8;
    const MonteCarlo<GBM, ArithmeticAsianPayoff> p1({100.0, OptionType::Call}, Maturity{kMaturity}, one);
    const MonteCarlo<GBM, ArithmeticAsianPayoff> p8({100.0, OptionType::Call}, Maturity{kMaturity}, many);
    const Estimate a = p1.price(kMarket, SeedKey{6});
    CHECK(a.value == p8.price(kMarket, SeedKey{6}).value);
    CHECK(a.value == p1.price(kMarket, SeedKey{6}).value);
    CHECK(a.value != p1.price(kMarket, SeedKey{7}).value);
}

TEST_CASE("Expired option prices at its intrinsic value with zero error", "[mc]") {
    const MonteCarlo<GBM, VanillaPayoff> expired({90.0, OptionType::Call}, Maturity{0.0}, MonteCarloConfig{.paths = 1000});
    const Estimate e = expired.price(kMarket, SeedKey{8});
    CHECK(e.value == 10.0);
    CHECK(e.std_error == 0.0);
}
