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

TEST_CASE("Terminal-spot control: unbiased, strong in the money, weak far out of it", "[mc][control]") {
    const MonteCarloConfig cfg{.paths = 200'000, .threads = 4};
    auto se_ratio = [&](double k) {
        const VanillaOption o{Strike{k}, Maturity{kMaturity}, OptionType::Call};
        const Estimate plain = MonteCarlo<GBM, VanillaPayoff>({k, OptionType::Call}, Maturity{kMaturity}, cfg).price(kMarket, SeedKey{10});
        const Estimate ctl = MonteCarlo<GBM, VanillaPayoff, TerminalSpotControl>({k, OptionType::Call}, Maturity{kMaturity}, cfg)
                                 .price(kMarket, SeedKey{10});
        check_within_4se(ctl, black_scholes_price(o, kMarket));
        return ctl.std_error / plain.std_error;
    };
    CHECK(se_ratio(70.0) < 0.3);  // deep in the money: the call is almost S_T - K
    CHECK(se_ratio(160.0) > 0.7); // far out of the money: correlation with S_T collapses
}

TEST_CASE("Geometric control on the arithmetic Asian: same price, far tighter", "[mc][control]") {
    const MonteCarloConfig cfg{.paths = 100'000, .steps = 12, .threads = 4};
    const Estimate plain = MonteCarlo<GBM, ArithmeticAsianPayoff>({100.0, OptionType::Call}, Maturity{kMaturity}, cfg)
                               .price(kMarket, SeedKey{11});
    const Estimate ctl = MonteCarlo<GBM, ArithmeticAsianPayoff, GeometricAsianControl>(
                             {100.0, OptionType::Call}, Maturity{kMaturity}, cfg, {100.0, OptionType::Call})
                             .price(kMarket, SeedKey{12});
    const double combined_se = std::hypot(plain.std_error, ctl.std_error);
    CHECK(std::abs(plain.value - ctl.value) < 4.0 * combined_se);
    CHECK(ctl.std_error < plain.std_error / 10.0);
}

TEST_CASE("A payoff used as its own control is priced exactly", "[mc][control]") {
    // Y = X: beta = 1 and every sample equals mu, whatever the paths.
    const MonteCarloConfig cfg{.paths = 20'000, .steps = 12, .threads = 2};
    const Estimate e = MonteCarlo<GBM, GeometricAsianPayoff, GeometricAsianControl>(
                           {100.0, OptionType::Put}, Maturity{kMaturity}, cfg, {100.0, OptionType::Put})
                           .price(kMarket, SeedKey{13});
    const double exact = geometric_asian_price(OptionType::Put, Strike{100}, Maturity{kMaturity}, 12, kMarket);
    CHECK(std::abs(e.value - exact) <= 1e-12 * exact);
    CHECK(e.std_error <= 1e-12 * exact);
}

TEST_CASE("Control and antithetic variates combine, deterministically", "[mc][control]") {
    MonteCarloConfig cfg{.paths = 100'000, .threads = 1, .antithetic = true};
    const MonteCarlo<GBM, VanillaPayoff, TerminalSpotControl> pricer({110.0, OptionType::Put}, Maturity{kMaturity}, cfg);
    const Estimate a = pricer.price(kMarket, SeedKey{14});
    check_within_4se(a, black_scholes_price(VanillaOption{Strike{110}, Maturity{kMaturity}, OptionType::Put}, kMarket));
    cfg.threads = 8;
    const Estimate b = MonteCarlo<GBM, VanillaPayoff, TerminalSpotControl>({110.0, OptionType::Put}, Maturity{kMaturity}, cfg)
                           .price(kMarket, SeedKey{14});
    CHECK(a.value == b.value);
    CHECK(a.std_error == b.std_error);
}

TEST_CASE("Randomized QMC is unbiased within 4 SE, with and without a bridge", "[mc][qmc]") {
    MonteCarloConfig cfg{.paths = 4096, .threads = 4, .sampling = Sampling::RandomizedQmc, .replications = 16};
    const VanillaOption call{Strike{100}, Maturity{kMaturity}, OptionType::Call};
    check_within_4se(MonteCarlo<GBM, VanillaPayoff>({100.0, OptionType::Call}, Maturity{kMaturity}, cfg).price(kMarket, SeedKey{20}),
                     black_scholes_price(call, kMarket));
    check_within_4se(MonteCarlo<GBM, DigitalPayoff>({110.0, OptionType::Put}, Maturity{kMaturity}, cfg).price(kMarket, SeedKey{21}),
                     digital_price(VanillaOption{Strike{110}, Maturity{kMaturity}, OptionType::Put}, kMarket));
    cfg.steps = 12;
    for (bool bridge : {false, true}) {
        INFO("bridge " << bridge);
        cfg.brownian_bridge = bridge;
        check_within_4se(MonteCarlo<GBM, GeometricAsianPayoff>({100.0, OptionType::Call}, Maturity{kMaturity}, cfg)
                             .price(kMarket, SeedKey{22}),
                         geometric_asian_price(OptionType::Call, Strike{100}, Maturity{kMaturity}, 12, kMarket));
    }
    const Estimate e = MonteCarlo<GBM, VanillaPayoff>({100.0, OptionType::Call}, Maturity{kMaturity}, cfg).price(kMarket, SeedKey{20});
    CHECK(e.samples == 16u * 4096u);
}

TEST_CASE("Randomized QMC beats pseudo-random at equal cost on a smooth problem", "[mc][qmc]") {
    // 16 x 4096 scrambled Sobol points against 65,536 pseudo-random paths.
    const MonteCarloConfig qmc{.paths = 4096, .threads = 4, .sampling = Sampling::RandomizedQmc, .replications = 16};
    const MonteCarloConfig prng{.paths = 16 * 4096, .threads = 4};
    const Estimate q = MonteCarlo<GBM, VanillaPayoff>({100.0, OptionType::Call}, Maturity{kMaturity}, qmc).price(kMarket, SeedKey{23});
    const Estimate p = MonteCarlo<GBM, VanillaPayoff>({100.0, OptionType::Call}, Maturity{kMaturity}, prng).price(kMarket, SeedKey{23});
    CHECK(q.std_error < p.std_error / 10.0);
}

TEST_CASE("The Brownian bridge sharpens QMC on a path payoff and is harmless with pseudo-random", "[mc][qmc][bridge]") {
    MonteCarloConfig cfg{.paths = 4096, .steps = 12, .threads = 4, .sampling = Sampling::RandomizedQmc, .replications = 16};
    const ArithmeticAsianPayoff asian{100.0, OptionType::Call};
    const Estimate plain = MonteCarlo<GBM, ArithmeticAsianPayoff>(asian, Maturity{kMaturity}, cfg).price(kMarket, SeedKey{24});
    cfg.brownian_bridge = true;
    const Estimate bridged = MonteCarlo<GBM, ArithmeticAsianPayoff>(asian, Maturity{kMaturity}, cfg).price(kMarket, SeedKey{24});
    CHECK(bridged.std_error < plain.std_error / 2.0);

    // Pseudo-random with a bridge is the same distribution: still unbiased.
    const MonteCarloConfig prng{.paths = 100'000, .steps = 12, .threads = 4, .brownian_bridge = true};
    check_within_4se(MonteCarlo<GBM, GeometricAsianPayoff>({100.0, OptionType::Put}, Maturity{kMaturity}, prng).price(kMarket, SeedKey{25}),
                     geometric_asian_price(OptionType::Put, Strike{100}, Maturity{kMaturity}, 12, kMarket));
}

TEST_CASE("Randomized QMC is thread-invariant and composes with antithetic and control variates", "[mc][qmc]") {
    MonteCarloConfig cfg{.paths = 2048, .steps = 12, .threads = 1, .antithetic = true,
                         .sampling = Sampling::RandomizedQmc, .replications = 8, .brownian_bridge = true};
    const ArithmeticAsianPayoff asian{100.0, OptionType::Call};
    const GeometricAsianControl control{100.0, OptionType::Call};
    const Estimate one = MonteCarlo<GBM, ArithmeticAsianPayoff, GeometricAsianControl>(asian, Maturity{kMaturity}, cfg, control)
                             .price(kMarket, SeedKey{26});
    cfg.threads = 8;
    const Estimate eight = MonteCarlo<GBM, ArithmeticAsianPayoff, GeometricAsianControl>(asian, Maturity{kMaturity}, cfg, control)
                               .price(kMarket, SeedKey{26});
    CHECK(one.value == eight.value);
    CHECK(one.std_error == eight.std_error);
    // Unbiased against a large pseudo-random reference with the same control.
    const MonteCarloConfig ref_cfg{.paths = 400'000, .steps = 12, .threads = 4};
    const Estimate ref = MonteCarlo<GBM, ArithmeticAsianPayoff, GeometricAsianControl>(asian, Maturity{kMaturity}, ref_cfg, control)
                             .price(kMarket, SeedKey{27});
    CHECK(std::abs(one.value - ref.value) < 4.0 * std::hypot(one.std_error, ref.std_error));
}
