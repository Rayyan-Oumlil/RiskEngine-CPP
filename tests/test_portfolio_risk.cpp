#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

#include "riskengine/core/rng/random_stream.hpp"
#include "riskengine/risk/approximations.hpp"
#include "riskengine/risk/portfolio.hpp"

using namespace riskengine;

namespace {

const MarketState kMarket{Spot{100}, Rate{0.05}, Rate{0.0}, Vol{0.2}};
constexpr double kDay = 1.0 / 252.0;

} // namespace

TEST_CASE("The hedged short straddle is delta-neutral, short gamma, long theta", "[portfolio]") {
    const Portfolio p = delta_hedged_short_straddle(Strike{100}, Maturity{30.0 / 365}, kMarket);
    const Greeks g = p.greeks(kMarket);
    CHECK(std::abs(g.delta) < 1e-15);
    CHECK(g.gamma < 0.0);
    CHECK(g.vega < 0.0);
    CHECK(g.theta > 0.0);
    // Aging: value after `elapsed` equals the closed forms at the shorter maturity.
    const double elapsed = 5.0 / 365;
    const double expected = p.stock * 100 - black_scholes_price({Strike{100}, Maturity{25.0 / 365}, OptionType::Call}, kMarket) -
                            black_scholes_price({Strike{100}, Maturity{25.0 / 365}, OptionType::Put}, kMarket);
    CHECK(std::abs(p.value(kMarket, elapsed) - expected) < 1e-13);
    CHECK(p.value(kMarket, 1.0) == p.stock * 100.0); // every option expired at the money: worth 0
}

TEST_CASE("Delta-gamma loss moments are exact for a normal risk factor", "[portfolio][approximations]") {
    // Against 2 x 10^6 sampled losses of the same quadratic, with a delta so that every term matters.
    const Greeks g{0.4, -0.3, 0.0, 3.0, 0.0}; // gamma large enough that the delta-gamma cross terms matter
    const double spot = 100, s = 0.015, dt = kDay;
    const LossMoments m = delta_gamma_moments(g, spot, s, dt);
    RandomStream rng(SeedKey{41}, 0);
    constexpr int n = 2'000'000;
    std::vector<double> l(n);
    double mean = 0.0;
    for (double& v : l) {
        const double x = s * rng.normal();
        v = -(g.delta * spot * x + 0.5 * g.gamma * spot * spot * x * x + g.theta * dt);
        mean += v / n;
    }
    double c2 = 0, c3 = 0, c4 = 0;
    for (double v : l) {
        const double d = v - mean;
        c2 += d * d / n;
        c3 += d * d * d / n;
        c4 += d * d * d * d / n;
    }
    CHECK(std::abs(mean - m.mean) < 4.0 * m.sd / std::sqrt(double(n)));
    CHECK(std::abs(std::sqrt(c2) / m.sd - 1.0) < 0.005);
    CHECK(std::abs(c3 / std::pow(c2, 1.5) - m.skew) < 0.05);
    CHECK(std::abs(c4 / (c2 * c2) - 3.0 - m.excess_kurtosis) < 0.3);
    // A pure quadratic is a scaled chi-square(1): skew sqrt(8), excess kurtosis 12.
    const LossMoments q = delta_gamma_moments(Greeks{0, -0.05, 0, 0, 0}, spot, s, dt);
    CHECK(std::abs(q.skew - std::sqrt(8.0)) < 1e-12);
    CHECK(std::abs(q.excess_kurtosis - 12.0) < 1e-12);
}

TEST_CASE("Delta-normal VaR of a delta-hedged book is zero", "[portfolio][approximations]") {
    const Portfolio p = delta_hedged_short_straddle(Strike{100}, Maturity{30.0 / 365}, kMarket);
    const RiskMeasures r = delta_normal(p.greeks(kMarket), 100, 0.2 * std::sqrt(kDay), 0.99);
    CHECK(std::abs(r.var) < 1e-12);
    CHECK(std::abs(r.es) < 1e-12);
    // Delta-gamma sees the short gamma: a strictly positive loss quantile.
    CHECK(delta_gamma_normal(p.greeks(kMarket), 100, 0.2 * std::sqrt(kDay), kDay, 0.99).var > 0.1);
}

TEST_CASE("Bootstrap intervals are reproducible and cover the true normal VaR", "[var][bootstrap]") {
    std::vector<double> losses(20'000);
    RandomStream rng(SeedKey{42}, 0);
    for (double& l : losses) l = rng.normal();
    const RiskMeasuresWithCI a = bootstrap_var_es(losses, 0.99, 200, SeedKey{43});
    const RiskMeasuresWithCI b = bootstrap_var_es(losses, 0.99, 200, SeedKey{43});
    CHECK(a.var.low == b.var.low);
    CHECK(a.var.low < a.estimate.var);
    CHECK(a.estimate.var < a.var.high);
    const RiskMeasures truth = normal_var_es(0, 1, 0.99);
    CHECK(a.var.low < truth.var);
    CHECK(truth.var < a.var.high);
    CHECK(a.es.low < truth.es);
    CHECK(truth.es < a.es.high);
}
