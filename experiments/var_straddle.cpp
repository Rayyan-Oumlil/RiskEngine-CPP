// Report 7.1-7.3: the risk of a delta-hedged short ATM straddle over one trading day, by seven
// methods. Delta-normal sees nothing (the book has no delta); delta-gamma sees the gamma but
// forces a normal shape onto a skewed quadratic loss; full revaluation (Monte Carlo or
// historical) reprices every scenario exactly. The spot-and-vol variants add the vega risk that a
// spot-only VaR ignores. VaR is at 99 %, ES at 97.5 % (the FRTB pair); every empirical measure
// carries a 95 % bootstrap interval.
//
// Reads data/raw/NASDAQCOM.csv and data/raw/VIXCLS.csv (frozen, see tools/fetch_data.py): the last
// 2,520 aligned trading days give the historical scenarios and calibrate the Monte Carlo vol shocks.
// A last variant rescales the historical returns to the 20 % vol of the other methods, to separate
// fat tails from a different vol level.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "harness/experiment.hpp"
#include "riskengine/core/rng/random_stream.hpp"
#include "riskengine/core/stats/welford.hpp"
#include "riskengine/io/time_series.hpp"
#include "riskengine/risk/approximations.hpp"
#include "riskengine/risk/portfolio.hpp"

using namespace riskengine;

namespace {

constexpr double kVarLevel = 0.99, kEsLevel = 0.975;
constexpr double kHorizon = 1.0 / 252.0; // one trading day, for diffusion and time decay alike
constexpr std::size_t kHistoryDays = 2520;
constexpr std::uint64_t kScenarios = 1'000'000;
constexpr std::uint32_t kBootstrap = 200;

struct Scenario {
    double log_return; // of the underlying over the horizon
    double vol_change; // additive change of the implied vol (decimal)
};

// Loss of the book over the horizon in each scenario, by full revaluation.
std::vector<double> revalue(const Portfolio& book, const MarketState& m, const std::vector<Scenario>& scenarios) {
    const double v0 = book.value(m);
    std::vector<double> losses;
    losses.reserve(scenarios.size());
    for (const Scenario& s : scenarios) {
        MarketState shocked = m;
        shocked.spot.value = m.spot.value * std::exp(s.log_return);
        shocked.vol.value = std::max(0.01, m.vol.value + s.vol_change);
        losses.push_back(v0 - book.value(shocked, kHorizon));
    }
    return losses;
}

} // namespace

int main(int argc, char** argv) {
    return harness::run("var_straddle", argc, argv, [](harness::Experiment& exp) {
        const MarketState market{Spot{100}, Rate{0.05}, Rate{0.0}, Vol{0.2}};
        const Portfolio book = delta_hedged_short_straddle(Strike{100}, Maturity{30.0 / 365}, market);
        const Greeks g = book.greeks(market);
        const double daily_sd = market.vol.value * std::sqrt(kHorizon);

        // History: aligned NASDAQ closes and VIX levels, last kHistoryDays + 1 days.
        const AlignedPair hist = align(read_fred_csv("data/raw/NASDAQCOM.csv"), read_fred_csv("data/raw/VIXCLS.csv"));
        if (hist.dates.size() <= kHistoryDays) throw std::runtime_error("not enough aligned history");
        const std::size_t first = hist.dates.size() - kHistoryDays - 1;
        std::vector<Scenario> history;
        for (std::size_t i = first + 1; i < hist.dates.size(); ++i)
            history.push_back({std::log(hist.a[i] / hist.a[i - 1]), (hist.b[i] - hist.b[i - 1]) / 100.0});
        // Calibration of the Monte Carlo vol shock: sd of daily vol changes and their correlation
        // with the daily log return, over the same window.
        Welford r_stats, v_stats;
        for (const Scenario& s : history) {
            r_stats.add(s.log_return);
            v_stats.add(s.vol_change);
        }
        double cov = 0.0;
        for (const Scenario& s : history) cov += (s.log_return - r_stats.mean()) * (s.vol_change - v_stats.mean());
        cov /= static_cast<double>(history.size() - 1);
        const double vol_shock_sd = std::sqrt(v_stats.variance());
        const double rho = cov / std::sqrt(r_stats.variance() * v_stats.variance());

        exp.param("position", "short 1 call + 1 put, K = 100, T = 30/365, delta-hedged with the underlying");
        exp.param("market", "S = 100, r = 0.05, q = 0, sigma = 0.2");
        exp.param("horizon_years", kHorizon);
        exp.param("var_level", kVarLevel);
        exp.param("es_level", kEsLevel);
        exp.param("book_delta", g.delta);
        exp.param("book_gamma", g.gamma);
        exp.param("book_vega", g.vega);
        exp.param("book_theta", g.theta);
        exp.param("mc_scenarios", kScenarios);
        exp.param("mc_seed", "SeedKey{2026}: spot normals from stream 0, vol normals from stream 1");
        exp.param("history", "NASDAQCOM and VIXCLS, aligned, " + hist.dates[first] + " .. " + hist.dates.back());
        exp.param("history_days", kHistoryDays);
        const double realized_vol = std::sqrt(r_stats.variance() / kHorizon);
        exp.param("history_realized_vol", realized_vol);
        exp.param("vol_shock_sd", vol_shock_sd);
        exp.param("spot_vol_correlation", rho);
        exp.param("bootstrap", "200 resamples, 95 % percentile interval, SeedKey{7}");

        auto csv = exp.csv({"method", "risk_factors", "var99", "var99_low", "var99_high", "es975", "es975_low",
                            "es975_high"});
        const double nan = std::numeric_limits<double>::quiet_NaN();
        auto analytic = [&](const char* method, double var, double es) {
            csv.row(method, "spot", var, nan, nan, es, nan, nan);
        };
        auto empirical = [&](const char* method, const char* factors, const std::vector<double>& losses) {
            const RiskMeasuresWithCI v = bootstrap_var_es(losses, kVarLevel, kBootstrap, SeedKey{7});
            const RiskMeasuresWithCI e = bootstrap_var_es(losses, kEsLevel, kBootstrap, SeedKey{8});
            csv.row(method, factors, v.estimate.var, v.var.low, v.var.high, e.estimate.es, e.es.low, e.es.high);
        };

        analytic("delta_normal", delta_normal(g, 100, daily_sd, kVarLevel).var, delta_normal(g, 100, daily_sd, kEsLevel).es);
        analytic("delta_gamma_normal", delta_gamma_normal(g, 100, daily_sd, kHorizon, kVarLevel).var,
                 delta_gamma_normal(g, 100, daily_sd, kHorizon, kEsLevel).es);
        analytic("delta_gamma_cornish_fisher", delta_gamma_cornish_fisher_var(g, 100, daily_sd, kHorizon, kVarLevel), nan);

        // Monte Carlo: GBM log-returns over the horizon, then correlated vol shocks.
        RandomStream spot_rng(SeedKey{2026, 0}, 0), vol_rng(SeedKey{2026, 1}, 0);
        const double drift = (market.rate.value - 0.5 * market.vol.value * market.vol.value) * kHorizon;
        std::vector<Scenario> mc_spot, mc_joint;
        mc_spot.reserve(kScenarios);
        mc_joint.reserve(kScenarios);
        for (std::uint64_t i = 0; i < kScenarios; ++i) {
            const double z1 = spot_rng.normal(), z2 = vol_rng.normal();
            const double x = drift + daily_sd * z1;
            mc_spot.push_back({x, 0.0});
            mc_joint.push_back({x, vol_shock_sd * (rho * z1 + std::sqrt(1.0 - rho * rho) * z2)});
        }
        empirical("full_revaluation_mc", "spot", revalue(book, market, mc_spot));
        empirical("full_revaluation_mc", "spot+vol", revalue(book, market, mc_joint));

        std::vector<Scenario> hist_spot = history;
        for (Scenario& s : hist_spot) s.vol_change = 0.0;
        empirical("historical", "spot", revalue(book, market, hist_spot));
        empirical("historical", "spot+vol", revalue(book, market, history));
        // The same history with its returns rescaled to the 20 % vol of the other methods: what is
        // left of the gap to the Monte Carlo VaR is the shape of the tails, not the level of vol.
        std::vector<Scenario> rescaled = hist_spot;
        for (Scenario& s : rescaled) s.log_return *= market.vol.value / realized_vol;
        empirical("historical_rescaled_to_20pct", "spot", revalue(book, market, rescaled));
    });
}
