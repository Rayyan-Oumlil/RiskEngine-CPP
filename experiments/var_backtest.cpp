// Report 7.4: backtest of four one-day 99 % VaR methods on the hedged short straddle, every trading
// day from 1992 to today. Each day t the book is set up afresh (short one 30-day ATM straddle,
// delta-hedged; spot normalized to 100, vol = VIX_t, rate = 3-month T-bill as of t), its VaR is
// forecast, and the realized loss is the full revaluation at t + 1 (the NASDAQ move, VIX_{t+1}, the
// rate as of t + 1, one trading day of decay). Writes the daily series (var_backtest.csv) and the
// test statistics per method (var_backtest_summary.csv).
//
// Methods: delta-normal; delta-gamma normal at the day's implied vol; Monte Carlo full revaluation
// (spot only, GBM at the day's implied vol, 20,000 common scenarios); historical full revaluation on
// the preceding 500 days of NASDAQ returns and VIX changes (spot and vol).

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "harness/experiment.hpp"
#include "riskengine/core/rng/random_stream.hpp"
#include "riskengine/io/time_series.hpp"
#include "riskengine/risk/approximations.hpp"
#include "riskengine/risk/backtest.hpp"
#include "riskengine/risk/portfolio.hpp"

using namespace riskengine;

namespace {

constexpr double kLevel = 0.99;
constexpr double kHorizon = 1.0 / 252.0;
constexpr std::size_t kWindow = 500;
constexpr std::size_t kMcScenarios = 20'000;
const Maturity kTenor{30.0 / 365.0};

// Value of the most recent observation on or before `date` (as-of join); the series must start
// before the first date asked for.
class AsOf {
public:
    explicit AsOf(std::vector<Observation> s) : s_(std::move(s)) {}
    double operator()(const std::string& date) {
        while (i_ + 1 < s_.size() && s_[i_ + 1].date <= date) ++i_;
        if (s_[i_].date > date) throw std::runtime_error("as-of series starts after " + date);
        return s_[i_].value;
    }

private:
    std::vector<Observation> s_;
    std::size_t i_ = 0;
};

double loss_after(const Portfolio& book, const MarketState& today, const MarketState& tomorrow) {
    return book.value(today) - book.value(tomorrow, kHorizon);
}

struct Method {
    const char* name;
    std::vector<double> var;
};

} // namespace

int main(int argc, char** argv) {
    return harness::run("var_backtest", argc, argv, [](harness::Experiment& exp) {
        const AlignedPair px = align(read_fred_csv("data/raw/NASDAQCOM.csv"), read_fred_csv("data/raw/VIXCLS.csv"));
        AsOf rate(read_fred_csv("data/raw/DGS3MO.csv"));
        const std::size_t first = kWindow, last = px.dates.size() - 1; // forecast days first .. last - 1
        exp.param("position", "short 1 ATM call + put, T = 30/365, delta-hedged, spot normalized to 100");
        exp.param("vol_and_rate", "VIX_t / 100 and DGS3MO as of t / 100");
        exp.param("period", px.dates[first] + " .. " + px.dates[last - 1] + " (forecast days)");
        exp.param("var_level", kLevel);
        exp.param("horizon_years", kHorizon);
        exp.param("historical_window_days", kWindow);
        exp.param("mc_scenarios", kMcScenarios);
        exp.param("mc_seed", "SeedKey{2026}, the same scenarios every day (common random numbers)");

        std::vector<double> z(kMcScenarios);
        RandomStream rng(SeedKey{2026}, 0);
        for (double& v : z) v = rng.normal();

        std::vector<Method> methods{{"delta_normal", {}}, {"delta_gamma_normal", {}}, {"mc_full_revaluation", {}},
                                    {"historical_full_revaluation", {}}};
        std::vector<double> realized;
        std::vector<std::string> dates;
        auto daily = exp.csv({"date", "vix", "realized_loss", "var_delta_normal", "var_delta_gamma_normal",
                              "var_mc_full_revaluation", "var_historical_full_revaluation"});
        std::vector<double> scratch;
        for (std::size_t t = first; t < last; ++t) {
            const double r_t = rate(px.dates[t]) / 100.0;
            const MarketState today{Spot{100}, Rate{r_t}, Rate{0.0}, Vol{px.b[t] / 100.0}};
            const Portfolio book = delta_hedged_short_straddle(Strike{100}, kTenor, today);
            const Greeks g = book.greeks(today);
            const double sd = today.vol.value * std::sqrt(kHorizon);

            // Monte Carlo full revaluation at the day's implied vol.
            scratch.clear();
            const double drift = (r_t - 0.5 * today.vol.value * today.vol.value) * kHorizon;
            for (double zi : z) {
                MarketState s = today;
                s.spot.value = 100.0 * std::exp(drift + sd * zi);
                scratch.push_back(loss_after(book, today, s));
            }
            const double var_mc = empirical_var_es(scratch, kLevel).var;
            // Historical full revaluation: the preceding kWindow days' spot returns and VIX changes.
            scratch.clear();
            for (std::size_t i = t - kWindow + 1; i <= t; ++i) {
                MarketState s = today;
                s.spot.value = 100.0 * px.a[i] / px.a[i - 1];
                s.vol.value = std::max(0.01, today.vol.value + (px.b[i] - px.b[i - 1]) / 100.0);
                scratch.push_back(loss_after(book, today, s));
            }
            const double var_hist = empirical_var_es(scratch, kLevel).var;
            const double var_dn = delta_normal(g, 100, sd, kLevel).var;
            const double var_dg = delta_gamma_normal(g, 100, sd, kHorizon, kLevel).var;

            // What happened: tomorrow's move, vol and rate.
            const MarketState next{Spot{100.0 * px.a[t + 1] / px.a[t]}, Rate{rate(px.dates[t + 1]) / 100.0}, Rate{0.0},
                                   Vol{px.b[t + 1] / 100.0}};
            const double loss = loss_after(book, today, next);

            methods[0].var.push_back(var_dn);
            methods[1].var.push_back(var_dg);
            methods[2].var.push_back(var_mc);
            methods[3].var.push_back(var_hist);
            realized.push_back(loss);
            dates.push_back(px.dates[t]);
            daily.row(px.dates[t], px.b[t], loss, var_dn, var_dg, var_mc, var_hist);
        }

        harness::Experiment summary("var_backtest_summary", exp.out_dir());
        summary.param("source", "var_backtest (same run)");
        summary.param("period", px.dates[first] + " .. " + px.dates[last - 1]);
        summary.param("crisis_windows", "2008-09-01 .. 2009-03-31 and 2020-02-15 .. 2020-04-30");
        auto table = summary.csv({"method", "days", "exceptions", "exception_rate", "kupiec_lr", "kupiec_p",
                                  "christoffersen_lr", "christoffersen_p", "conditional_coverage_p",
                                  "red_zone_windows", "windows", "exceptions_2008", "exceptions_2020"});
        for (const Method& m : methods) {
            std::vector<char> hits(realized.size());
            std::size_t x = 0, crisis_2008 = 0, crisis_2020 = 0;
            for (std::size_t i = 0; i < realized.size(); ++i) {
                hits[i] = realized[i] > m.var[i];
                x += hits[i] ? 1 : 0;
                if (hits[i] && dates[i] >= "2008-09-01" && dates[i] <= "2009-03-31") ++crisis_2008;
                if (hits[i] && dates[i] >= "2020-02-15" && dates[i] <= "2020-04-30") ++crisis_2020;
            }
            // Non-overlapping 250-day windows, as a supervisor would count them year by year.
            std::size_t red = 0, windows = 0;
            for (std::size_t start = 0; start + 250 <= hits.size(); start += 250, ++windows) {
                std::size_t in_window = 0;
                for (std::size_t i = start; i < start + 250; ++i) in_window += hits[i] ? 1 : 0;
                red += basel_zone(in_window) == BaselZone::Red ? 1 : 0;
            }
            std::unique_ptr<bool[]> buf(new bool[hits.size()]);
            for (std::size_t i = 0; i < hits.size(); ++i) buf[i] = hits[i] != 0;
            const std::span<const bool> seq(buf.get(), hits.size());
            const TestResult pof = kupiec_pof(hits.size(), x, 1.0 - kLevel);
            const TestResult ind = christoffersen_independence(seq);
            const TestResult cc = christoffersen_conditional_coverage(seq, 1.0 - kLevel);
            table.row(m.name, hits.size(), x, static_cast<double>(x) / static_cast<double>(hits.size()), pof.statistic,
                      pof.p_value, ind.statistic, ind.p_value, cc.p_value, red, windows, crisis_2008, crisis_2020);
        }
        summary.write_metadata();
    });
}
