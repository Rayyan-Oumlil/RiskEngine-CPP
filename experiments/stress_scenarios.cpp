// Report 7.5: historical stress scenarios (docs/riskengine_research.md 9.2) applied to today's
// delta-hedged short straddle (the book of var_straddle: S = K = 100, 30 days, sigma = 20 %,
// r = 5 %). Each scenario is a joint close-to-close move of the NASDAQ Composite, of implied vol
// (VIX; VXO for 1987, before the VIX existed) and of the 3-month T-bill yield, with the time decay
// of its trading days. The loss is also split into a spot-only and a vol-only part: how much of a
// crisis loss a spot-only risk measure could ever see.

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "harness/experiment.hpp"
#include "riskengine/io/time_series.hpp"
#include "riskengine/risk/portfolio.hpp"

using namespace riskengine;

namespace {

std::map<std::string, double> by_date(const std::vector<Observation>& s) {
    std::map<std::string, double> m;
    for (const Observation& o : s) m[o.date] = o.value;
    return m;
}

double at(const std::map<std::string, double>& m, const std::string& date, const char* name) {
    const auto it = m.find(date);
    if (it == m.end()) throw std::runtime_error(std::string(name) + " has no value on " + date);
    return it->second;
}

// Last value on or before `date`.
double as_of(const std::map<std::string, double>& m, const std::string& date) {
    auto it = m.upper_bound(date);
    if (it == m.begin()) throw std::runtime_error("no rate on or before " + date);
    return std::prev(it)->second;
}

} // namespace

int main(int argc, char** argv) {
    return harness::run("stress_scenarios", argc, argv, [](harness::Experiment& exp) {
        const auto nasdaq = read_fred_csv("data/raw/NASDAQCOM.csv");
        const auto nq = by_date(nasdaq), vix = by_date(read_fred_csv("data/raw/VIXCLS.csv")),
                   vxo = by_date(read_fred_csv("data/raw/VXOCLS.csv")), rates = by_date(read_fred_csv("data/raw/DGS3MO.csv"));
        const MarketState market{Spot{100}, Rate{0.05}, Rate{0.0}, Vol{0.2}};
        const Portfolio book = delta_hedged_short_straddle(Strike{100}, Maturity{30.0 / 365}, market);
        const double v0 = book.value(market);

        struct Scenario {
            const char* name;
            const char* start;
            const char* end;
            bool use_vxo;
        };
        const Scenario scenarios[] = {
            {"black_monday_1987", "1987-10-16", "1987-10-19", true},
            {"tarp_vote_2008", "2008-09-26", "2008-09-29", false},
            {"october_2008_week", "2008-10-03", "2008-10-10", false},
            {"volmageddon_2018", "2018-02-02", "2018-02-05", false},
            {"covid_2020_day", "2020-03-13", "2020-03-16", false},
            {"covid_2020_week", "2020-03-06", "2020-03-16", false},
        };
        exp.param("book", "short 1 ATM straddle, K = 100, T = 30/365, delta-hedged; S = 100, sigma = 0.2, r = 0.05");
        exp.param("shocks", "close-to-close NASDAQCOM log-return, VIX (VXO for 1987) change / 100, DGS3MO change / 100");
        exp.param("time_decay", "trading days in the window / 252");

        auto csv = exp.csv({"scenario", "start", "end", "trading_days", "spot_return", "vol_change", "rate_change",
                            "loss", "loss_spot_only", "loss_vol_only"});
        for (const Scenario& s : scenarios) {
            const auto first = std::find_if(nasdaq.begin(), nasdaq.end(), [&](const Observation& o) { return o.date == s.start; });
            const auto last = std::find_if(nasdaq.begin(), nasdaq.end(), [&](const Observation& o) { return o.date == s.end; });
            if (first == nasdaq.end() || last == nasdaq.end()) throw std::runtime_error(std::string("missing dates for ") + s.name);
            const auto days = static_cast<int>(last - first);
            const auto& vol = s.use_vxo ? vxo : vix;
            const double spot_return = std::log(at(nq, s.end, "NASDAQCOM") / at(nq, s.start, "NASDAQCOM"));
            const double vol_change = (at(vol, s.end, "vol") - at(vol, s.start, "vol")) / 100.0;
            const double rate_change = (as_of(rates, s.end) - as_of(rates, s.start)) / 100.0;
            const double elapsed = days / 252.0;

            auto shocked = [&](bool spot, bool vol_and_rate) {
                MarketState m = market;
                if (spot) m.spot.value = 100.0 * std::exp(spot_return);
                if (vol_and_rate) {
                    m.vol.value = std::max(0.01, market.vol.value + vol_change);
                    m.rate.value = market.rate.value + rate_change;
                }
                return v0 - book.value(m, elapsed);
            };
            csv.row(s.name, s.start, s.end, days, spot_return, vol_change, rate_change, shocked(true, true),
                    shocked(true, false), shocked(false, true));
        }
    });
}
