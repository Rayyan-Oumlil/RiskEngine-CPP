// Report 3.1: implied-vol round trip on a grid covering the wings, quoting both the in-the-money
// and the out-of-the-money option at each point.
//
// Every Black-Scholes price is a difference of two terms, a - b (for a call, a = S e^{-qT} N(d1)
// and b = K e^{-rT} N(d2)). Each term carries a relative rounding error of about eps * (1 + d^2):
// in the Gaussian tail, the rounding of d itself is amplified by d. The quote is therefore only
// known to eps * (1 + d^2) * (a + b), and no solver can pin the vol better than the noise bound
//     eps * (1 + (1 + d^2) * (a + b) / (sigma * vega)),   d = max(|d1|, |d2|),
// in relative terms (the leading 1 is the resolution of sigma itself). For an in-the-money quote
// a + b is of the order of S + K while the time value is tiny: the bound, and the error, explode.

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

#include "harness/experiment.hpp"
#include "riskengine/methods/analytic/implied_vol.hpp"

using namespace riskengine;

int main(int argc, char** argv) {
    return harness::run("iv_roundtrip", argc, argv, [](harness::Experiment& exp) {
        constexpr double eps = std::numeric_limits<double>::epsilon();
        const double s = 100.0, r = 0.03, q = 0.01;
        const double strikes[] = {50, 60, 75, 90, 100, 110, 125, 150, 200};
        const double maturities[] = {1.0 / 365, 7.0 / 365, 30.0 / 365, 0.25, 1, 2, 5};
        const double vols[] = {0.05, 0.1, 0.2, 0.5, 1.0};
        exp.param("spot", s);
        exp.param("rate", r);
        exp.param("dividend_yield", q);
        exp.param("strikes", "50 60 75 90 100 110 125 150 200");
        exp.param("maturities_years", "1/365 7/365 30/365 0.25 1 2 5");
        exp.param("vols", "0.05 0.1 0.2 0.5 1.0");

        auto csv = exp.csv({"strike", "maturity", "vol", "type", "in_the_money", "price", "vega", "status",
                            "implied_vol", "rel_vol_error", "noise_bound"});
        for (double k : strikes)
            for (double t : maturities)
                for (double sigma : vols) {
                    const MarketState m{Spot{s}, Rate{r}, Rate{q}, Vol{sigma}};
                    const double df_r = std::exp(-r * t), df_q = std::exp(-q * t);
                    const bool call_itm = s * df_q > k * df_r;
                    for (OptionType type : {OptionType::Call, OptionType::Put}) {
                        const VanillaOption o{Strike{k}, Maturity{t}, type};
                        const bool itm = (type == OptionType::Call) == call_itm;
                        const double price = black_scholes_price(o, m);
                        const double vega = black_scholes_greeks(o, m).vega;
                        const ImpliedVolResult iv = implied_vol(o, price, Spot{s}, Rate{r}, Rate{q});
                        const bool solved = iv.status == ImpliedVolStatus::Ok;
                        const double sst = sigma * std::sqrt(t);
                        const double d1 = (std::log(s / k) + (r - q) * t) / sst + 0.5 * sst, d2 = d1 - sst;
                        const double terms = type == OptionType::Call
                                                 ? s * df_q * norm_cdf(d1) + k * df_r * norm_cdf(d2)
                                                 : k * df_r * norm_cdf(-d2) + s * df_q * norm_cdf(-d1);
                        const double d = std::max(std::abs(d1), std::abs(d2));
                        const double noise = eps * (1.0 + (1.0 + d * d) * terms / (sigma * vega));
                        csv.row(k, t, sigma, type == OptionType::Call ? "call" : "put", itm ? 1 : 0, price, vega,
                                solved ? "ok" : "zero_time_value", iv.vol,
                                solved ? std::abs(iv.vol - sigma) / sigma : 1.0, noise);
                    }
                }
    });
}
