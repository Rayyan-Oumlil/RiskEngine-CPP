// Local HTTP server exposing the pricing engine as JSON endpoints, for the static
// dashboard in server/static/. Not a production service: no auth, no TLS, no live
// market data. Run it, then open http://localhost:8080 in a browser.
#include <cstdio>
#include <exception>

#include <httplib.h>
#include <nlohmann/json.hpp>

#include "riskengine/methods/analytic/black_scholes.hpp"
#include "riskengine/methods/analytic/implied_vol.hpp"
#include "riskengine/methods/montecarlo/engine.hpp"
#include "riskengine/methods/montecarlo/greeks.hpp"
#include "riskengine/methods/tree/binomial.hpp"
#include "riskengine/models/gbm.hpp"
#include "riskengine/models/heston.hpp"
#include "riskengine/models/merton.hpp"
#include "riskengine/payoffs/vanilla.hpp"
#include "riskengine/risk/var.hpp"

using json = nlohmann::json;
using namespace riskengine;

namespace {

OptionType parse_option_type(const std::string& s) {
    if (s == "put") return OptionType::Put;
    return OptionType::Call;
}

Exercise parse_exercise(const std::string& s) {
    if (s == "american") return Exercise::American;
    return Exercise::European;
}

json price_endpoint(const json& body) {
    const double spot = body.at("spot").get<double>();
    const double strike = body.at("strike").get<double>();
    const double rate = body.value("rate", 0.05);
    const double div = body.value("dividend", 0.0);
    const double vol = body.at("vol").get<double>();
    const double maturity = body.at("maturity").get<double>();
    const OptionType type = parse_option_type(body.value("type", std::string{"call"}));
    const Exercise exercise = parse_exercise(body.value("exercise", std::string{"european"}));

    const MarketState market{Spot{spot}, Rate{rate}, Rate{div}, Vol{vol}};
    const VanillaOption option{Strike{strike}, Maturity{maturity}, type};

    json out;

    // Closed-form Black-Scholes and Greeks (used as the European reference for every method).
    const double bs = black_scholes_price(option, market);
    const Greeks g = black_scholes_greeks(option, market);
    out["black_scholes"] = {
        {"price", bs},
        {"delta", g.delta},
        {"gamma", g.gamma},
        {"vega_per_vol_point", vega_per_vol_point(g)},
        {"theta_per_day", theta_per_calendar_day(g)},
        {"rho", g.rho},
    };

    // Trees. American exercise only makes sense off the European closed form as a reference,
    // so we price both exercise styles here regardless of what the request asked for.
    const unsigned n = static_cast<unsigned>(body.value("tree_steps", 1000));
    const double crr = binomial_price(TreeMethod::Crr, option, Exercise::European, market, n);
    const unsigned lr_n = n % 2 == 0 ? n + 1 : n; // Leisen-Reimer requires an odd step count
    const double lr = binomial_price(TreeMethod::LeisenReimer, option, Exercise::European, market, lr_n);
    out["trees"] = {
        {"crr", {{"price", crr}, {"steps", n}, {"error_vs_bs", crr - bs}}},
        {"leisen_reimer", {{"price", lr}, {"steps", lr_n}, {"error_vs_bs", lr - bs}}},
    };
    if (exercise == Exercise::American) {
        const unsigned an = static_cast<unsigned>(body.value("american_steps", 2000));
        out["trees"]["american_bbs_richardson"] = binomial_price(TreeMethod::BbsRichardson, option, Exercise::American, market, an);
    }

    // Monte Carlo: price with an honest standard error, plus a pathwise delta.
    const auto paths = static_cast<std::size_t>(body.value("mc_paths", 1 << 18));
    const MonteCarlo<GBM, VanillaPayoff> mc(VanillaPayoff{strike, type}, Maturity{maturity},
                                            MonteCarloConfig{.paths = paths, .threads = 4, .antithetic = true});
    const Estimate mc_est = mc.price(market, SeedKey{2026});
    out["monte_carlo"] = {
        {"price", mc_est.value},
        {"std_error", mc_est.std_error},
        {"paths", paths},
        {"standard_errors_from_bs", mc_est.std_error > 0 ? (mc_est.value - bs) / mc_est.std_error : 0.0},
    };
    if (type == OptionType::Call) { // pathwise delta is defined for the vanilla call payoff here
        const Estimate delta_est =
            mc_greek(Greek::Delta, GreekMethod::Pathwise, VanillaPayoff{strike, type}, Maturity{maturity}, market,
                    SeedKey{7}, GreekConfig{.paths = paths, .threads = 4});
        out["monte_carlo"]["pathwise_delta"] = {{"value", delta_est.value}, {"std_error", delta_est.std_error}};
    }

    // Model risk: the same option under Heston (stochastic vol) and Merton (jumps), for
    // comparison against the GBM-based methods above. Illustrative default parameters.
    const HestonParams heston{vol * vol, 2.0, vol * vol, 0.3, -0.7};
    const MertonParams merton{0.25, -0.2, 0.1};
    out["model_risk"] = {
        {"heston_cf", heston_price(type, spot, strike, rate, div, maturity, heston)},
        {"merton_series", merton_price(type, spot, strike, rate, div, vol, maturity, merton)},
    };

    // Implied vol round trip on the Black-Scholes price, as a sanity check on the inputs.
    const ImpliedVolResult iv = implied_vol(option, bs, market.spot, market.rate, market.div);
    out["implied_vol_check"] = {{"recovered_vol", iv.vol}, {"status", to_string(iv.status)}};

    return out;
}

json var_endpoint(const json& body) {
    const double mu = body.value("mean_loss", 0.0);
    const double sd = body.at("std_loss").get<double>();
    const double alpha = body.value("confidence", 0.99);

    const RiskMeasures m = normal_var_es(mu, sd, alpha);
    return {
        {"confidence", alpha},
        {"var", m.var},
        {"expected_shortfall", m.es},
    };
}

} // namespace

int main() {
    httplib::Server svr;

    svr.set_default_headers({{"Access-Control-Allow-Origin", "*"}});

    svr.Post("/price", [](const httplib::Request& req, httplib::Response& res) {
        try {
            const json result = price_endpoint(json::parse(req.body));
            res.set_content(result.dump(), "application/json");
        } catch (const std::exception& e) {
            res.status = 400;
            res.set_content(json{{"error", e.what()}}.dump(), "application/json");
        }
    });

    svr.Post("/var", [](const httplib::Request& req, httplib::Response& res) {
        try {
            const json result = var_endpoint(json::parse(req.body));
            res.set_content(result.dump(), "application/json");
        } catch (const std::exception& e) {
            res.status = 400;
            res.set_content(json{{"error", e.what()}}.dump(), "application/json");
        }
    });

    svr.set_mount_point("/", RISKENGINE_STATIC_DIR);

    std::printf("RiskEngine-CPP dashboard: http://localhost:8080\n");
    svr.listen("0.0.0.0", 8080);
    return 0;
}
