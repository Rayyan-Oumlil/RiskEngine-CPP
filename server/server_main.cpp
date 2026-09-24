// Local HTTP server exposing the pricing engine as JSON endpoints, for the static
// dashboard in server/static/. Not a production service: no auth, no TLS, no live
// market data. Run it, then open http://localhost:8080 in a browser.
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <exception>
#include <stdexcept>
#include <string>

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

// Reads an integer field, rejecting negative values and clamping to [1, max_value] so a
// malicious or malformed request can't force an O(n^2) tree induction or a multi-gigabyte
// Monte Carlo run on the server. Throws (caught by the request handler, returned as a 400)
// on a value that isn't an integer at all.
std::size_t bounded_size(const json& body, const char* key, std::size_t default_value, std::size_t max_value) {
    if (!body.contains(key)) return default_value;
    const json& v = body.at(key);
    if (!v.is_number_integer() || v.get<long long>() < 1)
        throw std::invalid_argument(std::string(key) + " must be a positive integer");
    return std::min(static_cast<std::size_t>(v.get<long long>()), max_value);
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
    // Capped at 20,000: the tree's backward induction is O(n^2), so an unbounded n from the
    // request could hang the server (docs/riskengine_research.md 3.2 uses 20,000 as its own
    // reference size, so this ceiling still covers every documented result).
    constexpr std::size_t kMaxTreeSteps = 20'000;
    const unsigned n = static_cast<unsigned>(bounded_size(body, "tree_steps", 1000, kMaxTreeSteps));
    const double crr = binomial_price(TreeMethod::Crr, option, Exercise::European, market, n);
    const unsigned lr_n = n % 2 == 0 ? n + 1 : n; // Leisen-Reimer requires an odd step count
    const double lr = binomial_price(TreeMethod::LeisenReimer, option, Exercise::European, market, lr_n);
    out["trees"] = {
        {"crr", {{"price", crr}, {"steps", n}, {"error_vs_bs", crr - bs}}},
        {"leisen_reimer", {{"price", lr}, {"steps", lr_n}, {"error_vs_bs", lr - bs}}},
    };
    if (exercise == Exercise::American) {
        const unsigned an = static_cast<unsigned>(bounded_size(body, "american_steps", 2000, kMaxTreeSteps));
        out["trees"]["american_bbs_richardson"] = binomial_price(TreeMethod::BbsRichardson, option, Exercise::American, market, an);
    }

    // Monte Carlo: price with an honest standard error, plus a pathwise delta. Capped at 2^22
    // paths (~4M) so one request can't force a multi-gigabyte, multi-minute simulation.
    const std::size_t paths = bounded_size(body, "mc_paths", 1 << 18, std::size_t{1} << 22);
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

    // Loopback only: this is a local tool with no authentication. Binding wider would expose
    // an unauthenticated service (and, combined with permissive CORS, invite DNS-rebinding
    // attacks from any page open in a browser on this machine) to the whole LAN.
    std::printf("RiskEngine-CPP dashboard: http://localhost:8080\n");
    svr.listen("127.0.0.1", 8080);
    return 0;
}
