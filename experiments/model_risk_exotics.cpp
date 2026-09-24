// Report 8.1: one at-the-money price, three models, diverging exotics. Black-Scholes (20 % vol),
// Heston and Merton are calibrated to the same one-year at-the-money call (S = K = 100, r = 5 %,
// q = 0); each model keeps a fixed shape (Heston: kappa = 1.5, xi = 0.6, rho = -0.7, v0 = theta
// solved; Merton: lambda = 0.25 per year, log jumps N(-0.2, 0.1^2), diffusion vol solved). Then:
//   model_risk_smile.csv    implied vol of the one-year smile from K = 60 to 150 (closed forms);
//   model_risk_exotics.csv  out-of-the-money vanillas, a digital, two daily-monitored barriers and a
//                           monthly Asian, by Monte Carlo (2^20 paths, 252 steps, SeedKey{2026}),
//                           with the closed form next to it where one exists.

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <span>
#include <string>

#include "harness/experiment.hpp"
#include "riskengine/methods/analytic/black_scholes.hpp"
#include "riskengine/methods/analytic/digital.hpp"
#include "riskengine/methods/analytic/implied_vol.hpp"
#include "riskengine/methods/montecarlo/engine.hpp"
#include "riskengine/models/gbm.hpp"
#include "riskengine/models/heston.hpp"
#include "riskengine/models/merton.hpp"
#include "riskengine/payoffs/digital.hpp"
#include "riskengine/payoffs/vanilla.hpp"

using namespace riskengine;

namespace {

constexpr double kSpot = 100, kRate = 0.05, kT = 1.0, kVol = 0.2;
constexpr std::uint32_t kSteps = 252;

// Path payoffs on the 252 daily fixings.
struct UpAndOutCall {
    double strike, barrier;
    double operator()(std::span<const double> f) const {
        for (double x : f)
            if (x >= barrier) return 0.0;
        return std::max(f.back() - strike, 0.0);
    }
};
struct DownAndInPut {
    double strike, barrier;
    double operator()(std::span<const double> f) const {
        const bool hit = std::any_of(f.begin(), f.end(), [&](double x) { return x <= barrier; });
        return hit ? std::max(strike - f.back(), 0.0) : 0.0;
    }
};
struct MonthlyAsianCall {
    double strike;
    double operator()(std::span<const double> f) const {
        double sum = 0.0;
        for (std::size_t j = 1; j <= 12; ++j) sum += f[21 * j - 1];
        return std::max(sum / 12.0 - strike, 0.0);
    }
};

// Root of an increasing function on [lo, hi] by bisection, to machine precision.
double solve(const std::function<double(double)>& f, double lo, double hi) {
    for (int i = 0; i < 200 && hi - lo > 1e-15 * hi; ++i) {
        const double mid = 0.5 * (lo + hi);
        (f(mid) < 0.0 ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

template <class Model, class Payoff, class... Params>
Estimate simulate(const Payoff& payoff, double vol, const Params&... params) {
    const MarketState m{Spot{kSpot}, Rate{kRate}, Rate{0.0}, Vol{vol}};
    const MonteCarlo<Model, Payoff> mc(payoff, Maturity{kT}, MonteCarloConfig{.paths = 1u << 20, .steps = kSteps, .threads = 4},
                                       params...);
    return mc.price(m, SeedKey{2026});
}

} // namespace

int main(int argc, char** argv) {
    return harness::run("model_risk_exotics", argc, argv, [](harness::Experiment& exp) {
        const VanillaOption atm{Strike{100}, Maturity{kT}, OptionType::Call};
        const double target = black_scholes_price(atm, MarketState{Spot{kSpot}, Rate{kRate}, Rate{0.0}, Vol{kVol}});

        // Calibration to the at-the-money price.
        const double heston_var = solve([&](double v) {
            return heston_price(OptionType::Call, kSpot, 100, kRate, 0, kT, {v, 1.5, v, 0.6, -0.7}) - target;
        }, 0.001, 0.2);
        const HestonParams heston{heston_var, 1.5, heston_var, 0.6, -0.7};
        const MertonParams merton{0.25, -0.2, 0.1};
        const double merton_vol = solve([&](double s) {
            return merton_price(OptionType::Call, kSpot, 100, kRate, 0, s, kT, merton) - target;
        }, 0.01, 0.5);

        exp.param("market", "S = 100, r = 0.05, q = 0, T = 1; calibrated to the ATM call at 20 % Black-Scholes vol");
        exp.param("atm_call_price", target);
        exp.param("heston", "kappa = 1.5, xi = 0.6, rho = -0.7, v0 = theta = heston_v0_theta");
        exp.param("heston_v0_theta", heston_var);
        exp.param("heston_feller", heston.feller());
        exp.param("merton", "lambda = 0.25, jump_mean = -0.2, jump_sd = 0.1, diffusion vol = merton_vol");
        exp.param("merton_vol", merton_vol);
        exp.param("monte_carlo", "2^20 paths, 252 daily steps, SeedKey{2026}; Heston by QE");

        // One-year smiles, from the closed forms.
        harness::Experiment smile("model_risk_smile", exp.out_dir());
        smile.param("source", "model_risk_exotics (same calibration)");
        auto sc = smile.csv({"model", "strike", "call_price", "implied_vol"});
        for (int k = 60; k <= 150; k += 5) {
            const VanillaOption o{Strike{double(k)}, Maturity{kT}, OptionType::Call};
            const double prices[3] = {
                black_scholes_price(o, MarketState{Spot{kSpot}, Rate{kRate}, Rate{0.0}, Vol{kVol}}),
                heston_price(OptionType::Call, kSpot, k, kRate, 0, kT, heston),
                merton_price(OptionType::Call, kSpot, k, kRate, 0, merton_vol, kT, merton)};
            const char* names[3] = {"black_scholes", "heston", "merton"};
            for (int i = 0; i < 3; ++i) {
                const ImpliedVolResult iv = implied_vol(o, prices[i], Spot{kSpot}, Rate{kRate}, Rate{0.0});
                sc.row(names[i], k, prices[i], iv.ok() ? iv.vol : std::numeric_limits<double>::quiet_NaN());
            }
        }
        smile.write_metadata();

        // Exotics by Monte Carlo, closed forms alongside where they exist.
        const double nan = std::numeric_limits<double>::quiet_NaN();
        auto csv = exp.csv({"product", "model", "price", "std_error", "closed_form"});
        auto three = [&](const char* product, auto payoff, double bs_cf, double heston_cf, double merton_cf) {
            using P = decltype(payoff);
            const Estimate b = simulate<GBM, P>(payoff, kVol);
            const Estimate h = simulate<HestonQE, P>(payoff, 0.0, heston);
            const Estimate m = simulate<Merton, P>(payoff, merton_vol, merton);
            csv.row(product, "black_scholes", b.value, b.std_error, bs_cf);
            csv.row(product, "heston", h.value, h.std_error, heston_cf);
            csv.row(product, "merton", m.value, m.std_error, merton_cf);
        };
        const MarketState bs_market{Spot{kSpot}, Rate{kRate}, Rate{0.0}, Vol{kVol}};
        auto vanilla = [&](double k, OptionType type, const char* name) {
            three(name, VanillaPayoff{k, type},
                  black_scholes_price(VanillaOption{Strike{k}, Maturity{kT}, type}, bs_market),
                  heston_price(type, kSpot, k, kRate, 0, kT, heston),
                  merton_price(type, kSpot, k, kRate, 0, merton_vol, kT, merton));
        };
        vanilla(100, OptionType::Call, "call_k100");
        vanilla(80, OptionType::Put, "put_k80");
        vanilla(120, OptionType::Call, "call_k120");
        three("digital_call_k100", DigitalPayoff{100, OptionType::Call},
              digital_price(VanillaOption{Strike{100}, Maturity{kT}, OptionType::Call}, bs_market),
              heston_digital_call(kSpot, 100, kRate, 0, kT, heston),
              merton_digital_call(kSpot, 100, kRate, 0, merton_vol, kT, merton));
        three("up_and_out_call_k100_b130", UpAndOutCall{100, 130}, nan, nan, nan);
        three("down_and_in_put_k100_b80", DownAndInPut{100, 80}, nan, nan, nan);
        three("asian_call_k100_monthly", MonthlyAsianCall{100}, nan, nan, nan);
    });
}
