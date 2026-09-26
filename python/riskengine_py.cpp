// Python bindings (nanobind) over the header-only engine: `import riskengine`. Every function calls
// the same C++ code as the tests and experiments, so a Python price is the C++ price, bit for bit.
// This file is only the boundary: it validates arguments (a Python caller can pass anything; the
// engine itself checks its preconditions with debug-build asserts) and converts types.
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include "riskengine/core/simd.hpp"
#include "riskengine/methods/analytic/black_scholes.hpp"
#include "riskengine/methods/analytic/implied_vol.hpp"
#include "riskengine/methods/montecarlo/engine.hpp"
#include "riskengine/methods/montecarlo/longstaff_schwartz.hpp"
#include "riskengine/methods/tree/binomial.hpp"
#include "riskengine/models/gbm.hpp"
#include "riskengine/models/heston.hpp"
#include "riskengine/models/merton.hpp"
#include "riskengine/payoffs/vanilla.hpp"
#include "riskengine/risk/var.hpp"

namespace nb = nanobind;
using namespace nb::literals;
using namespace riskengine;

namespace {

OptionType option_type(std::string_view kind) {
    if (kind == "call") return OptionType::Call;
    if (kind == "put") return OptionType::Put;
    throw std::invalid_argument("kind must be 'call' or 'put', not '" + std::string(kind) + "'");
}

Exercise exercise_style(std::string_view exercise) {
    if (exercise == "european") return Exercise::European;
    if (exercise == "american") return Exercise::American;
    throw std::invalid_argument("exercise must be 'european' or 'american', not '" + std::string(exercise) + "'");
}

TreeMethod tree_method(std::string_view method) {
    if (method == "crr") return TreeMethod::Crr;
    if (method == "crr_averaged") return TreeMethod::CrrAveraged;
    if (method == "leisen_reimer") return TreeMethod::LeisenReimer;
    if (method == "bbs") return TreeMethod::Bbs;
    if (method == "bbs_richardson") return TreeMethod::BbsRichardson;
    throw std::invalid_argument("method must be one of crr, crr_averaged, leisen_reimer, bbs, bbs_richardson");
}

// The engine's own preconditions, checked here in every build.
void check_market(double spot, double strike, double maturity, double vol) {
    if (!(spot > 0.0)) throw std::invalid_argument("spot must be > 0");
    if (!(strike > 0.0)) throw std::invalid_argument("strike must be > 0");
    if (!(maturity >= 0.0)) throw std::invalid_argument("maturity must be >= 0");
    if (!(vol >= 0.0)) throw std::invalid_argument("vol must be >= 0");
}

void check_paths(std::uint64_t paths, unsigned threads) {
    if (paths < 2) throw std::invalid_argument("paths must be >= 2");
    if (threads < 1) throw std::invalid_argument("threads must be >= 1");
}

MarketState market(double spot, double rate, double dividend, double vol) {
    return MarketState{Spot{spot}, Rate{rate}, Rate{dividend}, Vol{vol}};
}

const char* simd_width() {
#if defined(__AVX512F__)
    return "avx512";
#elif defined(__AVX2__)
    return "avx2";
#else
    return "scalar";
#endif
}

} // namespace

NB_MODULE(riskengine, m) {
    m.doc() = "RiskEngine-CPP: options pricing and risk in C++20. Every function calls the C++ engine directly; "
              "stochastic results are reproducible bit for bit from their seed, whatever the thread count.";

    nb::class_<Estimate>(m, "Estimate", "A Monte Carlo result: value, standard error, known scheme bias, samples.")
        .def_ro("value", &Estimate::value)
        .def_ro("std_error", &Estimate::std_error)
        .def_ro("discretization", &Estimate::discretization, "Known or estimated scheme bias (e.g. LSM regression).")
        .def_ro("samples", &Estimate::samples)
        .def("__repr__", [](const Estimate& e) {
            return "Estimate(value=" + std::to_string(e.value) + ", std_error=" + std::to_string(e.std_error) +
                   ", discretization=" + std::to_string(e.discretization) + ", samples=" + std::to_string(e.samples) + ")";
        });

    m.def("simd_width", &simd_width, "Vector width this module was built with: 'avx512', 'avx2' or 'scalar'.");

    m.def(
        "black_scholes",
        [](double spot, double strike, double maturity, double vol, double rate, double dividend, const std::string& kind) {
            check_market(spot, strike, maturity, vol);
            return black_scholes_price(VanillaOption{Strike{strike}, Maturity{maturity}, option_type(kind)},
                                       market(spot, rate, dividend, vol));
        },
        "spot"_a, "strike"_a, "maturity"_a, "vol"_a, "rate"_a = 0.0, "dividend"_a = 0.0, "kind"_a = "call",
        "Closed-form Black-Scholes-Merton price with a continuous dividend yield.");

    m.def(
        "greeks",
        [](double spot, double strike, double maturity, double vol, double rate, double dividend, const std::string& kind) {
            check_market(spot, strike, maturity, vol);
            const Greeks g = black_scholes_greeks(VanillaOption{Strike{strike}, Maturity{maturity}, option_type(kind)},
                                                  market(spot, rate, dividend, vol));
            nb::dict d;
            d["delta"] = g.delta;
            d["gamma"] = g.gamma;
            d["vega"] = g.vega;
            d["theta"] = g.theta;
            d["rho"] = g.rho;
            return d;
        },
        "spot"_a, "strike"_a, "maturity"_a, "vol"_a, "rate"_a = 0.0, "dividend"_a = 0.0, "kind"_a = "call",
        "Closed-form Greeks in raw units: vega per 1.00 of vol, theta per year, rho per 1.00 of rate.");

    m.def(
        "implied_vol",
        [](double price, double spot, double strike, double maturity, double rate, double dividend, const std::string& kind) {
            check_market(spot, strike, maturity, 0.0);
            const ImpliedVolResult r = implied_vol(VanillaOption{Strike{strike}, Maturity{maturity}, option_type(kind)},
                                                   price, Spot{spot}, Rate{rate}, Rate{dividend});
            if (!r.ok()) throw std::invalid_argument(std::string("no implied vol: ") + to_string(r.status));
            return r.vol;
        },
        "price"_a, "spot"_a, "strike"_a, "maturity"_a, "rate"_a = 0.0, "dividend"_a = 0.0, "kind"_a = "call",
        "Implied volatility by a bracketed solver; raises ValueError when the price violates no-arbitrage bounds.");

    m.def(
        "binomial",
        [](double spot, double strike, double maturity, double vol, double rate, double dividend, const std::string& kind,
           const std::string& exercise, unsigned steps, const std::string& method) {
            check_market(spot, strike, maturity, vol);
            return binomial_price(tree_method(method), VanillaOption{Strike{strike}, Maturity{maturity}, option_type(kind)},
                                  exercise_style(exercise), market(spot, rate, dividend, vol), steps);
        },
        "spot"_a, "strike"_a, "maturity"_a, "vol"_a, "rate"_a = 0.0, "dividend"_a = 0.0, "kind"_a = "call",
        "exercise"_a = "european", "steps"_a = 1000, "method"_a = "crr",
        "Binomial tree price. Invalid lattices (probabilities outside (0, 1), wrong parity of steps) raise ValueError.");

    m.def(
        "monte_carlo",
        [](double spot, double strike, double maturity, double vol, double rate, double dividend, const std::string& kind,
           std::uint64_t paths, std::uint64_t seed, unsigned threads, bool antithetic) {
            check_market(spot, strike, maturity, vol);
            check_paths(paths, threads);
            if (antithetic && paths % 2 != 0) throw std::invalid_argument("antithetic needs an even number of paths");
            const MonteCarlo<GBM, VanillaPayoff> mc(VanillaPayoff{strike, option_type(kind)}, Maturity{maturity},
                                                    MonteCarloConfig{.paths = paths, .threads = threads, .antithetic = antithetic});
            nb::gil_scoped_release release; // pure C++ from here on: other Python threads may run
            return mc.price(market(spot, rate, dividend, vol), SeedKey{seed});
        },
        "spot"_a, "strike"_a, "maturity"_a, "vol"_a, "rate"_a = 0.0, "dividend"_a = 0.0, "kind"_a = "call",
        "paths"_a = 1u << 20, "seed"_a = 2026, "threads"_a = 1, "antithetic"_a = false,
        "Monte Carlo price of a European option under GBM, with its standard error. "
        "The result depends on the seed and the path count, never on the thread count.");

    m.def(
        "longstaff_schwartz",
        [](double spot, double strike, double maturity, double vol, double rate, double dividend, const std::string& kind,
           std::uint64_t paths, std::uint32_t steps, std::uint64_t seed, unsigned threads) {
            check_market(spot, strike, maturity, vol);
            check_paths(paths, threads);
            if (steps < 1) throw std::invalid_argument("steps must be >= 1");
            const LongstaffSchwartz<GBM> lsm(VanillaOption{Strike{strike}, Maturity{maturity}, option_type(kind)},
                                             LongstaffSchwartzConfig{.paths = paths, .steps = steps, .threads = threads});
            nb::gil_scoped_release release;
            return lsm.price(market(spot, rate, dividend, vol), SeedKey{seed});
        },
        "spot"_a, "strike"_a, "maturity"_a, "vol"_a, "rate"_a = 0.0, "dividend"_a = 0.0, "kind"_a = "put",
        "paths"_a = 1u << 17, "steps"_a = 50, "seed"_a = 2026, "threads"_a = 1,
        "American option by least-squares Monte Carlo. Estimate.discretization is the regression's low bias.");

    m.def(
        "heston",
        [](double spot, double strike, double maturity, double rate, double dividend, double v0, double kappa, double theta,
           double xi, double rho, const std::string& kind) {
            check_market(spot, strike, maturity, 0.0);
            return heston_price(option_type(kind), spot, strike, rate, dividend, maturity, HestonParams{v0, kappa, theta, xi, rho});
        },
        "spot"_a, "strike"_a, "maturity"_a, "rate"_a, "dividend"_a, "v0"_a, "kappa"_a, "theta"_a, "xi"_a, "rho"_a,
        "kind"_a = "call", "Heston stochastic-volatility price (characteristic function, little-trap form).");

    m.def(
        "merton",
        [](double spot, double strike, double maturity, double vol, double rate, double dividend, double jump_intensity,
           double jump_mean, double jump_sd, const std::string& kind) {
            check_market(spot, strike, maturity, vol);
            return merton_price(option_type(kind), spot, strike, rate, dividend, vol, maturity,
                                MertonParams{jump_intensity, jump_mean, jump_sd});
        },
        "spot"_a, "strike"_a, "maturity"_a, "vol"_a, "rate"_a, "dividend"_a, "jump_intensity"_a, "jump_mean"_a,
        "jump_sd"_a, "kind"_a = "call", "Merton jump-diffusion price (series of Black-Scholes prices).");

    m.def(
        "normal_var_es",
        [](double mean_loss, double sd_loss, double confidence) {
            if (!(sd_loss >= 0.0)) throw std::invalid_argument("sd_loss must be >= 0");
            if (!(confidence > 0.0 && confidence < 1.0)) throw std::invalid_argument("confidence must be in (0, 1)");
            const RiskMeasures r = normal_var_es(mean_loss, sd_loss, confidence);
            nb::dict d;
            d["var"] = r.var;
            d["es"] = r.es;
            return d;
        },
        "mean_loss"_a, "sd_loss"_a, "confidence"_a = 0.99, "Value-at-Risk and Expected Shortfall of a normal loss.");

    // Zero-copy: the C++ buffer becomes the NumPy array's memory; the capsule frees it with the array.
    m.def(
        "normals",
        [](std::size_t n, std::uint64_t seed, std::uint32_t stream, std::uint32_t block) {
            if (stream & kReservedStreamBits) throw std::invalid_argument("stream must be < 2**30 (top bits are reserved)");
            auto data = std::make_unique<double[]>(n);
            {
                nb::gil_scoped_release release;
                RandomStream rng(SeedKey{seed, stream}, block);
                rng.normals(std::span<double>(data.get(), n));
            }
            double* raw = data.release();
            nb::capsule owner(raw, [](void* p) noexcept { delete[] static_cast<double*>(p); });
            return nb::ndarray<nb::numpy, double, nb::ndim<1>>(raw, {n}, owner);
        },
        "n"_a, "seed"_a, "stream"_a = 0, "block"_a = 0,
        "n standard normals from the engine's reproducible generator (Philox 4x32-10 + AS241), as a NumPy array. "
        "The same (seed, stream, block) always gives the same numbers, on any machine and at any vector width.");
}
