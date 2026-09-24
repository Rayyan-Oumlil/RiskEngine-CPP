// Report 8.3: discretization bias against statistical error under Heston. The exact price comes
// from the characteristic function; Monte Carlo prices the same call with Andersen's QE scheme and
// with full-truncation Euler, for time steps from one year down to 1/64 of a year, on 2^21 paths
// each. Two regimes: a Feller-satisfied set (T = 1) and Andersen's (2008) case I, which violates
// Feller badly (2 kappa theta = 0.04 against xi^2 = 1) over T = 10, at the money and at K = 140.
//
// Each step size gets its own seed (SeedKey{2026, steps}), so that the bias estimates at different
// steps are independent: with one seed for all, the runs share most of their normals and every row
// inherits the same sampling error. Within a step size, QE and Euler share the seed.

#include <cmath>
#include <cstdint>
#include <string>

#include "harness/experiment.hpp"
#include "riskengine/methods/montecarlo/engine.hpp"
#include "riskengine/models/heston.hpp"
#include "riskengine/payoffs/vanilla.hpp"

using namespace riskengine;

namespace {

struct Case {
    const char* name;
    double rate, t, strike;
    HestonParams params;
    std::uint32_t steps_per_year_max;
};

template <class Scheme>
Estimate price(const Case& c, std::uint32_t steps) {
    const MarketState m{Spot{100}, Rate{c.rate}, Rate{0.0}, Vol{0.0}};
    const MonteCarlo<Scheme, VanillaPayoff> mc(VanillaPayoff{c.strike, OptionType::Call}, Maturity{c.t},
                                               MonteCarloConfig{.paths = 1u << 21, .steps = steps, .threads = 4},
                                               c.params);
    return mc.price(m, SeedKey{2026, steps});
}

} // namespace

int main(int argc, char** argv) {
    return harness::run("heston_discretization", argc, argv, [](harness::Experiment& exp) {
        const HestonParams feller{0.04, 2.0, 0.04, 0.3, -0.7};
        const HestonParams andersen{0.04, 0.5, 0.04, 1.0, -0.9};
        const Case cases[] = {
            {"feller_t1_k100", 0.05, 1.0, 100.0, feller, 64},
            {"andersen_t10_k100", 0.0, 10.0, 100.0, andersen, 32},
            {"andersen_t10_k140", 0.0, 10.0, 140.0, andersen, 32},
        };
        exp.param("spot", 100);
        exp.param("feller_set", "v0 = 0.04, kappa = 2, theta = 0.04, xi = 0.3, rho = -0.7, r = 0.05, T = 1");
        exp.param("andersen_case1", "v0 = 0.04, kappa = 0.5, theta = 0.04, xi = 1, rho = -0.9, r = 0, T = 10");
        exp.param("paths", 1u << 21);
        exp.param("seed", "SeedKey{2026, steps}: independent across step sizes, shared by the two schemes");
        exp.param("exact", "characteristic function (little trap), models/heston.hpp");

        auto csv = exp.csv({"case", "scheme", "steps", "dt", "price", "std_error", "exact", "bias", "bias_over_se"});
        for (const Case& c : cases) {
            const double exact = heston_price(OptionType::Call, 100, c.strike, c.rate, 0.0, c.t, c.params);
            for (std::uint32_t per_year = 1; per_year <= c.steps_per_year_max; per_year *= 2) {
                const auto steps = static_cast<std::uint32_t>(per_year * c.t);
                for (const char* scheme : {"qe", "euler"}) {
                    const Estimate e = std::string(scheme) == "qe" ? price<HestonQE>(c, steps) : price<HestonEuler>(c, steps);
                    csv.row(c.name, scheme, steps, c.t / steps, e.value, e.std_error, exact, e.value - exact,
                            (e.value - exact) / e.std_error);
                }
            }
        }
    });
}
