// Report 3.2: error of finite-difference delta and gamma against the closed form, as a function
// of the relative bump h. Truncation error falls as h^2 to the right; rounding error in the price
// rises as eps/h (delta) and eps/h^2 (gamma) to the left. The optimum sits near eps^{1/3}
// (delta) and eps^{1/4} (gamma): a smaller bump is not a more accurate one.

#include <cmath>
#include <format>
#include <limits>

#include "harness/experiment.hpp"
#include "riskengine/methods/analytic/black_scholes.hpp"

using namespace riskengine;

int main(int argc, char** argv) {
    return harness::run("fd_vcurve", argc, argv, [](harness::Experiment& exp) {
        const VanillaOption call{Strike{100.0}, Maturity{1.0}, OptionType::Call};
        const double s = 100.0, r = 0.05, q = 0.0, sigma = 0.2;
        const auto price_at = [&](double spot) {
            return black_scholes_price(call, MarketState{Spot{spot}, Rate{r}, Rate{q}, Vol{sigma}});
        };
        const Greeks exact = black_scholes_greeks(call, MarketState{Spot{s}, Rate{r}, Rate{q}, Vol{sigma}});

        constexpr double eps = std::numeric_limits<double>::epsilon();
        constexpr int per_decade = 8, first_exponent = -14, last_exponent = -1;
        exp.param("spot", s);
        exp.param("strike", 100.0);
        exp.param("rate", r);
        exp.param("dividend_yield", q);
        exp.param("vol", sigma);
        exp.param("maturity", 1.0);
        exp.param("exact_delta", exact.delta);
        exp.param("exact_gamma", exact.gamma);
        exp.param("machine_epsilon", eps);
        exp.param("h_grid", std::format("10^(k/{}) for k = {}..{}", per_decade, first_exponent * per_decade,
                                        last_exponent * per_decade));

        auto csv = exp.csv({"h_rel", "delta_fd", "delta_abs_error", "gamma_fd", "gamma_abs_error"});
        const double p0 = price_at(s);
        for (int k = first_exponent * per_decade; k <= last_exponent * per_decade; ++k) {
            const double h_rel = std::pow(10.0, static_cast<double>(k) / per_decade);
            const double up = s * (1.0 + h_rel), down = s * (1.0 - h_rel);
            // Divide by the step actually realized in floating point, not the intended one, so the
            // curve shows the rounding of the price rather than the rounding of the bump.
            const double half_step = 0.5 * (up - down);
            const double p_up = price_at(up), p_down = price_at(down);
            const double delta = (p_up - p_down) / (up - down);
            const double gamma = (p_up - 2.0 * p0 + p_down) / (half_step * half_step);
            csv.row(h_rel, delta, std::abs(delta - exact.delta), gamma, std::abs(gamma - exact.gamma));
        }
    });
}
