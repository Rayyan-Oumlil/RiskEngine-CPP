#pragma once

#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>

#include "riskengine/core/market.hpp"
#include "riskengine/core/normal.hpp"
#include "riskengine/methods/analytic/black_scholes.hpp"
#include "riskengine/payoffs/vanilla.hpp"

namespace riskengine {

// Merton (1976) jump diffusion under the risk-neutral measure: GBM with volatility sigma (the market
// state's vol) plus a compound Poisson process of intensity lambda whose jumps multiply the spot by
// e^Y, Y ~ N(jump_mean, jump_sd^2). The drift is compensated by lambda k, k = E[e^Y] - 1, so the
// discounted spot stays a martingale.
struct MertonParams {
    double lambda;    // jump intensity, per year
    double jump_mean; // mean of the log jump
    double jump_sd;   // standard deviation of the log jump

    double mean_jump() const { return std::exp(jump_mean + 0.5 * jump_sd * jump_sd) - 1.0; } // k

    void validate() const {
        if (!(lambda >= 0.0 && jump_sd >= 0.0)) throw std::invalid_argument("Merton: needs lambda, jump_sd >= 0");
    }
};

// Merton's series: conditional on n jumps the log spot is Gaussian, so the price is a Poisson
// mixture of Black-Scholes prices with vol sigma_n^2 = sigma^2 + n jump_sd^2 / T and rate
// r_n = r - lambda k + n ln(1 + k) / T, weighted by the Poisson(lambda (1 + k) T) probabilities.
// Summed until the remaining Poisson mass is below 1e-16.
inline double merton_price(OptionType type, double s, double k, double r, double q, double sigma, double t,
                           const MertonParams& p) {
    p.validate();
    const double kbar = p.mean_jump();
    const double intensity = p.lambda * (1.0 + kbar) * t; // lambda' T
    double weight = std::exp(-intensity), mass = 0.0, price = 0.0;
    for (int n = 0; n < 10000; ++n) {
        if (n > 0) weight *= intensity / n;
        const double sigma_n = std::sqrt(sigma * sigma + n * p.jump_sd * p.jump_sd / t);
        const double r_n = r - p.lambda * kbar + n * std::log1p(kbar) / t;
        price += weight * detail::bs_price(type, s, k, r_n, q, sigma_n, t);
        mass += weight;
        if (n > intensity && 1.0 - mass < 1e-16) return price;
    }
    throw std::runtime_error("merton_price: Poisson series did not converge");
}

// Exact simulation over any step: the number of jumps N ~ Poisson(lambda dt), their summed log size
// N jump_mean + jump_sd sqrt(N) Z, and the diffusion. Three normals per step: diffusion, the Poisson
// draw (through U = N(z) and inversion of the Poisson distribution), and the jump sizes. State: ln S.
class Merton {
public:
    using State = double;
    using Params = MertonParams;
    static constexpr std::size_t factors = 3;

    Merton(const MarketState& m, const Params& p)
        : spot_(m.spot.value),
          drift_(m.rate.value - m.div.value - p.lambda * p.mean_jump() - 0.5 * m.vol.value * m.vol.value),
          vol_(m.vol.value),
          p_(p) {
        p.validate();
    }

    State initial_state() const { return std::log(spot_); }

    State step(State x, double dt, std::span<const double> z) const {
        const int jumps = poisson_inverse(norm_cdf(z[1]), p_.lambda * dt);
        return x + drift_ * dt + vol_ * std::sqrt(dt) * z[0] + jumps * p_.jump_mean +
               p_.jump_sd * std::sqrt(static_cast<double>(jumps)) * z[2];
    }

    double spot(State x) const { return std::exp(x); }

    // Smallest n with P(N <= n) >= u, for N ~ Poisson(mean).
    static int poisson_inverse(double u, double mean) {
        double pmf = std::exp(-mean), cdf = pmf;
        int n = 0;
        while (cdf < u && n < 1000) {
            ++n;
            pmf *= mean / n;
            cdf += pmf;
        }
        return n;
    }

private:
    double spot_, drift_, vol_;
    Params p_;
};

} // namespace riskengine
