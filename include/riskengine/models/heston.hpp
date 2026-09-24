#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include <cstddef>
#include <span>
#include <stdexcept>

#include "riskengine/core/market.hpp"
#include "riskengine/core/normal.hpp"
#include "riskengine/core/quadrature.hpp"
#include "riskengine/payoffs/vanilla.hpp"

namespace riskengine {

// Heston (1993) stochastic volatility under the risk-neutral measure:
//   dS = (r - q) S dt + sqrt(v) S dW_S,   dv = kappa (theta - v) dt + xi sqrt(v) dW_v,   d<W_S, W_v> = rho dt.
// The market state supplies spot, rate and dividend yield; its vol is not used (v0 replaces it).
struct HestonParams {
    double v0;    // initial variance
    double kappa; // mean-reversion speed
    double theta; // long-run variance
    double xi;    // volatility of variance
    double rho;   // spot-variance correlation

    // Feller condition 2 kappa theta >= xi^2: the variance process never reaches zero.
    bool feller() const { return 2.0 * kappa * theta >= xi * xi; }

    void validate() const {
        if (!(v0 >= 0.0 && kappa > 0.0 && theta > 0.0 && xi > 0.0 && rho >= -1.0 && rho <= 1.0))
            throw std::invalid_argument("Heston: needs v0 >= 0, kappa, theta, xi > 0 and |rho| <= 1");
    }
};

namespace detail {

// Characteristic function E[exp(i u ln S_T)] in the "little Heston trap" form of Albrecher, Mayer,
// Schoutens and Tistaert (2007): with the principal square root, g = (beta - d)/(beta + d) keeps
// |g exp(-dT)| < 1, so the complex logarithm never crosses its branch cut, whatever the maturity.
// The original formulation is discontinuous in u for long maturities and gives wrong prices there.
inline std::complex<double> heston_cf(std::complex<double> u, double log_forward, double t, const HestonParams& p) {
    using C = std::complex<double>;
    const C i(0.0, 1.0);
    const double xi2 = p.xi * p.xi;
    const C beta = p.kappa - p.rho * p.xi * i * u;
    const C d = std::sqrt(beta * beta + xi2 * (i * u + u * u));
    const C g = (beta - d) / (beta + d);
    const C e = std::exp(-d * t);
    const C big_d = (beta - d) / xi2 * (1.0 - e) / (1.0 - g * e);
    const C big_c = p.kappa * p.theta / xi2 * ((beta - d) * t - 2.0 * std::log((1.0 - g * e) / (1.0 - g)));
    return std::exp(big_c + big_d * p.v0 + i * u * log_forward);
}

} // namespace detail

// Price of a European call or put under Heston, from the characteristic function (Gil-Pelaez):
//   C = S e^{-qT} P1 - K e^{-rT} P2,  P_j = 1/2 + (1/pi) int_0^inf Re[e^{-iu ln K} f_j(u) / (iu)] du,
// with f_2(u) = phi(u) and f_1(u) = phi(u - i) / phi(-i); the put follows by parity. Both
// integrals are taken together, on panels scaled to the total variance, until they have decayed.
inline double heston_price(OptionType type, double s, double k, double r, double q, double t, const HestonParams& p) {
    p.validate();
    if (!(s > 0.0 && k > 0.0 && t > 0.0)) throw std::invalid_argument("Heston: needs spot, strike, maturity > 0");
    using C = std::complex<double>;
    const C i(0.0, 1.0);
    const double log_forward = std::log(s) + (r - q) * t, log_k = std::log(k);
    const double df_r = std::exp(-r * t), df_q = std::exp(-q * t);
    // e^{-rT} phi(u - i) - K e^{-rT} phi(u), over iu: the P1 and P2 integrands with their weights.
    auto integrand = [&](double u) {
        const C phase = std::exp(-i * u * log_k);
        const C f = df_r * (detail::heston_cf(C(u, -1.0), log_forward, t, p) - k * detail::heston_cf(C(u, 0.0), log_forward, t, p));
        return std::real(phase * f / (i * u));
    };
    const double total_var = std::max(p.v0, p.theta) * t;
    const double h = 0.5 / std::sqrt(std::max(total_var, 1e-6));
    const double call = 0.5 * (s * df_q - k * df_r) + integrate_half_line(integrand, h) / std::numbers::pi;
    return type == OptionType::Call ? call : call - s * df_q + k * df_r;
}

// Cash-or-nothing digital call paying 1 if S_T > K: e^{-rT} P2.
inline double heston_digital_call(double s, double k, double r, double q, double t, const HestonParams& p) {
    p.validate();
    using C = std::complex<double>;
    const C i(0.0, 1.0);
    const double log_forward = std::log(s) + (r - q) * t, log_k = std::log(k);
    auto integrand = [&](double u) {
        return std::real(std::exp(-i * u * log_k) * detail::heston_cf(C(u, 0.0), log_forward, t, p) / (i * u));
    };
    const double h = 0.5 / std::sqrt(std::max(std::max(p.v0, p.theta) * t, 1e-6));
    return std::exp(-r * t) * (0.5 + integrate_half_line(integrand, h) / std::numbers::pi);
}

// Path simulation, state (ln S, v); two normals per step (variance first, then spot).
struct HestonState {
    double log_spot;
    double var;
};

// Andersen's (2008) quadratic-exponential scheme. The variance step matches the first two moments
// of the exact non-central chi-squared transition: a squared Gaussian when the variance is large
// (psi <= 1.5), a mass at zero plus an exponential tail when it is near zero, so the variance is
// never negative and the process can reach zero, as it does when Feller fails. The log spot uses the
// exact integral of the variance's dynamics with trapezoidal weights (gamma1 = gamma2 = 1/2); no
// martingale correction, so the discounted spot carries a small O(dt) drift (measured in report 8.3).
class HestonQE {
public:
    using State = HestonState;
    using Params = HestonParams;
    static constexpr std::size_t factors = 2;

    HestonQE(const MarketState& m, const Params& p) : spot_(m.spot.value), carry_(m.rate.value - m.div.value), p_(p) {
        p.validate();
    }

    State initial_state() const { return {std::log(spot_), p_.v0}; }

    State step(State st, double dt, std::span<const double> z) const {
        const double e = std::exp(-p_.kappa * dt);
        const double m = p_.theta + (st.var - p_.theta) * e;
        const double s2 = st.var * p_.xi * p_.xi * e * (1.0 - e) / p_.kappa +
                          p_.theta * p_.xi * p_.xi * (1.0 - e) * (1.0 - e) / (2.0 * p_.kappa);
        const double psi = s2 / (m * m);
        double v_next;
        if (psi <= 1.5) {
            const double two_over_psi = 2.0 / psi;
            const double b2 = two_over_psi - 1.0 + std::sqrt(two_over_psi) * std::sqrt(two_over_psi - 1.0);
            const double a = m / (1.0 + b2);
            const double y = std::sqrt(b2) + z[0];
            v_next = a * y * y;
        } else {
            const double prob_zero = (psi - 1.0) / (psi + 1.0);
            const double beta = (1.0 - prob_zero) / m;
            // U = N(z); 1 - U = N(-z), computed directly to keep its precision in the upper tail.
            v_next = norm_cdf(z[0]) <= prob_zero ? 0.0 : std::log((1.0 - prob_zero) / norm_cdf(-z[0])) / beta;
        }
        const double rho_over_xi = p_.rho / p_.xi, one_minus_rho2 = 1.0 - p_.rho * p_.rho;
        const double k0 = -rho_over_xi * p_.kappa * p_.theta * dt;
        const double k1 = 0.5 * dt * (p_.kappa * rho_over_xi - 0.5) - rho_over_xi;
        const double k2 = 0.5 * dt * (p_.kappa * rho_over_xi - 0.5) + rho_over_xi;
        const double k34 = 0.5 * dt * one_minus_rho2 * (st.var + v_next);
        return {st.log_spot + carry_ * dt + k0 + k1 * st.var + k2 * v_next + std::sqrt(k34) * z[1], v_next};
    }

    double spot(State st) const { return std::exp(st.log_spot); }

private:
    double spot_, carry_;
    Params p_;
};

// Full-truncation Euler (Lord, Koekkoek and van Dijk 2010): the simplest scheme that keeps the
// square roots defined, v+ = max(v, 0) wherever v enters the dynamics. Its bias is the counterexample
// of report 8.3.
class HestonEuler {
public:
    using State = HestonState;
    using Params = HestonParams;
    static constexpr std::size_t factors = 2;

    HestonEuler(const MarketState& m, const Params& p)
        : spot_(m.spot.value), carry_(m.rate.value - m.div.value), p_(p) {
        p.validate();
    }

    State initial_state() const { return {std::log(spot_), p_.v0}; }

    State step(State st, double dt, std::span<const double> z) const {
        const double vp = std::max(st.var, 0.0);
        const double sd = std::sqrt(vp * dt);
        const double z_spot = p_.rho * z[0] + std::sqrt(1.0 - p_.rho * p_.rho) * z[1];
        return {st.log_spot + (carry_ - 0.5 * vp) * dt + sd * z_spot,
                st.var + p_.kappa * (p_.theta - vp) * dt + p_.xi * sd * z[0]};
    }

    double spot(State st) const { return std::exp(st.log_spot); }

private:
    double spot_, carry_;
    Params p_;
};

} // namespace riskengine
