#pragma once

#include <cmath>

#include "riskengine/greeks/greeks.hpp"
#include "riskengine/risk/var.hpp"

// Local approximations of the one-period loss of a book driven by a single normal log-return
// x ~ N(0, s^2) of the underlying S over a horizon dt (report 7.2). They only need the book's
// Greeks; full revaluation (risk/var.hpp on repriced scenarios) is the reference they approximate.
namespace riskengine {

// Delta-normal: loss = -delta S x. Normal with mean 0 and sd |delta| S s. Blind to gamma and time.
inline RiskMeasures delta_normal(const Greeks& g, double spot, double s, double alpha) {
    return normal_var_es(0.0, std::abs(g.delta) * spot * s, alpha);
}

// Moments of the delta-gamma-theta loss L = -(delta S x + gamma S^2 x^2 / 2 + theta dt), exact for
// a normal x: the quadratic form's mean, variance, skewness and excess kurtosis.
struct LossMoments {
    double mean, sd, skew, excess_kurtosis;
};

inline LossMoments delta_gamma_moments(const Greeks& g, double spot, double s, double dt) {
    // L = c + b x + a x^2 with a = -gamma S^2 / 2, b = -delta S, c = -theta dt, x = s Z.
    const double a = -0.5 * g.gamma * spot * spot, b = -g.delta * spot, c = -g.theta * dt;
    const double s2 = s * s;
    const double mean = c + a * s2;
    const double var = b * b * s2 + 2.0 * a * a * s2 * s2;
    // Third and fourth central moments of b x + a (x^2 - s^2).
    const double m3 = 6.0 * a * b * b * s2 * s2 + 8.0 * a * a * a * s2 * s2 * s2;
    // E[(b x + a (x^2 - s^2))^4] = 3 b^4 s^4 + 60 a^2 b^2 s^6 + 60 a^4 s^8, written from 3 var^2.
    const double m4 = 3.0 * var * var + 48.0 * a * a * b * b * s2 * s2 * s2 + 48.0 * a * a * a * a * s2 * s2 * s2 * s2;
    const double sd = std::sqrt(var);
    return {mean, sd, var > 0.0 ? m3 / (var * sd) : 0.0, var > 0.0 ? m4 / (var * var) - 3.0 : 0.0};
}

// Delta-gamma-normal: the quadratic loss replaced by a normal with its exact mean and variance.
inline RiskMeasures delta_gamma_normal(const Greeks& g, double spot, double s, double dt, double alpha) {
    const LossMoments m = delta_gamma_moments(g, spot, s, dt);
    return normal_var_es(m.mean, m.sd, alpha);
}

// Delta-gamma Cornish-Fisher: the same moments plus skewness and excess kurtosis. Only the VaR is
// defined by the expansion; ES is left to the methods that have a distribution.
inline double delta_gamma_cornish_fisher_var(const Greeks& g, double spot, double s, double dt, double alpha) {
    const LossMoments m = delta_gamma_moments(g, spot, s, dt);
    return m.mean + m.sd * cornish_fisher_quantile(alpha, m.skew, m.excess_kurtosis);
}

} // namespace riskengine
