#pragma once

namespace riskengine {

// Raw sensitivities, per unit of the input (see docs/conventions.md):
// vega per 1.00 of vol, rho per 1.00 of rate, theta = dV/dt per year.
struct Greeks {
    double delta;
    double gamma;
    double vega;
    double theta;
    double rho;
};

inline double vega_per_vol_point(const Greeks& g) { return g.vega / 100.0; }
inline double rho_per_rate_point(const Greeks& g) { return g.rho / 100.0; }
inline double theta_per_calendar_day(const Greeks& g) { return g.theta / 365.0; }

} // namespace riskengine
