#pragma once

#include <cmath>

namespace riskengine {

// Forward-mode automatic differentiation: a value and its derivative along one direction.
// Evaluating a scalar-generic function on Dual{x, 1} returns f(x) and f'(x), exact to rounding,
// with no bump. Used to cross-check the hand-written pathwise Greeks (docs/riskengine_research.md
// 2.6): models and payoffs are templates on the scalar type, so the same code runs on double for
// pricing and on Dual for derivatives.
struct Dual {
    double v; // value
    double d; // derivative

    constexpr Dual(double value = 0.0, double derivative = 0.0) : v(value), d(derivative) {}

    static constexpr Dual variable(double x) { return {x, 1.0}; }
};

constexpr Dual operator+(Dual a, Dual b) { return {a.v + b.v, a.d + b.d}; }
constexpr Dual operator-(Dual a, Dual b) { return {a.v - b.v, a.d - b.d}; }
constexpr Dual operator-(Dual a) { return {-a.v, -a.d}; }
constexpr Dual operator*(Dual a, Dual b) { return {a.v * b.v, a.d * b.v + a.v * b.d}; }
constexpr Dual operator/(Dual a, Dual b) { return {a.v / b.v, (a.d * b.v - a.v * b.d) / (b.v * b.v)}; }

// Comparisons act on values: a branch of the payoff is chosen by the value, and its derivative is
// the derivative of that branch (so a kink or a jump contributes nothing, as in pathwise Greeks).
constexpr bool operator<(Dual a, Dual b) { return a.v < b.v; }
constexpr bool operator>(Dual a, Dual b) { return a.v > b.v; }
constexpr bool operator<=(Dual a, Dual b) { return a.v <= b.v; }
constexpr bool operator>=(Dual a, Dual b) { return a.v >= b.v; }

inline Dual exp(Dual a) {
    const double e = std::exp(a.v);
    return {e, e * a.d};
}
inline Dual log(Dual a) { return {std::log(a.v), a.d / a.v}; }
inline Dual sqrt(Dual a) {
    const double s = std::sqrt(a.v);
    return {s, a.d / (2.0 * s)};
}

} // namespace riskengine
