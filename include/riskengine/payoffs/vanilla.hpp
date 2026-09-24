#pragma once

#include "riskengine/core/units.hpp"

namespace riskengine {

enum class OptionType { Call, Put };

struct VanillaOption {
    Strike strike;
    Maturity maturity;
    OptionType type;
};

// Terminal payoffs are templates on the scalar type (double for pricing, Dual for automatic
// differentiation) and expose derivative(s), the almost-everywhere derivative dPayoff/dS used by
// pathwise Greeks.

// max(0, phi (S - K)).
struct VanillaPayoff {
    double strike;
    OptionType type;

    template <class T>
    T operator()(const T& s) const {
        const T x = type == OptionType::Call ? s - strike : T(strike) - s;
        return x > 0.0 ? x : T(0.0);
    }

    double derivative(double s) const {
        if (type == OptionType::Call) return s > strike ? 1.0 : 0.0;
        return s < strike ? -1.0 : 0.0;
    }
};

// |S - K|: a call plus a put. Nearly even in the driving normal at the money, which is exactly
// where antithetic variates stop helping.
struct StraddlePayoff {
    double strike;

    template <class T>
    T operator()(const T& s) const {
        return s > strike ? s - strike : T(strike) - s;
    }

    double derivative(double s) const { return s > strike ? 1.0 : -1.0; }
};

} // namespace riskengine
