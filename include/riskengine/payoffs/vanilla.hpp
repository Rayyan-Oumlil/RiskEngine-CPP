#pragma once

#include <algorithm>

#include "riskengine/core/units.hpp"

namespace riskengine {

enum class OptionType { Call, Put };

struct VanillaOption {
    Strike strike;
    Maturity maturity;
    OptionType type;
};

// max(0, phi (S - K)).
struct VanillaPayoff {
    double strike;
    OptionType type;

    double operator()(double s) const {
        return std::max(0.0, type == OptionType::Call ? s - strike : strike - s);
    }
};

// |S - K|: a call plus a put. Nearly even in the driving normal at the money, which is exactly
// where antithetic variates stop helping.
struct StraddlePayoff {
    double strike;

    double operator()(double s) const { return s > strike ? s - strike : strike - s; }
};

} // namespace riskengine
