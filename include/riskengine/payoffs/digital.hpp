#pragma once

#include "riskengine/payoffs/vanilla.hpp"

namespace riskengine {

// Cash-or-nothing: pays 1 if the option finishes in the money. Discontinuous at the strike.
struct DigitalPayoff {
    double strike;
    OptionType type;

    double operator()(double s) const {
        return (type == OptionType::Call ? s > strike : s < strike) ? 1.0 : 0.0;
    }
};

} // namespace riskengine
