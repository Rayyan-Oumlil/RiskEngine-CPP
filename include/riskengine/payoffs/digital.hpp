#pragma once

#include "riskengine/payoffs/vanilla.hpp"

namespace riskengine {

// Cash-or-nothing: pays 1 if the option finishes strictly in the money (S > K for a call, S < K for
// a put), 0 at S == K. Discontinuous at the strike; the convention at the strike itself does not
// affect prices under a continuous distribution of S.
struct DigitalPayoff {
    double strike;
    OptionType type;

    double operator()(double s) const {
        return (type == OptionType::Call ? s > strike : s < strike) ? 1.0 : 0.0;
    }
};

} // namespace riskengine
