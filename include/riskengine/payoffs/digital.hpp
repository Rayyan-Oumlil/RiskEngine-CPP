#pragma once

#include "riskengine/payoffs/vanilla.hpp"

namespace riskengine {

// Cash-or-nothing: pays 1 if the option finishes strictly in the money (S > K for a call, S < K for
// a put), 0 at S == K. Discontinuous at the strike; the convention at the strike itself does not
// affect prices under a continuous distribution of S.
struct DigitalPayoff {
    double strike;
    OptionType type;

    template <class T>
    T operator()(const T& s) const {
        return T((type == OptionType::Call ? s > strike : s < strike) ? 1.0 : 0.0);
    }

    // Zero almost everywhere: the jump at the strike is invisible to a pathwise derivative, which
    // is why the pathwise digital delta converges, confidently, to 0 (report 6).
    double derivative(double) const { return 0.0; }
};

} // namespace riskengine
