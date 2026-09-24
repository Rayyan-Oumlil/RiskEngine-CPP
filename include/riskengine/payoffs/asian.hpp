#pragma once

#include <algorithm>
#include <cmath>
#include <span>

#include "riskengine/payoffs/vanilla.hpp"

namespace riskengine {

// Asian options on the fixings S(t_1), ..., S(t_n): max(0, phi (A - K)) with A the average.

struct ArithmeticAsianPayoff {
    double strike;
    OptionType type;

    double operator()(std::span<const double> fixings) const {
        double sum = 0.0;
        for (double s : fixings) sum += s;
        const double avg = sum / static_cast<double>(fixings.size());
        return std::max(0.0, type == OptionType::Call ? avg - strike : strike - avg);
    }
};

// Pays 1 if the arithmetic average finishes strictly in the money. Its discontinuity, the surface
// {A = K}, is not aligned with any coordinate of the path, which is what makes it hard for
// quasi-Monte Carlo (report 5.3).
struct ArithmeticAsianDigitalPayoff {
    double strike;
    OptionType type;

    double operator()(std::span<const double> fixings) const {
        double sum = 0.0;
        for (double s : fixings) sum += s;
        const double avg = sum / static_cast<double>(fixings.size());
        return (type == OptionType::Call ? avg > strike : avg < strike) ? 1.0 : 0.0;
    }
};

// Geometric average, computed in log space. Has a closed form under GBM
// (methods/analytic/geometric_asian.hpp), hence its use as a control variate.
struct GeometricAsianPayoff {
    double strike;
    OptionType type;

    double operator()(std::span<const double> fixings) const {
        double log_sum = 0.0;
        for (double s : fixings) log_sum += std::log(s);
        const double avg = std::exp(log_sum / static_cast<double>(fixings.size()));
        return std::max(0.0, type == OptionType::Call ? avg - strike : strike - avg);
    }
};

} // namespace riskengine
