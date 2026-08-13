#pragma once

#include <cstddef>

namespace riskengine {

// Universal return type: no price ever travels as a bare double.
struct Estimate {
    double value;
    double std_error      = 0.0; // statistical error (0 for deterministic methods)
    double discretization = 0.0; // known or estimated scheme bias
    std::size_t samples   = 0;
};

} // namespace riskengine
