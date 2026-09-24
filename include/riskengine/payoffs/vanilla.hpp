#pragma once

#include "riskengine/core/units.hpp"

namespace riskengine {

enum class OptionType { Call, Put };

struct VanillaOption {
    Strike strike;
    Maturity maturity;
    OptionType type;
};

} // namespace riskengine
