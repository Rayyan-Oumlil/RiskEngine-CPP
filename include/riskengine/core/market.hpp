#pragma once

#include "riskengine/core/units.hpp"

namespace riskengine {

struct MarketState {
    Spot spot;
    Rate rate;
    Rate div; // continuous dividend yield
    Vol vol;
};

} // namespace riskengine
