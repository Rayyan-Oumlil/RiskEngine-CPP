#pragma once

#include <concepts>

#include "riskengine/types.hpp"

namespace riskengine {

template <typename P>
concept Pricer = requires(const P p, const OptionSpec& o, const MarketData& m) {
    { p.price(o, m) } -> std::same_as<PriceResult>;
};

} // namespace riskengine
