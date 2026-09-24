#pragma once

#include <concepts>

#include "riskengine/core/estimate.hpp"
#include "riskengine/core/market.hpp"

namespace riskengine {

// A pricer owns its instrument and prices it against a market state.
template <class P>
concept Pricer = requires(const P& p, const MarketState& m) {
    { p.price(m) } -> std::same_as<Estimate>;
};

} // namespace riskengine
