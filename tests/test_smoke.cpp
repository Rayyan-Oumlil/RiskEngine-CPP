#include <catch2/catch_test_macros.hpp>

#include "riskengine/core/market.hpp"
#include "riskengine/payoffs/vanilla.hpp"

TEST_CASE("Option and market construct with strong typedefs", "[smoke]") {
    using namespace riskengine;
    const VanillaOption option{Strike{100.0}, Maturity{1.0}, OptionType::Call};
    const MarketState market{Spot{100.0}, Rate{0.05}, Rate{0.0}, Vol{0.2}};
    REQUIRE(option.strike.value == 100.0);
    REQUIRE(option.type == OptionType::Call);
    REQUIRE(market.vol.value == 0.2);
}
