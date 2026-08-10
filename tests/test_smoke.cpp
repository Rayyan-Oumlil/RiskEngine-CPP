#include <catch2/catch_test_macros.hpp>

#include "riskengine/types.hpp"

TEST_CASE("OptionSpec constructs with strong typedefs", "[smoke]") {
    using namespace riskengine;
    OptionSpec spec{Strike{100.0}, Maturity{1.0}, OptionType::Call};
    REQUIRE(spec.strike.value == 100.0);
    REQUIRE(spec.type == OptionType::Call);
}
