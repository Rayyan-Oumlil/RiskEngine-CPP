#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <span>
#include <string>

#include "riskengine/risk/backtest.hpp"

using namespace riskengine;

// Reference statistics and p-values computed independently with scipy.stats.chi2 (Python, from the
// likelihood-ratio definitions).

TEST_CASE("Kupiec POF against independent reference values", "[backtest]") {
    const TestResult red = kupiec_pof(250, 10, 0.01);
    CHECK(std::abs(red.statistic - 12.955491062356) < 1e-9);
    CHECK(std::abs(red.p_value - 0.000318984508) < 1e-11);
    const TestResult fine = kupiec_pof(250, 2, 0.01);
    CHECK(std::abs(fine.statistic - 0.108435216237) < 1e-9);
    CHECK(std::abs(fine.p_value - 0.741932700953) < 1e-11);
    const TestResult exact = kupiec_pof(1000, 10, 0.01); // exactly the expected rate
    CHECK(std::abs(exact.statistic) < 1e-9);
    CHECK(std::abs(exact.p_value - 1.0) < 1e-9);
    CHECK(kupiec_pof(250, 0, 0.01).statistic > 0.0); // no exception at all is also informative
}

TEST_CASE("Christoffersen independence rejects clustered exceptions", "[backtest]") {
    // Plain arrays: std::vector<bool> is not contiguous and cannot be viewed as a span<const bool>.
    bool buf[1000] = {};
    for (int d : {100, 101, 102, 103, 500}) buf[d] = true;
    const TestResult c = christoffersen_independence(std::span<const bool>(buf, 1000));
    CHECK(std::abs(c.statistic - 27.387651854265) < 1e-9);
    CHECK(std::abs(c.p_value - 1.664909609507e-07) < 1e-15);

    // The same five exceptions spread out: no evidence against independence.
    bool spread[1000] = {};
    for (int d : {100, 300, 500, 700, 900}) spread[d] = true;
    CHECK(christoffersen_independence(std::span<const bool>(spread, 1000)).p_value > 0.5);
    // Conditional coverage adds the POF statistic, on two degrees of freedom.
    const TestResult cc = christoffersen_conditional_coverage(std::span<const bool>(buf, 1000), 0.01);
    CHECK(std::abs(cc.statistic - (kupiec_pof(1000, 5, 0.01).statistic + c.statistic)) < 1e-12);
    CHECK(std::abs(cc.p_value - std::exp(-0.5 * cc.statistic)) < 1e-15);
}

TEST_CASE("Basel traffic-light zones", "[backtest]") {
    CHECK(basel_zone(0) == BaselZone::Green);
    CHECK(basel_zone(4) == BaselZone::Green);
    CHECK(basel_zone(5) == BaselZone::Yellow);
    CHECK(basel_zone(9) == BaselZone::Yellow);
    CHECK(basel_zone(10) == BaselZone::Red);
    CHECK(std::string(to_string(BaselZone::Red)) == "red");
}
