#pragma once

namespace riskengine {

// Strong typedefs to kill the sigma-vs-sigma^2 class of bugs (see docs/riskengine_research.md 1.1).
struct Spot { double value; };
struct Strike { double value; };
struct Rate { double value; };
struct DividendYield { double value; };
struct Vol { double value; };
struct Maturity { double value; };

enum class OptionType { Call, Put };
enum class Exercise { European, American };

struct OptionSpec {
    Strike strike;
    Maturity maturity;
    OptionType type;
    Exercise exercise = Exercise::European;
};

struct MarketData {
    Spot spot;
    Rate rate;
    DividendYield dividend_yield;
    Vol vol;
};

struct PriceResult {
    double price;
    double std_error = 0.0; // nonzero for MC pricers
};

struct Greeks {
    double delta;
    double gamma;
    double vega;
    double theta;
    double rho;
};

} // namespace riskengine
