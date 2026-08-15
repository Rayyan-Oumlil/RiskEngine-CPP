#pragma once

namespace riskengine {

// Strong typedefs to kill the sigma-vs-sigma^2 and percent-vs-decimal classes of bugs.
// Units are fixed by docs/conventions.md.
struct Spot     { double value; };
struct Strike   { double value; };
struct Vol      { double value; }; // decimal, annualized (0.20 = 20 %)
struct Rate     { double value; }; // decimal, continuously compounded, annualized
struct Maturity { double value; }; // years

} // namespace riskengine
