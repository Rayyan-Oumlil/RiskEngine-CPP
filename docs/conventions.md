# Conventions

Factor-of-100 and factor-of-365 errors are the most common bugs in this domain. Every number the
library takes or returns follows the units below; the strong typedefs in
`include/riskengine/core/units.hpp` make the inputs hard to mix up.

## Inputs

| Quantity | Type | Unit | Example |
|---|---|---|---|
| Spot, strike | `Spot`, `Strike` | currency units | `Spot{100.0}` |
| Volatility | `Vol` | decimal, annualized | `Vol{0.20}` is 20 % |
| Interest rate | `Rate` | decimal, annualized, **continuously compounded** | `Rate{0.05}` |
| Dividend yield | `Rate` (`MarketState::div`) | decimal, annualized, continuous | `Rate{0.02}` |
| Maturity | `Maturity` | years | `Maturity{0.25}` |

**Time.** Maturities are year fractions. When converting from days, use ACT/365 for calendar days
(`days / 365.0`). The 252-trading-day convention is reserved for scaling daily risk horizons
(Phase 6) and must be stated wherever it is used.

## Greeks

`Greeks` holds **raw** sensitivities, i.e. partial derivatives with respect to the input in the
units above:

| Greek | Definition | Unit |
|---|---|---|
| `delta` | ∂V/∂S | per 1 unit of spot |
| `gamma` | ∂²V/∂S² | per unit of spot, squared |
| `vega` | ∂V/∂σ | per **1.00** of vol (not per vol point) |
| `theta` | ∂V/∂t = −∂V/∂T | per **year**, calendar time |
| `rho` | ∂V/∂r | per **1.00** of rate |

Display conversions are explicit functions, never silent rescaling:

| Function | Meaning | Canonical value (S = K = 100, r = 5 %, σ = 20 %, T = 1) |
|---|---|---|
| `vega_per_vol_point(g)` | vega / 100: value change for σ + 1 point | 0.3752 |
| `rho_per_rate_point(g)` | rho / 100: value change for r + 1 point | 0.5323 |
| `theta_per_calendar_day(g)` | theta / 365 | −0.01757 |

## Degenerate inputs

Preconditions: spot > 0, strike > 0, σ ≥ 0, T ≥ 0 (checked with `assert` in debug builds).

When σ·√T = 0 (zero vol, expired option, or underflow of the product), there is no diffusion
left. The pricer then returns the **discounted forward intrinsic value**
max(0, φ(S e^{−qT} − K e^{−rT})), with φ = +1 for a call and −1 for a put. It never divides by zero.
In that branch:

- gamma = vega = 0;
- delta, theta and rho are those of the discounted forward intrinsic value;
- exactly at the money forward (S e^{−qT} = K e^{−rT}), N(d1) and N(d2) take their limit 1/2, so
  the call delta is e^{−qT}/2.

## Numerical conventions

- **Normal CDF:** N(x) = `0.5 * erfc(-x / sqrt(2))`. N(−x) is computed the same way, never as
  1 − N(x), which loses all relative precision in the tails.
- **Implied vol** is solved on the out-of-the-money quote (put-call parity) with Brent's method in
  log price. An in-the-money quote can only determine the vol to about ε·price / vega, because the
  time value sits next to a much larger intrinsic value. A price below the discounted forward
  intrinsic value (beyond rounding) or at or above the upper bound is reported through
  `ImpliedVolStatus`, never as a number.
- **No `-ffast-math`** in reference builds. Tests compare floating-point results with relative
  tolerances, never exact equality, except where the result is exact by construction.

## Reference values

`tools/bs_reference.py` regenerates the reference vectors used in `tests/test_black_scholes.cpp`
with mpmath at 50 digits. The Greeks there are numerical derivatives of the high-precision price,
so they check the closed-form formulas independently.
