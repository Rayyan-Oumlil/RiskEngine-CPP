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
(`days / 365.0`). The 252-trading-day convention is reserved for risk horizons: a one-day
risk horizon is 1/252 of a year, for the diffusion (σ√Δt) and for the time decay of the book
alike, so that a one-day loss never mixes trading and calendar time. Option maturities stay in
calendar time (a 30-day option has T = 30/365).

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
  log price. A quote only determines the vol to about ε·(1 + d²)·(a + b) / (σ·vega) relative, where
  a − b is the price formula; for an in-the-money quote the time value sits next to a much larger
  intrinsic value and this bound becomes large (report §3.1). A price below the discounted forward
  intrinsic value (beyond rounding) or at or above the upper bound is reported through
  `ImpliedVolStatus`, never as a number.
- **No `-ffast-math`** in reference builds. Tests compare floating-point results with relative
  tolerances, never exact equality, except where the result is exact by construction.

## Randomness and reproducibility

- **Generator:** Philox 4x32-10 (`core/rng/philox.hpp`), checked against the Random123
  known-answer vectors at compile time. Philox call *j* of block *b* uses key = seed and
  counter = (j, b, stream); its 128 bits give uniforms 2*j* and 2*j* + 1 of the block. There is no
  generator state to share or advance.
- **`SeedKey{seed, stream}`** identifies a random experiment. A stochastic pricer is a pure
  function of (market, key): re-running with the same key after a bump gives exact common random
  numbers. Use a different `stream` for an independent sub-simulation. Callers use streams below
  2³⁰: the top two bits are reserved for sub-streams the engine derives from a caller's key (control-
  variate pilot runs, QMC scrambling seeds), so they can never collide with a caller's own streams.
- **Uniforms** are the midpoints (k + ½)·2⁻⁵² of a 2⁵² grid: never 0 or 1, and 1 − u is exact.
- **Normals** come from the inverse CDF (Wichura AS241, `core/rng/normal_icdf.hpp`), never from
  `std::normal_distribution` (implementation-defined) or Box-Muller (breaks quasi-random points).
  Because the grid is symmetric, `normal(1 − u) == −normal(u)` exactly.
- **Blocks:** simulations are cut into a fixed number of blocks (`BlockPlan::blocks`, default 64),
  merged with Welford/Chan in block order. The result depends on the seed, the sample count and the
  block count, and **never** on the thread count.
- **Randomized QMC:** Sobol points (Joe-Kuo direction numbers, up to 1,024 dimensions) in Gray-code
  order, one coordinate per normal, Owen-scrambled per replication and mapped to (0, 1) as
  (x + ½)·2⁻³². The estimate is the mean of `replications` independent scramblings and its standard
  error their standard deviation over √replications.
- **What is bit-identical where:** uniforms, on every platform (integer arithmetic only). Normals
  and everything downstream, across thread counts always, and across GCC and Clang on Linux (same
  libm; `-ffp-contract=off` forbids silent FMA fusion, and `-fno-builtin-{exp,log,erfc,pow}` forbids
  the compiler from evaluating those functions itself on constants, with a rounding that differs from
  the run-time library's). MSVC's `log`/`exp` may differ in the last
  ulp, so on MSVC the golden normals are checked to a few ulps and the exact simulation golden
  value is skipped (its statistical check still runs).

## Experiments

Each executable in `experiments/` produces one figure or table of the report:
`data/results/<id>.csv` (full-precision `%.17g`, `\n` line endings) and `<id>.meta.json` (git
commit and dirty flag, compiler, build type, flags, parameters, UTC time). Results committed to the
repository must come from a Release build of a clean, committed tree (`"git_dirty": false`).
`tools/make_figures.py` turns them into `docs/figures/<id>.svg`.

## Market data

Historical series (NASDAQ Composite, VIX, VXO, 3-month T-bill yield) come from FRED and are
**frozen** in `data/raw/` by `tools/fetch_data.py`, byte for byte as served, with
`data/raw/manifest.json` recording the URL, extraction time, SHA-256 and date range of each file
(checked by `tests/test_time_series.cpp`). Experiments never download anything; refreshing the data
is a deliberate commit that changes the manifest. Days without a value (FRED's `.`) are skipped;
the index and vol series are joined on common dates and the rate is taken as of each date, since
its calendar differs; a rate or vol quoted in percent is divided by 100 at the point of use. Implied volatility is proxied by the VIX (VXO before 1990): real historical option
chains are not freely available, so the option books of report §7 are hypothetical.

## Reference values

Reference values are generated independently of the code under test, with mpmath
(`pip install -r tools/requirements.txt`, which pins the exact versions used for the committed
references and figures):

- `tools/bs_reference.py`: Black-Scholes prices and Greeks at 50 digits for
  `tests/test_black_scholes.cpp`. The Greeks are numerical derivatives of the high-precision price,
  so they check the closed-form formulas independently.
- `tools/normal_icdf_reference.py`: the inverse normal CDF at 60 digits by bisection, for
  `tests/test_rng.cpp`.
