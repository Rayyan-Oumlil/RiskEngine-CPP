# RiskEngine-CPP — Numerical-method, risk-measure and model risk report

> **Status:** in progress. §2, §3 and §5 are written (Phases 1, 2 and 4). The other sections are
> filled in as the phases of the plan land ([`riskengine_research.md`](riskengine_research.md) §12).
>
> **Editorial rules.** Every figure has a self-contained caption (*what it shows*, then *why it
> matters*). Every number points to the CSV or test that produced it. Every stochastic estimate is
> given with its standard error. Every negative claim ("X fails") comes with its mechanism and its
> remedy.

## Executive summary

*To be written in Phase 8, once the quantified findings of §3 to §8 are established.*

---

## 1. Introduction and risk taxonomy

*To be written in Phase 8.* Framework: [`riskengine_research.md`](riskengine_research.md) §0.

---

## 2. Setup, conventions and reproducibility

### 2.1 Dynamics, reference parameters, unit conventions

Units (continuously compounded decimal rates, annualized decimal vol, maturities in years, raw
Greeks per unit of input) and the treatment of degenerate cases are fixed in
[`conventions.md`](conventions.md). The reference vector of this report is S = K = 100, r = 5 %,
q = 0, σ = 20 %, T = 1 year.

### 2.2 Protocol: generator, seeds, efficiency metric, regeneration

**Generator.** Philox 4×32-10, a counter-based generator: call *j* of block *b* is
Philox(key = seed, counter = (j, b, stream)), and its 128 bits give uniforms 2*j* and 2*j* + 1.
There is no state to share or advance. The
implementation is checked at compile time against the Random123 known-answer vectors. Uniforms
are the midpoints (k + ½)·2⁻⁵² of a 2⁵²-point grid: never 0 or 1, and 1 − u is exact. Normals
come from inversion (Wichura AS241, relative error ≤ 10⁻¹⁵ against a 60-digit reference), never
from `std::normal_distribution`, whose algorithm depends on the standard library, nor from
Box-Muller, which destroys the structure of quasi-random points. Because the grid is symmetric,
`normal(1 − u) == −normal(u)` exactly.

**Bit-for-bit reproducibility.** A simulation is cut into a fixed number of blocks, independent of
the thread count. Each block owns its random stream and its Welford accumulator; the partial
results are merged (Chan) in block order. Floating-point addition is not associative, so it is this
fixed order that makes the result bit-identical with 1, 2, 3, 8 or 64 threads. Automatic FMA
contraction is disabled (`-ffp-contract=off`), so GCC and Clang produce the same bits, including
with `-march=native`. So is compile-time evaluation of `exp`, `log`, `erfc` and `pow`
(`-fno-builtin-*`): GCC evaluates them on constants with correctly rounded arithmetic while the
run-time library may differ in the last ulp, so the bits of an expression depended on what the
optimizer happened to inline. This was caught when adding unrelated code to the engine moved the
exact Black-Scholes price recorded by an experiment in its 15th digit, and it changed the
rounding-dominated half of the V-curve of §3.2 (its figure and table are from the fixed build). A Monte Carlo estimate is therefore a pure function of (seed, number of
draws, number of blocks).

**Validation gates (Phase 2)**, all in `tests/test_rng.cpp` and `tests/test_simulation.cpp`, with
fixed seeds (deterministic tests):

| Gate | Result |
|---|---|
| Philox against the Random123 vectors | 3/3, checked by `static_assert` |
| Uniforms: mean, variance, Kolmogorov-Smirnov (10⁶ draws) | within 4 standard deviations; √n·D < 1.95 (99.9 % level) |
| Normals: mean, variance, skewness, kurtosis, Kolmogorov-Smirnov (10⁶ draws) | same |
| Serial (lag-1) and cross-block correlation | < 4/√n |
| Same result with 1/2/3/8/64 threads | exact equality |
| Same result on GCC and Clang | exact equality on golden values |
| Standard-error calibration, 200 independent replications | ≥ 90 % of \|z\| ≤ 1.96, no \|z\| > 4 |

**Efficiency.** Variance-reduction techniques will be compared by
efficiency = 1 / (variance × CPU time) (Glasserman), from Phase 4 on.

**Regeneration.** Every figure or table is produced by an executable in `experiments/`, which
writes `data/results/<id>.csv` and `<id>.meta.json`: git commit and dirty-tree flag, compiler,
configuration, flags, parameters, date. Committed results come from a Release build of a clean tree
(`"git_dirty": false`). Figures are rendered by `tools/make_figures.py`:

```
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release
for e in fd_vcurve iv_roundtrip mc_convergence mc_coverage mc_efficiency qmc_convergence; do build-release/experiments/$e; done
python3 tools/make_figures.py        # pip install -r tools/requirements.txt
```

---

## 3. Ground truth: analytic Black-Scholes

Every method in the following sections (trees, Monte Carlo, stochastic Greeks) is judged against
the closed form. This section establishes that the closed form itself is exact to the precision
shown, and documents the two numerical limits it already imposes: inversion to implied volatility,
and finite differences.

### 3.1 Validation, invariants, implied volatility

**Reference values.** Ten cases (canonical vector, continuous dividend, Hull's example 15.6, and
two wing strikes whose out-of-the-money price is ~10⁻¹²) are compared with a 50-digit mpmath
reference ([`tools/bs_reference.py`](../tools/bs_reference.py)). The reference Greeks are
*numerical derivatives of the high-precision price*, not the closed forms, so they independently
check the formulas, their signs and the dividend handling. The price and all five Greeks agree to
10⁻⁹ relative in all ten cases, including the 10⁻¹² wing prices (`tests/test_black_scholes.cpp`).

| Canonical vector | Value |
|---|---|
| Call / put | 10.450584 / 5.573526 |
| Call delta | 0.636831 |
| Gamma | 0.018762 |
| Vega | 37.524 per 1.00 of vol (0.375 per point) |
| Call theta | −6.414 per year (−0.01757 per calendar day) |
| Call rho | 53.232 per 1.00 of rate |

**Invariants**, on a 315-point grid (strike 50 to 200 for S = 100, maturity 1 day to 5 years,
vol 5 % to 100 %): put-call parity to 10⁻¹⁴ relative to S + K; no-arbitrage bounds; ∂C/∂K < 0 and
convexity in strike; closed-form Greeks consistent with finite differences of the price. As T → 0
and σ → 0, the pricer returns the discounted forward intrinsic value through an explicit branch,
never through a division by zero; no NaN appears down to T = σ = 10⁻³⁰⁰.

**Implied volatility.** The solver maps the quote to the out-of-the-money option by parity, then
solves log price(σ) = log target with Brent's method on a bracket. Working in log price keeps the
problem well conditioned in the wings, where the price spans hundreds of orders of magnitude.
Newton alone is not used: it diverges where vega → 0.

![Relative error of the recovered implied vol against the noise bound of the quote](figures/iv_roundtrip.svg)

*Figure 1 — Price → vol → price round trip on the 315-point grid, for the out-of-the-money (blue)
and in-the-money (orange) quote at each point. On the x-axis, the noise bound
ε·(1 + (1 + d²)(a + b)/(σ·vega)), where a − b is the price formula: the relative vol error that
the rounding of the quote causes on its own. **What it shows:** every point lies below the
diagonal; the solver's error never exceeds the noise in the quote, and the in-the-money degradation
tracks that noise exactly. **Why it matters:** an inaccurate implied vol taken from an
in-the-money quote is not a solver defect but information lost in the quote itself; no solver can
recover it.* Source: [`data/results/iv_roundtrip.csv`](../data/results/iv_roundtrip.csv).

| Quote | Solved | Max relative error | Above 10⁻⁹ | Max error / bound |
|---|---|---|---|---|
| Out of the money | 289 / 315 | 1.3 × 10⁻¹³ | 0 | 0.68 |
| In the money | 226 / 315 | 5.7 × 10⁻² | 19 | 0.64 |

The 26 unsolved out-of-the-money quotes have a price that underflows to zero (for example 1 day,
K/S = 2, σ = 5 %): there is no vol to recover. The 89 unsolved in-the-money quotes have a time
value below the rounding of their intrinsic value; the solver reports this (`ZeroTimeValue`)
instead of returning a number. The worst solved case, a put with K = 150, 3 months and σ = 10 %,
recovers the vol to within 5.7 %.

**Mechanism.** A Black-Scholes price is a difference a − b of two terms (for a call,
a = S e^{−qT} N(d₁) and b = K e^{−rT} N(d₂)). In the Gaussian tail the rounding of d is amplified
by a factor d, so each term carries a relative error of about ε·(1 + d²). The quote is therefore
only known to within ε·(1 + d²)·(a + b). Out of the money, a + b stays of the order of the price
and the vol is recovered to ~10⁻¹³. In the money, a + b ≈ S + K while the time value is tiny: the
loss reaches several percent.

**Remedy.** Always invert the out-of-the-money quote; an implied vol taken from a deep in-the-money
quote must come with its noise bound. Jäckel's formulation (*Let's Be Rational*, 2015) remains the
target for reducing the pricer's own share of the error in the wings.

### 3.2 Floating-point rounding and the finite-difference V-curve

![Absolute error of finite-difference delta and gamma against the relative bump](figures/fd_vcurve.svg)

*Figure 2 — Absolute error of central-difference delta and gamma on the Black-Scholes price
(canonical vector), against the closed form, for a relative bump h from 10⁻¹⁴ to 10⁻¹. Dashed lines
mark the theoretical optima ε^{1/3} (delta) and ε^{1/4} (gamma). **What it shows:** the error is
V-shaped; on the right, truncation falls as h²; on the left, price rounding grows as ε/h for delta
and ε/h² for gamma. **Why it matters:** a smaller bump does not make a Greek more accurate; below
the optimum it makes it wrong, and for gamma catastrophically so.*
Source: [`data/results/fd_vcurve.csv`](../data/results/fd_vcurve.csv).

In the rounding regime the error at any single h is a noisy draw (it changes sign and can cancel by
luck), so the table gives the median over the seven grid points within ±⅜ of a decade of each h.

| Relative bump h | Delta error | Gamma error (relative) |
|---|---|---|
| 10⁻³ | 8.6 × 10⁻⁷ | 1.2 × 10⁻⁶ |
| 10⁻⁴ | 8.6 × 10⁻⁹ | 2.4 × 10⁻⁸ |
| 10⁻⁵ | 8.1 × 10⁻¹¹ | 1.4 × 10⁻⁶ |
| 10⁻⁶ | 1.1 × 10⁻¹¹ | 1.1 × 10⁻⁴ |
| 10⁻⁸ | 2.0 × 10⁻⁹ | 130 % |
| 10⁻¹⁴ | 9.7 × 10⁻⁴ | 5.6 × 10¹¹ |

The delta error floor is about 10⁻¹¹ for h between 3 × 10⁻⁷ and ε^{1/3} ≈ 6 × 10⁻⁶; the gamma
floor is about 4 × 10⁻¹⁰ (2 × 10⁻⁸ relative) near h ≈ 5 × 10⁻⁵, just below ε^{1/4} ≈ 1.2 × 10⁻⁴.
Isolated downward spikes (for example 4 × 10⁻¹⁵ for delta at h = 2.4 × 10⁻⁶) are sign changes of
the error, not attainable accuracy. From h = 10⁻⁸ down, all 46 grid points give a gamma that is
100 % wrong or worse.

**Mechanism.** For the central difference, the truncation error is ~(h²S²/6)·∂³V/∂S³. The rounding
error is not ε·V: as in §3.1, the price is a difference a − b of two terms that each carry a relative
error of a few ε, so it is only known to ~ε·(a + b) ≈ 2.6 × 10⁻¹⁴ here (a + b ≈ 117 for a price of
10.45). The rounding error is therefore ~ε(a + b)/(hS) for delta and ~ε(a + b)/(hS)² for gamma, and
their sum with the truncation error is smallest at h* ∝ ε^{1/3} and ε^{1/4} respectively. At
h = 10⁻⁵ this predicts a gamma error of 2.6 × 10⁻¹⁴ / 10⁻⁶ ≈ 1.4 × 10⁻⁶ relative, which is what the
table shows. The experiment divides by the step actually realized in floating point, not the intended
one, so that the curve shows only the rounding of the price.

**Remedy.** Use a relative bump, never an absolute one, close to the optimum for the order of the
derivative. A relative bump of 10⁻⁴, a common choice, is safe for both Greeks on a deterministic
pricer. On a Monte Carlo pricer, statistical noise replaces rounding and moves the optimum by
several orders of magnitude: that is the subject of §6.

---

## 4. Trees: convergence and pathologies

*To be written in Phase 3.*

## 5. Monte Carlo: convergence and variance reduction

The engine (`methods/montecarlo/engine.hpp`) simulates any `PathModel` (GBM so far, with an exact
log-space step) for terminal and path-dependent payoffs, discounts, and returns an `Estimate` with its
standard error. Every experiment below uses the canonical market of §2.1 and its closed-form
prices as the truth. The engine's tests hold at 4 standard errors for calls, puts and digitals
with 1 and 12 time steps, and for the straddle and the geometric Asian with 12 steps
(`tests/test_monte_carlo.cpp`).

### 5.1 The N^{−1/2} rate and interval coverage

![Standard error and actual error of the Monte Carlo price against the number of paths](figures/mc_convergence.svg)

*Figure 3 — Monte Carlo price of the canonical ATM call with N = 2¹⁰ to 2²² paths (an independent
seed for each N): standard error and actual error against Black-Scholes. **What it shows:** the
standard error falls as N^{−1/2} (fitted slope −0.5004), and the actual error scatters below and
around it, within 1.55 standard errors at every N. **Why it matters:** the rate is the textbook one,
so any faster apparent convergence in later sections has to come from variance reduction, not luck.*
Source: [`data/results/mc_convergence.csv`](../data/results/mc_convergence.csv).

A slope proves the standard error has the right *rate*. It does not prove it has the right *size*.
For that, 1,000 independent pricings (10,000 paths each) of the ATM call and of an out-of-the-money
digital (K = 130) were compared with their exact prices.

![Histograms of the z-scores of 1,000 independent Monte Carlo prices against the standard normal density](figures/mc_coverage.svg)

*Figure 4 — z-scores (estimate − exact) / SE of 1,000 independent pricings, against the N(0, 1)
density. **What it shows:** 94.5 % (call) and 95.7 % (digital) of the 95 % confidence intervals
contain the exact price, against 95 % ± 0.69 % expected; the z-scores have standard deviation 1.04
and 1.00. **Why it matters:** the error bars reported with every Monte Carlo number in this report
are neither optimistic nor conservative, including for the skewed Bernoulli payoff of a digital.*
Source: [`data/results/mc_coverage.csv`](../data/results/mc_coverage.csv).

### 5.2 Efficiency table, failures included

Techniques are compared by **efficiency = 1 / (variance per path × CPU time per path)**, relative to
plain Monte Carlo on the same payoff. "Variance per path" charges each technique for every path it
simulates (an antithetic pair counts as two). Times are single-threaded, the minimum of three runs,
and include the pilot run that estimates the control-variate coefficient. 10⁶ paths per row.

| Payoff | Technique | Variance per path | ns per path | Efficiency vs plain |
|---|---|---|---|---|
| ATM call | plain | 215.8 | 39.6 | 1 |
| | antithetic | 108.0 | 27.1 | 2.9 |
| | S_T control | 31.5 | 40.1 | 6.8 |
| | antithetic + S_T control | 7.6 | 27.4 | **41** |
| Deep ITM call, K = 70 | plain | 397.4 | 26.3 | 1 |
| | antithetic | 23.1 | 17.8 | 25 |
| | S_T control | 0.98 | 30.1 | 355 |
| | antithetic + S_T control | 0.43 | 20.5 | **1,200** |
| Far OTM call, K = 160 | plain | 3.43 | 26.2 | 1 |
| | antithetic | 3.51 | 17.3 | 1.5 |
| | S_T control | 3.11 | 26.5 | **1.1** |
| ATM straddle | plain | 174.6 | 37.6 | 1 |
| | antithetic | 287.7 | 25.4 | **0.90** |
| Arithmetic Asian, ATM, 12 fixings | plain | 72.4 | 241.4 | 1 |
| | antithetic | 34.8 | 145.8 | 3.4 |
| | geometric control | 0.056 | 297.8 | **1,039** |
| | antithetic + geometric control | 0.061 | 189.4 | 1,511 |

Source: [`data/results/mc_efficiency.csv`](../data/results/mc_efficiency.csv) (timings on the CPU
recorded in its metadata; all estimates are within 2 standard errors of the exact price where one
exists).

**Antithetic variates** pair each path with its reflection Z → −Z. The variance per path becomes
σ²(1 + ρ), where ρ is the correlation between the two halves of a pair. For a payoff monotone in Z,
ρ < 0: the ATM call halves its variance, and the deep ITM call, which is almost linear in S_T, cuts it
17-fold. They also save time: a pair reuses its normals, and the inverse normal is the dominant cost.
**Where they fail:** the ATM straddle |S_T − K| is nearly even in Z, so ρ > 0 and the variance per path
*rises* by 65 %; only the cheaper draws keep the efficiency at 0.90. Far out of the money, almost no
path pays off in either half of a pair, ρ ≈ 0, and the only gain left is the saved draws (1.5×).

**Control variates** subtract β(X − E[X]) with β estimated on an independent pilot run. The variance
falls by a factor 1 − ρ²_YX. The terminal spot is an excellent control for a deep ITM call, which is
nearly S_T − K (355×). **Where it fails:** for the call struck at 160, the payoff is zero on most paths
where S_T varies most, the correlation collapses, and the variance falls by only 9 % (1.1×). The
geometric-average Asian, whose price is known in closed form, is so close to the arithmetic one
that the variance falls by a factor of 1,300 and the efficiency by more than 1,000 even after paying
for the pilot and the extra logarithms.

**Variance reduction beats parallelism.** On this 4-thread machine, perfect parallel scaling would
give at most 4×. The combined techniques give 41× on the ATM call and over 1,000× on the ITM call and
the Asian, from mathematics alone, and they compose with threading since the engine's result is
independent of the thread count.

*Timing note.* The plain ATM call costs 40 ns per path but the deep ITM and far OTM calls only
26 ns, for the same arithmetic: at the money, the branch in max(0, S_T − K) is a coin flip and
mispredicts half the time. This is why efficiencies are only compared within a payoff.

### 5.3 Randomized quasi-Monte Carlo

Quasi-Monte Carlo replaces pseudo-random points with a low-discrepancy sequence: here Sobol points
(Joe-Kuo direction numbers, checked bit for bit against SciPy) with hash-based Owen scrambling.
Scrambling keeps the net structure that makes the points accurate (tested: every one-dimensional
projection and the first two-dimensional one remain exact nets) while making each point uniformly
distributed, so each scrambling gives an unbiased estimate and independent scramblings give an
error bar. Path payoffs can be built by Brownian bridge, which spends the first, best-distributed
coordinates on the terminal value and then the coarse shape of the path.

![Spread of one estimate against the number of points, pseudo-random against randomized QMC, for four payoffs](figures/qmc_convergence.svg)

*Figure 5 — Standard deviation of a single N-point estimate (over 64 independent seeds or
scramblings) against N, for four problems of increasing difficulty; each line is labeled with its
fitted slope. **What it shows:** pseudo-random Monte Carlo converges at N^{−1/2} everywhere;
randomized QMC reaches N^{−1} on the one-dimensional problems, including the digital, drops to
N^{−0.69} (N^{−0.80} with a bridge) on the 12-dimensional Asian, and barely beats N^{−1/2} on a digital of
the 12-dimensional average. **Why it matters:** QMC is not a free speed-up; its gain is set by the
smoothness and effective dimension of the payoff, and it collapses exactly where a risk system needs
it on exotic, discontinuous payoffs.* Source: [`data/results/qmc_convergence.csv`](../data/results/qmc_convergence.csv).

| Problem | Pseudo-random slope | RQMC slope | RQMC + bridge slope | Spread at N = 2¹⁶: pseudo-random → best RQMC |
|---|---|---|---|---|
| ATM call (1-D, kink) | −0.49 | −1.00 | — | 5.0 × 10⁻² → 1.4 × 10⁻⁴ (**344×**) |
| Digital, K = 130 (1-D, jump) | −0.50 | −1.00 | — | 1.0 × 10⁻³ → 6.6 × 10⁻⁶ (158×) |
| Arithmetic Asian, 12 fixings (12-D, kink) | −0.49 | −0.69 | −0.80 | 3.1 × 10⁻² → 4.2 × 10⁻⁴ (74×) |
| Digital on the 12-fixing average (12-D, jump) | −0.48 | −0.56 | −0.56 | 1.9 × 10⁻³ → 4.1 × 10⁻⁴ (**4.7×**) |

A spread 344 times smaller is worth 344² ≈ 118,000 times more pseudo-random paths. The randomized
QMC means agree with the pseudo-random ones (and with the closed forms where they exist) within
their error bars.

**Mechanism.** In one dimension a kink or a jump sits in a single interval of the stratification, so
only that interval contributes variance and the spread falls as N^{−1}, jump included: the textbook
statement that QMC fails on discontinuous payoffs is wrong in one dimension. In d dimensions a smooth
jump surface such as {average = K}, which is not aligned with the axes, cuts about N^{1−1/d} of the N
elementary boxes of the net, each contributing variance of order N^{−2}. The variance then falls only as
N^{−1−1/d}, a spread slope of −1/2 − 1/(2d) = −0.54 for d = 12, against −0.56 measured. The kinked
Asian sits in between: its difficulty is its effective dimension, which the Brownian bridge lowers
(−0.69 → −0.80, and a 4.3× smaller spread at N = 2¹⁶); for the digital the bridge helps the constant
(2.1×) but not the rate, because the jump surface stays oblique whatever the path construction.

**Remedy.** Use randomized QMC (never unscrambled points, which give no error bar) with a Brownian
bridge for smooth path payoffs; for discontinuous multi-dimensional payoffs, expect little more than
a constant factor and prefer smoothing the payoff (conditional expectation) or plain Monte Carlo with
a control variate.

## 6. Greek estimation under noise *(flagship section)*

*To be written in Phase 5.*

## 7. Risk measures on non-linear positions

*To be written in Phase 6.*

## 8. Model risk beyond GBM

*Phase 7, optional ([`riskengine_research.md`](riskengine_research.md) §12.1).*

## 9. Engineering and performance notes

*To be written in Phase 8.*

## 10. Limitations and future work

*To be written in Phase 8.*

## 11. Conclusion: operational recommendations

*To be written in Phase 8.*

---

## Appendices

*A. derivations of the pathwise and LR estimators; B. parameter tables; C. reference values and
sources; D. reproduction procedure (see §2.2 for the current state).*
