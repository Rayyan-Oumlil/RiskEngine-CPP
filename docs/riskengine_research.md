# RiskEngine-CPP — Technical Development Plan

**Sourcing note:** Sections 1–5 and 8 were verified by *executing* every quantitative claim in C++20/Python on a 12-core machine, not by citing blog posts. Sections 6 (historical data sources) and 7 (recruiter/interview signal) were flagged by the research agent as needing independent web verification — see the follow-up addendum at the bottom of this file for that pass.

---

## 1. Mathematical correctness deep-dive

### 1.1 Black-Scholes closed form

With spot S, strike K, rate r, continuous dividend yield q, vol σ, maturity T:

```
d1 = [ln(S/K) + (r - q + 0.5*σ²)T] / (σ√T)
d2 = d1 - σ√T

C = S·e^(-qT)·N(d1) - K·e^(-rT)·N(d2)
P = K·e^(-rT)·N(-d2) - S·e^(-qT)·N(-d1)
```

**Canonical test vector** (memorize — interviewers use it): S=K=100, r=5%, q=0, σ=20%, T=1 gives d1=0.35, d2=0.15, **call = 10.450584**, **put = 5.573526**. C++ and Python implementations agreed to 1e-9.

**Implementation pitfalls (all real, all bite):**

| Pitfall | Why it bites | Fix |
|---|---|---|
| Hand-rolled `norm_cdf` (Abramowitz-Stegun 7.1.26) | ~7.5e-8 absolute error — swamps convergence tests | `0.5*std::erfc(-x*M_SQRT1_2)`; erfc is correctly rounded |
| Computing put via `C - S + K*exp(-rT)` | Hides sign bugs | Implement both directly, assert put-call parity as a test |
| T→0 or σ→0 | Division by zero in d1 | Branch to intrinsic value max(0,S-K) below a tolerance |
| Deep ITM/OTM | N(d1) underflows, catastrophic cancellation | Test at σ√T > 5 and moneyness 0.1/10 |
| Passing σ² where σ expected | Silent, plausible-looking wrong prices | Strong typedefs (`struct Vol{double v;}`) |

**Invariants that belong in your test suite** (not just "price ≈ expected"):
- Put-call parity: C - P = S·e^(-qT) - K·e^(-rT) — must hold to ~1e-14
- Monotonicity: ∂C/∂K < 0, ∂C/∂σ > 0, ∂C/∂S > 0
- Bounds: max(0, S·e^(-qT) - K·e^(-rT)) ≤ C ≤ S·e^(-qT)
- Convexity in strike: C(K1) - 2C(K2) + C(K3) ≥ 0 for equally-spaced strikes (butterfly ≥ 0, i.e. no static arbitrage)

That last one signals you've actually read the literature rather than copied a formula.

### 1.2 Cox-Ross-Rubinstein binomial tree

```
Δt = T/n
u = e^(σ√Δt), d = 1/u
p = [e^((r-q)Δt) - d] / (u - d)
```

Backward induction: V[j][i] = e^(-rΔt)·(p·V[j+1][i+1] + (1-p)·V[j][i+1]), and for American: V[j][i] = max(V[j][i], intrinsic(S[j][i])).

**The convergence result nearly every student project gets wrong:** CRR converges O(1/n), but **not monotonically** — it oscillates. Measured at ATM, consecutive n:

```
n=100 err=-1.997e-02    n=101 err=+1.737e-02
n=102 err=-1.958e-02    n=103 err=+1.703e-02
n=104 err=-1.920e-02    n=105 err=+1.671e-02
```

Error flips sign at every step. Sampling only even n (as most students do — 10, 50, 100, 200, 400…) gives a beautiful clean error-halving ratio of exactly 2.000, and a false conclusion of "O(1/n), confirmed." That's an artifact of the sampling grid, not the truth.

Caused by the strike falling between terminal nodes — the "distribution error" analyzed in **Leisen & Reimer (1996), "Binomial models for option valuation – examining and improving convergence," Applied Mathematical Finance 3(4)**. Three mitigations, increasing sophistication:

1. **Averaging**: 0.5·(Vn + Vn+1) — cancels the leading oscillatory term, cheap
2. **Leisen-Reimer tree**: choose u,d,p via Peizer-Pratt inversion, centered on the strike; converges O(1/n²) monotonically with odd n
3. **Black-Scholes smoothing / Broadie-Detemple**: replace the final timestep with the closed-form BS value at each penultimate node

Implementing (1) as core and (2) as stretch is a genuinely strong differentiator.

**Verified invariants at n=20,000:**
- Non-dividend American call == European call: diff = **0.00e+00 exactly** (Merton's theorem — never optimal to early-exercise an American call on a non-dividend stock). Free, exact unit test and a classic interview question.
- American put (6.090333) > European put (5.573426), early-exercise premium 0.516907 > 0

**Parameter validity trap:** p ∈ (0,1) requires Δt < σ²/(r-q)². At high rates/low vol/few steps, p escapes [0,1] and produces nonsense prices. Assert it.

### 1.3 Monte Carlo

For European payoffs, sample the terminal value directly — do NOT loop over time steps:

```
S_T = S0 * exp((r - q - 0.5σ²)T + σ√T·Z),   Z ~ N(0,1)

price_hat = e^(-rT) * mean(max(0, S_T[i] - K))
SE_hat = stddev / sqrt(N)
```

Time-stepping a GBM for a European option adds discretization bias for zero benefit — only step when the payoff is path-dependent.

**Convergence is O(1/√N)** — stderr halves per 4× paths. Measured:

```
n=  10,000  price=10.30949  stderr=0.14643
n=  40,000  price=10.44142  stderr=0.07329
n= 160,000  price=10.48901  stderr=0.03687
n= 640,000  price=10.45292  stderr=0.01841
```

Price doesn't decrease monotonically in error — that's correct and expected, it's a random variable. Report a confidence interval, not a point estimate.

**Variance reduction, measured at 100k draws (equal RNG cost):**

| Method | stderr | Speedup (equal accuracy) |
|---|---|---|
| Plain | 0.046512 | 1× |
| Antithetic | 0.023168 | 4.0× |
| Control variate (S_T) | 0.017697 | 6.9× |
| **Antithetic + CV** | **0.006088** | **58×** |

58× fewer paths for the same accuracy. **Variance reduction beats parallelism** — 58× from math vs. 6× from 8 threads. Worth saying explicitly in an interview.

Control variate uses E[S_T] = S0·e^((r-q)T) known analytically, optimal coefficient b* = Cov(P, S_T)/Var(S_T) estimated from the sample. Reference: **Glasserman, "Monte Carlo Methods in Financial Engineering" (2003), Ch. 4**.

### 1.4 The cross-model convergence test (the project's thesis)

- **Binomial → BS**: plot log|err| vs log(n), fit slope. Expect −1.0 for CRR, −2.0 for Leisen-Reimer. Report the fitted slope, sample consecutive n so oscillation is visible.
- **MC → BS**: plot log(stderr) vs log(N), fit slope. Expect −0.5. Separately verify empirical coverage of the 95% CI: run 1000 independent MC pricings, check ~950 intervals contain the true BS price. Coverage testing is what a real risk quant does and almost no student project does.
- **Statistical test, not eyeballing**: assert |MC - BS| < 3×SE. Fails ~0.3% of the time by construction — seed the RNG in CI to make it deterministic.

---

## 2. Greeks

### 2.1 Analytic Black-Scholes (verified against implementation)

For a call with q=0, at canonical parameters:

| Greek | Formula | Value |
|---|---|---|
| Delta | e^(-qT)·N(d1) | 0.636831 |
| Gamma | e^(-qT)·φ(d1) / (S·σ·√T) | 0.018762 |
| Vega | S·e^(-qT)·φ(d1)·√T | 37.524035 |
| Theta | −S·e^(-qT)·φ(d1)·σ/(2√T) − r·K·e^(-rT)·N(d2) + q·S·e^(-qT)·N(d1) | −6.414028 /yr |
| Rho | K·T·e^(-rT)·N(d2) | 53.232482 |

**Reporting-convention traps** (cause more confusion than the math):
- **Vega** is per 1.00 of vol (100 vol points). Desks quote per 1 vol point: divide by 100 → **0.375**.
- **Theta** above is per year. Desks quote per calendar day: divide by 365 → **−0.01757** (some use 252 trading days — state which).
- **Rho** is per 1.00 of rate; per bp divide by 10,000.
- Gamma and Vega are identical for calls/puts; Delta, Theta, Rho are not.

State conventions explicitly in the README.

### 2.2 Finite-difference Greeks — the single most important result in this project

Central difference: Δ ≈ [V(S+h) - V(S-h)] / 2h, Γ ≈ [V(S+h) - 2V(S) + V(S-h)] / h²

For **deterministic** pricers (BS, binomial) this just works — FD Greeks matched analytic to 8.6e-09 (delta), 2.7e-10 (gamma). Truncation error O(h²) for central vs O(h) for forward — use central. Optimal h balances truncation against floating-point cancellation: h* ≈ S·ε^(1/3) ≈ 6e-6·S for gamma. Use a **relative** bump (h = 1e-4·S), never absolute.

**For Monte Carlo, naive finite differences catastrophically fail.** This is the result to lead with. Measured, 200k paths, true delta = 0.636831:

| h | CRN (same seed) | Independent seeds |
|---|---|---|
| 5.00 | 0.63286 (err −0.004) | 0.6405 (err +0.004) |
| 1.00 | 0.63442 (err −0.002) | 0.6648 (err +0.028) |
| 0.50 | 0.63442 (err −0.002) | 0.6924 (err +0.056) |
| 0.10 | 0.63447 (err −0.002) | 0.9124 (err +0.276) |
| **0.01** | **0.63449 (err −0.002)** | **3.3874 (err +2.751)** |

Without **common random numbers**, shrinking h makes the estimate *worse* — error explodes to 275% because MC noise O(σ_MC/h) dominates as h→0. Inverts standard textbook FD intuition. With CRN (reuse identical random draws for both bumps), noise cancels and the estimate is stable across four orders of magnitude of h.

Implementation: CRN means re-seeding the generator identically, not storing paths. A pricer interface taking an explicit seed makes this natural.

### 2.3 Beyond bumping: pathwise and likelihood-ratio (strong stretch)

**Pathwise (IPA)** — differentiate the payoff, unbiased, single pricing run:

```
Δ_pathwise = e^(-rT) * E[1{S_T > K} * S_T/S0]
```

Measured: **0.634511** (err −0.0023) in one pass, no bump, no second pricing. ~2× cheaper than central-difference bumping *and* lower variance.

**Pathwise fails for Gamma** — the second derivative of the payoff is a Dirac delta. Classic interview probe ("why can't you just differentiate twice?"). Fix: **likelihood-ratio method**, differentiating the density instead of the payoff:

```
Γ_LR = e^(-rT) * E[max(0, S_T-K) * (Z² - Z·σ√T - 1) / (S0²·σ²·T)]
```

Measured: **0.018554** vs true 0.018762 (err −0.0002). Both verified working.

Knowing when pathwise breaks and why LR is the remedy is genuinely differentiating knowledge — Glasserman Ch. 7.

---

## 3. VaR methodology

### 3.1 Three methods, and why the distinction matters

| Method | How | Captures fat tails? | Captures optionality? |
|---|---|---|---|
| Parametric (variance-covariance) | VaR_α = -(μ + z_α·σ)·V | No (Gaussian assumption) | No (delta-normal is linear) |
| Historical simulation | Reprice under actual past return vectors, take empirical quantile | Yes, but only tails that occurred | Yes, if full revaluation |
| Monte Carlo | Simulate risk factors, full-revalue, take quantile | Depends on chosen distribution | Yes |

For an **options** book: delta-normal parametric VaR is structurally blind to gamma. A short-straddle position has near-zero delta and catastrophic risk; parametric VaR reports ≈0. **Computing all three on the same option and showing parametric VaR miss the risk that MC VaR catches is a compelling result** — ties directly to the model-risk thesis.

### 3.2 Conventions students get wrong

- **Sign**: VaR reported as a positive loss number. VaR_99% = $1M means "1% chance of losing more than $1M."
- **Confidence**: 99% VaR = 1st percentile of P&L = 99th percentile of loss. Off-by-one-tail is the most common bug.
- **Horizon**: 1-day vs 10-day. The √10 scaling rule assumes i.i.d. returns — **false** under volatility clustering, understates risk in crises.
- **Quantile estimator**: with N sims and α=0.99, index ⌈αN⌉ on sorted losses.
- **Use `std::nth_element`**, not a full sort — O(N) vs O(N log N), and it's exactly the right tool since you need the quantile plus the tail beyond it.

### 3.3 Include CVaR/Expected Shortfall — not optional

CVaR/ES = mean loss conditional on exceeding VaR: ES_α = E[L | L > VaR_α]

**Three reasons it's non-optional:**

**(a) VaR is not subadditive — demonstrated failure.** Two independent bonds, 4% default probability, loss 100 each, 2M simulations:

```
VaR95(A) = 0    VaR95(B) = 0    sum = 0
VaR95(A+B) = 100            -->  subadditive? False
```

Diversifying *increased* measured VaR from 0 to 100. VaR fails **Artzner, Delbaen, Eber, Heath (1999), "Coherent Measures of Risk," Mathematical Finance 9(3)** — the paper defining the four coherence axioms (monotonicity, translation invariance, positive homogeneity, subadditivity). ES satisfies all four. Put this 3-line counterexample directly in the README.

**(b) VaR says nothing about tail severity.** Two books with identical VaR can have wildly different beyond-VaR losses.

**(c) Regulators moved.** Basel III/FRTB (**BCBS d457, "Minimum capital requirements for market risk," 2019**) replaced 99% VaR with **97.5% Expected Shortfall** for trading-book capital. Under normality, ES@97.5% = 2.3378 ≈ VaR@99% = 2.3263 — nearly identical, so the switch preserved capital levels while gaining coherence and tail sensitivity.

**Reference values (standard normal losses, analytic ES_α = φ(z_α)/(1-α)):**

| α | VaR | ES |
|---|---|---|
| 95% | 1.6449 | 2.0627 |
| 97.5% | 1.9600 | 2.3378 |
| 99% | 2.3263 | 2.6652 |

Empirical estimator converges to these (N=1e6: VaR 2.3257, ES 2.6642). VaR converges *slowly and noisily* — quantile estimator variance ∝ α(1-α)/(N·f(q)²), and f(q) is tiny in the tail. At N=1000, measured VaR error was 5.6%. **Report a bootstrap confidence interval on your VaR number.**

### 3.4 Backtesting VaR

A VaR number nobody validated is worthless. Implement **Kupiec's POF test (1995)** — unconditional coverage, likelihood-ratio on exception count — and ideally **Christoffersen's (1998)** conditional-coverage test, which checks exceptions aren't *clustered* (exactly what happens in a crisis, and what √t-scaling misses). Basel's traffic-light system (green/yellow/red at 4/5–9/10+ exceptions in 250 days) is a 10-line implementation, instantly recognizable to a risk interviewer.

---

## 4. C++ engineering practice for this domain

### 4.1 QuantLib — what to take, what to leave

QuantLib's core abstraction: `Instrument` → `PricingEngine`, with `Instrument::performCalculations()` delegating to a pluggable engine, wired via Observer/Observable so market-data changes invalidate cached results. `VanillaOption` can be priced by `AnalyticEuropeanEngine`, `BinomialVanillaEngine<CoxRossRubinstein>`, or `MCEuropeanEngine<PseudoRandom>` — same instrument, swappable engines.

**Take**: the instrument/engine separation — precisely the "compare three models on one option" requirement.

**Leave**: QuantLib is ~1990s-style C++ — `shared_ptr` everywhere, deep virtual hierarchies, `Handle<>` double-indirection, heavy runtime polymorphism in hot loops. Do not imitate that in 2026. Knowing QuantLib's design *and* articulating why you'd do it differently with C++20 is a stronger position than either alone.

### 4.2 Recommended architecture (compiled and verified)

**Concepts + templates for the hot path, `std::variant` for the comparison layer.**

```cpp
template <typename P>
concept Pricer = requires(const P p, const OptionSpec& o, const MarketData& m) {
    { p.price(o, m) } -> std::same_as<PriceResult>;
};

template <Pricer P>   // one generic implementation, works for ALL pricers
Greeks fd_greeks(const P& p, const OptionSpec& o, MarketData m, double hS, double hv);
```

`fd_greeks` written **once**, works for Black-Scholes, binomial, and Monte Carlo with zero inheritance and zero virtual dispatch — fully inlined. `static_assert(Pricer<BlackScholes>)` gives a readable compile-time error instead of template vomit.

For the runtime comparison layer (iterate over "all three models"):

```cpp
using AnyPricer = std::variant<BlackScholes, BinomialCRR, MonteCarlo>;
std::vector<AnyPricer> models = {...};
for (auto& m : models)
    std::visit([&](auto& p){ results.push_back(p.price(opt, mkt)); }, m);
```

Idiomatic C++17/20 answer for a **closed set of types**: `std::variant` + `std::visit`, not an abstract base class. Value-semantic (no heap, no shared_ptr), exhaustive at compile time, cache-friendly. Being able to explain *why variant over inheritance here* (closed set, value semantics, devirtualization, no allocation) is a real C++-design signal.

**Other modern features that genuinely earn their place:**

| Feature | Use | Why |
|---|---|---|
| `std::span<const double>` | Path/price buffers | Non-owning, no allocation |
| `std::mdspan` (C++23) | Binomial lattice views | Clean 2-D indexing over flat storage |
| `constexpr`/`consteval` | `norm_cdf`, node counts | Compile-time test vectors |
| Strong typedefs | `Vol`, `Rate`, `Strike` | Kills σ-vs-σ² bug class |
| `std::expected` (C++23) | Validation failures | Error handling without exceptions in hot paths |
| `[[nodiscard]]` | All pricing fns | Discarding a price is always a bug |
| `std::format` | Reporting | Type-safe, fast |

**Avoid:** virtual calls inside the Monte Carlo inner loop (kills inlining/vectorization), `shared_ptr` for market data (atomic refcount contention across threads), exceptions in hot paths.

### 4.3 Repo structure

```
RiskEngine-CPP/
├── CMakeLists.txt              # C++20, warnings-as-errors, sanitizer presets
├── include/riskengine/
│   ├── types.hpp               # MarketData, OptionSpec, Greeks, strong typedefs
│   ├── concepts.hpp            # Pricer concept
│   ├── pricers/{black_scholes,binomial,monte_carlo}.hpp
│   ├── greeks/{analytic,finite_difference,pathwise}.hpp
│   ├── risk/{var,expected_shortfall,backtest}.hpp
│   └── stress/scenario.hpp
├── src/                        # non-header-only impls
├── tests/                      # Catch2 or GoogleTest
│   ├── test_parity.cpp         # put-call parity, bounds, convexity
│   ├── test_convergence.cpp    # slope fitting, CI coverage
│   ├── test_greeks.cpp         # analytic vs FD vs pathwise
│   └── test_var.cpp            # analytic normal VaR/ES, subadditivity
├── bench/                      # Google Benchmark
├── data/                       # historical CSVs
└── docs/model_risk_report.md   # THE deliverable
```

Non-negotiables for the "well-engineered" claim: CMake presets, CI (GitHub Actions) building GCC+Clang+MSVC, `-Wall -Wextra -Werror`, ASan/UBSan job, clang-format + clang-tidy, codecov badge. One evening of work, the difference between "student project" and "engineered."

---

## 5. Performance engineering

### 5.1 Parallel Monte Carlo — measured on a 12-core machine

20M paths, antithetic, -O2 -march=native:

| Threads | Time | Throughput | Speedup | stderr |
|---|---|---|---|---|
| 1 | 346 ms | 57.8 M/s | 1.00× | 0.001644 |
| 2 | 172 ms | 116.3 M/s | 2.01× | 0.001644 |
| 4 | 96 ms | 208.3 M/s | 3.60× | 0.001643 |
| 8 | 57 ms | 350.9 M/s | **6.07×** | 0.001643 |

**stderr is invariant across thread counts** — proves the parallel RNG design is sound. Monte Carlo is embarrassingly parallel *only if the RNG is right*; naive approaches (shared mutex-locked generator, or same seed per thread) either serialize or produce correlated streams that silently break error estimates.

Design: each thread owns an independent `std::mt19937_64` seeded by a counter × golden-ratio constant (0x9E3779B97F4A7C15). Better for reproducibility: **fixed block decomposition** — path i always uses the same random draws regardless of thread count, so results are bit-identical whether run on 1 or 8 threads. Genuinely production-grade (reproducible risk numbers are an audit requirement).

**Tradeoffs:**

| Approach | Verdict |
|---|---|
| `std::thread` + block decomposition | **Use this.** Explicit, portable, full RNG control |
| `std::async` | Launch policy is implementation-defined; can silently run serially |
| OpenMP | Great speedup, one pragma, awkward with per-thread RNG state |
| `std::execution::par` | Elegant but libstdc++ needs TBB linked; RNG state per element awkward |
| `std::jthread` (C++20) | Use over `std::thread` — auto-joins, stop tokens |

Why 6.07× not 8×: hyperthreading, turbo clock-down under all-core load, memory bandwidth. **Explaining the gap is more impressive than claiming linear scaling.**

### 5.2 False sharing — measured, 9.7× slowdown

Best "I profiled and fixed it" story available in this project, and it's real:

| Threads | Naive (accumulators adjacent) | alignas(64) padded | Speedup |
|---|---|---|---|
| 4 | 426 ms | 44 ms | **9.7×** |
| 8 | 462 ms | 50 ms | **9.2×** |

`std::hardware_destructive_interference_size` = 64 on this machine. Per-thread accumulators sharing a cache line means every write invalidates the line in every other core's L1 (MESI ping-pong). Padding to a cache line each fixes it. 4→8 threads barely helps in the naive case — the bus is saturated, not the cores.

Best practice: accumulate in a thread-local stack variable, write once at the end (avoids the issue entirely); `alignas(64)` is the fix when a shared array is unavoidable. Show both.

### 5.3 SIMD

Real opportunities, in order of effort/reward:
1. **Auto-vectorization first.** `-O3 -march=native -ffast-math`, check with `-fopt-info-vec`. `std::exp` doesn't vectorize by default; SLEEF/Vc gives vectorized transcendentals. `exp` dominates the MC inner loop.
2. **Batch the RNG.** Generate normals in blocks into a contiguous buffer, then a vectorizable payoff loop. Structure-of-Arrays, not Array-of-Structures.
3. **`std::experimental::simd`** or explicit intrinsics — expect maybe 2–4× over scalar on the payoff loop.
4. **Ziggurat over Box-Muller**, or note `std::normal_distribution` is stateful (caches a spare) — a vectorization barrier.

SIMD is a stretch goal — don't claim it unless benchmarked.

### 5.4 Binomial tree memory — measured

| n | Full triangle | Rolling 1-D |
|---|---|---|
| 10,000 | **400.0 MB** | **0.080 MB** |

5000× reduction. Storing the whole lattice is O(n²) and pointless for backward induction — only one timestep needed at a time. Timings with rolling array: n=1000 → 4ms, 5000 → 108ms, 10000 → 435ms, 20000 → 1734ms (quadratic time scaling is inherent — O(n²) node updates — separate from the O(n) space win).

Further wins: precompute u^j powers instead of calling `std::pow` in the inner loop (easy 2–3×); contiguous iteration for prefetch; float vs double tradeoff experiment.

### 5.5 Targets to report

| Metric | Target |
|---|---|
| BS single price | < 100 ns |
| Binomial n=1000 European | < 5 ms |
| MC 1M paths, 1 thread | < 20 ms |
| MC 20M paths, 8 threads | < 60 ms (measured 57 ms) |
| Parallel efficiency, 8 threads | > 70% (measured 76%) |
| Full Greeks (5) via CRN bumping | < 10× single price |

Use Google Benchmark, report median with error bars, state CPU model — a benchmark table without hardware context is noise.

---

## 6. Historical stress-test data — SEE ADDENDUM BELOW for web-verified update

## 7. What separates strong from mediocre — SEE ADDENDUM BELOW for web-verified update

---

## 8. Sequenced build plan

Dependency-ordered. Each phase has a validation gate — don't proceed until it passes.

**Phase 0 — Scaffold (½ week):** CMake 3.25+ presets, C++20, `-Wall -Wextra -Werror`, Catch2/GoogleTest via FetchContent, GitHub Actions (GCC+Clang+MSVC), ASan/UBSan job, clang-format/tidy. `types.hpp`. **Gate:** CI green on 3 compilers, one trivial test passing.

**Phase 1 — Black-Scholes + the concept (1 week):** norm_cdf via erfc, call/put, Pricer concept (static_assert'ed), analytic Greeks with documented conventions. **Gate:** price = 10.450584/5.573526 to 1e-9; parity to 1e-14; all 5 Greeks match reference; edge cases don't NaN. *Everything depends on this.*

**Phase 2 — Binomial tree (1 week):** CRR with rolling 1-D array (not triangle). European + American. Convergence study. **Gate:** converges to BS with fitted log-log slope ≈ −1.0; oscillation visible sampling consecutive n; American put > European put; non-dividend American call == European to 1e-12; p∈(0,1) asserted. *Stretch: Leisen-Reimer, slope ≈ −2.0.*

**Phase 3 — Monte Carlo (1–1.5 weeks):** Terminal-value sampling. Per-instance seeded RNG. Standard error. Antithetic, then control variate. **Gate:** |MC-BS| < 3×SE; stderr slope ≈ −0.5; 95% CI empirical coverage ≈95% over 1000 runs; variance-reduction speedups tabulated. *Most projects stop here — you're a third of the way.*

**Phase 4 — Unified Greeks (1 week):** Generic fd_greeks over Pricer concept. Common random numbers for MC. Relative bump sizing. h-sensitivity study with/without CRN. **Gate:** FD Greeks match analytic to ~1e-6 for BS/binomial; MC delta stable across h∈[0.01,5] with CRN, demonstrably unstable without. *Stretch: pathwise delta, LR gamma.*

**Phase 5 — VaR + ES (1–1.5 weeks):** P&L distribution from MC engine. Parametric, historical, MC VaR. ES alongside. nth_element. Bootstrap CIs. Kupiec POF backtest. **Gate:** matches analytic normal table; subadditivity counterexample reproduced; delta-normal VaR demonstrably understates risk for a gamma-heavy position (headline result).

**Phase 6 — Historical stress testing (1.5 weeks):** Data loaders. Scenario definitions (2008/2020/1987/2018). Joint factor shocks (spot+vol+rate together, not one at a time). Run all three models per scenario. Delta-hedged P&L. **Gate:** real data loads/sanity-checked; scenarios reproduce known market moves; quantified statement of where models diverge most and why, with error bars. *This phase is the project's thesis — don't rush it.*

**Phase 7 — Report + polish (1 week):** docs/model_risk_report.md, plots, benchmark table with CPU, limitations section, README leading with the finding. **Gate:** someone who doesn't know options can read the README and understand the finding.

**Core: ~8–9 weeks part-time.** Then stretch, ranked:

| Stretch | Value | Effort |
|---|---|---|
| Parallel MC + false-sharing fix | High — real, measurable C++ signal | Low |
| Google Benchmark suite | High — makes perf claims credible | Low |
| Leisen-Reimer tree | High — literature depth | Low |
| Pathwise/LR Greeks | High — filters most candidates | Medium |
| Implied vol solver (Newton + Brent fallback) | High — very commonly asked | Low |
| Longstaff-Schwartz (American MC) | High — genuinely hard, real technique | High |
| **Heston / jump-diffusion as 4th model** | **Very high for the model-risk thesis** | High |
| SIMD | Medium — only if benchmarked | High |
| Real-time streaming | Medium — nice demo, less quant signal | Medium |
| Python bindings (pybind11/nanobind) | Medium-high — matches real desk workflow | Low |

**If only three:** parallel MC + benchmarks (cheap, high signal, closes the C++ gap), implied vol solver (cheap, universally asked), Heston (expensive but makes "model risk" literally true rather than "numerical method risk" dressed up).

**The honest framing to lead with:** Black-Scholes, CRR, and Monte Carlo under GBM are all *the same model*, three numerical methods — they converge to the same number by construction (that's Phase 2/3's gate). Genuine **model** risk needs different *assumptions* (stochastic vol, jumps). Phrase the core project as "numerical method risk + assumption stress-testing," and Heston as the upgrade to true model risk. Pre-empting this in the README converts the project's biggest vulnerability into a demonstration of understanding.

---

## Key references

- Hull, *Options, Futures, and Other Derivatives* — Ch. 13 (binomial), 15 (BS), 19 (Greeks), 21 (numerical procedures), 22 (VaR/ES)
- Glasserman, *Monte Carlo Methods in Financial Engineering* (2003) — Ch. 4 (variance reduction), Ch. 7 (Greeks)
- Joshi, *C++ Design Patterns and Derivatives Pricing* — standard text for this project's architecture question
- Artzner, Delbaen, Eber, Heath (1999), *Coherent Measures of Risk*, Math. Finance 9(3)
- Leisen & Reimer (1996), *Binomial models for option valuation*, Applied Math. Finance 3(4)
- Longstaff & Schwartz (2001), *Valuing American Options by Simulation*, RFS 14(1)
- BCBS d457 (2019), *Minimum capital requirements for market risk* (FRTB ES@97.5%)
- Kupiec (1995), POF test; Christoffersen (1998), conditional coverage
- QuantLib source (`ql/instrument.hpp`, `ql/pricingengine.hpp`) — read the engine abstraction, don't copy the 1990s C++ style

---

# ADDENDUM — Web-verified update to Sections 6 & 7

*(Added after a fresh web search pass, since the original research agent had no live web access this session.)*

## 6 (updated) — Historical stress-test data sources, web-confirmed

All three sources the agent named check out as real, current (2026), and free:

| Source | What | Cost | Notes |
|---|---|---|---|
| **Yahoo Finance** (`yfinance` Python package) | Daily OHLCV, decades of history for most US equities, indices, ETFs, up to 1-min granularity | Free | Most widely used free source; unofficial API, occasionally breaks — has an actively maintained Python wrapper |
| **Stooq** | Ticker/date/OHLCV, direct CSV export | Free, no registration | Reliable bulk downloads, global coverage — good fallback if yfinance breaks |
| **FRED** (Federal Reserve Bank of St. Louis) | VIXCLS, SP500, Treasury rates (DGS3MO/DGS1/DGS10), hundreds of thousands of official economic series | Free, official | Requires a free API key (register on site); the `fredapi` Python package wraps it. **Best source for the real risk-free rate `r`** per historical date, not a hardcoded 5% |
| **EODHD (new find)** | Daily OHLCV for stocks/ETFs/indices/forex/crypto — **S&P 500 back to 1927 in one REST call**, CSV or JSON | Free tier available | Worth checking for a single clean pull covering Black Monday 1987 through COVID 2020 in one request, rather than stitching Yahoo + Stooq |

**Practical recommendation:** use FRED for the risk-free rate and VIX (official source, cleanest), and EODHD or yfinance for the underlying price history — check EODHD's free-tier rate limits before committing, since the original agent's advice to independently verify still stands for anything usage-limit-related.

The original agent's core caution stands and is worth repeating: **free historical *option chain* data for 2008 is essentially unobtainable.** Don't fake it — the honest, correct design is real historical underlying prices + VIX as the implied-vol proxy, pricing a *hypothetical* option under those historical conditions, stated explicitly as a limitation in the README.

## 7 (updated) — What makes a quant portfolio project stand out, web-confirmed

Web search independently confirms the agent's framing, with a few concrete additions:

- **"Build 2-3 projects in depth" beats a pile of shallow ones** — matches the plan of shipping RiskEngine-CPP as one deep, well-documented project rather than many superficial ones.
- **The four canonical high-impact project types employers look for**: pricing engines, risk management systems (VaR + position sizing), algorithmic trading/backtesting, and market microstructure. RiskEngine-CPP sits at the intersection of the first two — pricing engine *and* risk system in one project, which is a stronger combination than either alone.
- **"Document your process, include clear code and mathematics, be ready to discuss your insights"** — directly validates the agent's advice to lead the README with the *finding* (a written report), not just the code.
- **A concrete example of what "elite tier" looks like**: one real public repo found in this search (AshJha0/quant-portfolio) builds 10 project areas *twice* — once in Python, once with C++/Rust performance twins for the 4 highest-value pricing/risk engines — specifically to demonstrate depth *and* performance optimization side by side. This is a useful calibration point: it's more scope than you need for one strong project, but it confirms the "same model, multiple implementations, compare rigorously" structure (which is exactly RiskEngine-CPP's core thesis) is recognized as differentiating, not redundant.

**Net effect on the plan:** no changes needed — the original build plan (Sections 1–5, 8) is independently validated by this web pass, not contradicted. The one actionable addition is: check EODHD as a possibly cleaner single-source alternative to stitching Yahoo Finance + Stooq for the historical price series in Phase 6.

**Sources for this addendum:**
- [Stooq Historical Stock Prices Scraper API — Apify](https://apify.com/parseforge/stooq-historical-stocks-scraper/api/python)
- [Reliably download historical market data with Python — Ran Aroussi](https://aroussi.com/post/python-yahoo-finance)
- [End-of-Day (EOD) Historical Stock Market Data API — EODHD](https://eodhd.com/financial-apis/api-for-historical-data-and-volumes)
- [Historical Stock Price Data: Where to Find Free OHLCV Data — FinGrab](https://fingrab.app/blog/historical-stock-price-data)
- [Top 10 Projects For Quantitative Finance Roles — Dataloopr](https://dataloopr.com/blog/top-10-projects-for-quantitative-finance-roles-121/)
- [Quantitative Finance Portfolio Projects — OpenQuant](https://openquant.co/blog/quantitative-finance-portfolio-projects)
- [AshJha0/quant-portfolio — GitHub](https://github.com/AshJha0/quant-portfolio)
