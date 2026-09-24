# RiskEngine-CPP

A C++20 options pricing and risk engine that prices the same option three different ways — closed-form Black-Scholes, a Cox-Ross-Rubinstein binomial tree, and Monte Carlo simulation — and rigorously studies **where and why the methods disagree**.

This is not "a pricer." Black-Scholes, CRR, and Monte Carlo under GBM are three numerical methods converging to the same closed-form answer by construction — that convergence is a correctness gate, not the finding. The actual subject of this project is:

- **Numerical method risk** — how fast and how reliably each method converges, where naive implementations quietly break (e.g. finite-difference Greeks on noisy Monte Carlo estimators), and how to fix it.
- **Model risk**, properly scoped — genuine model risk requires different *assumptions* (stochastic volatility, jumps), not just different numerical schemes. This project is honest about that distinction rather than dressing up numerical-method risk as something bigger.
- **Where risk measures disagree** — parametric (delta-normal) VaR is structurally blind to gamma risk; historical and Monte Carlo VaR are not. Demonstrated on a real short-gamma position, not asserted.

The deliverable is [`docs/model_risk_report.md`](docs/model_risk_report.md) — the code exists to produce that report's numbers, not the other way around. For a guided tour of everything that was built, see [`docs/showcase.md`](docs/showcase.md); for a first run, `build/examples/quickstart` prices one option by every method in the library.

## Why three methods on purpose

| Method | Convergence | What it's good for | What it misses |
|---|---|---|---|
| Black-Scholes (closed-form) | Exact | Ground truth for European options under GBM | Doesn't generalize past GBM/European |
| CRR binomial tree | O(1/n), non-monotonic (oscillates) | American exercise, discrete dividends | Slow, biased at small n unless corrected |
| Monte Carlo | O(1/√N) | Path-dependent payoffs, high-dimensional problems | Noisy; naive Greeks via bumping are actively wrong without common random numbers |

Full methodology, formulas, and verified reference numbers: [`docs/riskengine_research.md`](docs/riskengine_research.md). The core project (numerical-method risk + risk-measure disagreement) is scoped to these three methods; the model-risk phase adds Heston and Merton to demonstrate genuine model risk (different assumptions, not just different numerics).

## Architecture

- Header-only C++20 library (`include/riskengine/`), built around concepts rather than inheritance: `Pricer` (closed form, trees), `StochasticPricer` (Monte Carlo), `PathModel` (GBM, Heston, Merton), terminal and path payoffs, control variates. No virtual dispatch in the hot path; a model or payoff that lacks the required interface is rejected at compile time, and an estimator applied to a Greek it cannot estimate is rejected at run time.
- Every Monte Carlo result is an `Estimate` carrying its standard error, and a stochastic pricer is a pure function of (market, seed): common random numbers come for free, and the fixed-block reduction makes results bit-identical for any thread count.
- Experiments (`experiments/`) produce every number and figure of the report, with their metadata (git commit, compiler, flags, CPU).

## Status

Every phase of the plan ([`docs/riskengine_research.md`](docs/riskengine_research.md) §12) is done, and the report [`docs/model_risk_report.md`](docs/model_risk_report.md) is complete: executive summary, reproducibility protocol, analytic ground truth, trees, Monte Carlo, Greeks under noise (the flagship section), risk measures on non-linear positions, model risk beyond GBM, performance notes, limitations and operational recommendations.

- **Phase 1, analytic ground truth:** Black-Scholes-Merton with continuous dividend yield, closed-form Greeks, bracketed implied-vol solver, digital and discrete geometric-Asian closed forms, [`docs/conventions.md`](docs/conventions.md). The report shows why an in-the-money quote cannot pin its implied vol and the finite-difference V-curve.
- **Phase 2, stochastic infrastructure:** Philox 4×32-10, AS241 inverse normal, Welford/Chan accumulators, fixed-block reduction that is bit-identical for any thread count and across GCC/Clang, and the experiment harness (`experiments/` → `data/results/*.csv` + `.meta.json` → `docs/figures/*.svg`).
- **Phase 3, trees:** CRR, averaged CRR, Leisen-Reimer, BBS and BBS-Richardson, European and American, with delta, gamma and theta read off an extended tree. CRR's error flips sign between even and odd n (−2.00/n against +1.75/n at the money); Leisen-Reimer is the only second-order method (3.6 × 10⁻⁹ at n = 10,000); BBS-Richardson is second order only for n divisible by 4, and early exercise brings every method back to first order.
- **Phase 4, Monte Carlo:** a generic engine over the `PathModel` concept, antithetic and control variates, randomized QMC (Sobol, Owen scrambling, Brownian bridge). Variance reduction reaches 1,500× and QMC a 344× smaller error on an ATM call; the report shows where each technique fails.
- **Phase 5, Greeks under noise:** finite differences (independent seeds, common random numbers), pathwise (checked by forward automatic differentiation), likelihood ratio and mixed estimators, on a call and a digital. Measured convergence rates match theory (e.g. −0.40 against −2/5 for the CRN digital delta); the pathwise digital delta converges, with zero standard error, to 0 instead of 0.0188; the report ends with a payoff × regime × Greek recommendation matrix.
- **Phase 6, risk measures on non-linear positions:** on a delta-hedged short straddle, delta-normal VaR is 0, delta-gamma normal understates the full-revaluation VaR by 41 %, and adding the implied-vol factor doubles it. A 34-year backtest on frozen FRED data (NASDAQ, VIX, T-bill) puts the linear and quadratic methods in the Basel red zone in every 250-day window; historical simulation passes coverage on average but fails Christoffersen's independence test (clustered crisis exceptions). Historical stress scenarios (1987, 2008, 2018, 2020) cost up to 11 times the historical 99 % VaR.
- **Phase 8, performance and write-up:** a Google Benchmark suite (Black-Scholes 34 ns, n = 1,000 tree 0.15 ms, 93 % parallel efficiency on 4 threads), two performance pathologies found and measured (subnormal numbers slowing the tree 9×, fixed with bit-identical results; false sharing 20× on adjacent atomics), and one command, `cmake --build <dir> --target report`, that regenerates every result and figure.
- **Phase 7, model risk:** Heston (little-trap characteristic function; Andersen's QE and full-truncation Euler) and Merton (series; exact simulation), checked against 30-digit references. Calibrated to the same at-the-money price as Black-Scholes, they disagree by up to 146 % on an 80 put and 108 % on an up-and-out barrier; a Black-Scholes delta hedge stalls at a P&L sd of 2.8–2.9 per option in their worlds, whatever the frequency; under Heston with Feller violated, Euler at 320 steps is still 27 standard errors off where QE is within noise.

## Building

```
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

Requires a C++20 compiler (GCC, Clang, or MSVC) and CMake 3.25+.

To regenerate the report's results and figures (Release build of a clean tree; Python tools need `pip install -r tools/requirements.txt`):

```
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --target report   # every experiment, then every figure (~5 min)
```

Everything except the timing columns reproduces bit for bit. Benchmarks: see [`bench/README.md`](bench/README.md).

## License

MIT — see [LICENSE](LICENSE).
