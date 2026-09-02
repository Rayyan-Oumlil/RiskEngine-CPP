# RiskEngine-CPP

A C++20 options pricing and risk engine that prices the same option three different ways — closed-form Black-Scholes, a Cox-Ross-Rubinstein binomial tree, and Monte Carlo simulation — and rigorously studies **where and why the methods disagree**.

This is not "a pricer." Black-Scholes, CRR, and Monte Carlo under GBM are three numerical methods converging to the same closed-form answer by construction — that convergence is a correctness gate, not the finding. The actual subject of this project is:

- **Numerical method risk** — how fast and how reliably each method converges, where naive implementations quietly break (e.g. finite-difference Greeks on noisy Monte Carlo estimators), and how to fix it.
- **Model risk**, properly scoped — genuine model risk requires different *assumptions* (stochastic volatility, jumps), not just different numerical schemes. This project is honest about that distinction rather than dressing up numerical-method risk as something bigger.
- **Where risk measures disagree** — parametric (delta-normal) VaR is structurally blind to gamma risk; historical and Monte Carlo VaR are not. Demonstrated on a real short-gamma position, not asserted.

The deliverable is [`docs/model_risk_report.md`](docs/model_risk_report.md) (written up as each phase lands) — the code exists to produce that report's numbers, not the other way around.

## Why three methods on purpose

| Method | Convergence | What it's good for | What it misses |
|---|---|---|---|
| Black-Scholes (closed-form) | Exact | Ground truth for European options under GBM | Doesn't generalize past GBM/European |
| CRR binomial tree | O(1/n), non-monotonic (oscillates) | American exercise, discrete dividends | Slow, biased at small n unless corrected |
| Monte Carlo | O(1/√N) | Path-dependent payoffs, high-dimensional problems | Noisy; naive Greeks via bumping are actively wrong without common random numbers |

Full methodology, formulas, and verified reference numbers: [`docs/riskengine_research.md`](docs/riskengine_research.md). The core project (numerical-method risk + risk-measure disagreement) is scoped to these three methods; an optional extension phase adds Heston and Merton to demonstrate genuine model risk (different assumptions, not just different numerics) — see the plan's descoping section for what's core vs. stretch.

## Architecture

- A C++20 `Pricer` concept means Greeks-computation code (finite-difference, pathwise) is written **once** and works generically across all three pricing models — no inheritance, no virtual dispatch in the hot path.
- `std::variant<BlackScholes, BinomialCRR, MonteCarlo>` + `std::visit` for the runtime comparison layer, since the set of pricing models is closed and fixed. Value-semantic, no heap allocation, exhaustive at compile time.

## Status

Phases 1, 2, 4 and 5 of the plan ([`docs/riskengine_research.md`](docs/riskengine_research.md) §12) are done. The report [`docs/model_risk_report.md`](docs/model_risk_report.md) has §2 (reproducibility protocol), §3 (analytic ground truth), §5 (Monte Carlo) and §6 (Greeks under noise, the flagship section).

- **Phase 1, analytic ground truth:** Black-Scholes-Merton with continuous dividend yield, closed-form Greeks, bracketed implied-vol solver, digital and discrete geometric-Asian closed forms, [`docs/conventions.md`](docs/conventions.md). The report shows why an in-the-money quote cannot pin its implied vol and the finite-difference V-curve.
- **Phase 2, stochastic infrastructure:** Philox 4×32-10, AS241 inverse normal, Welford/Chan accumulators, fixed-block reduction that is bit-identical for any thread count and across GCC/Clang, and the experiment harness (`experiments/` → `data/results/*.csv` + `.meta.json` → `docs/figures/*.svg`).
- **Phase 4, Monte Carlo:** a generic engine over the `PathModel` concept, antithetic and control variates, randomized QMC (Sobol, Owen scrambling, Brownian bridge). Variance reduction reaches 1,500× and QMC a 344× smaller error on an ATM call; the report shows where each technique fails.
- **Phase 5, Greeks under noise:** finite differences (independent seeds, common random numbers), pathwise (checked by forward automatic differentiation), likelihood ratio and mixed estimators, on a call and a digital. Measured convergence rates match theory (e.g. −0.40 against −2/5 for the CRN digital delta); the pathwise digital delta converges, with zero standard error, to 0 instead of 0.0188; the report ends with a payoff × regime × Greek recommendation matrix.
- **Next:** Phase 6, risk measures on non-linear positions (delta-normal VaR blind to a short straddle, delta-gamma, full revaluation, Expected Shortfall, backtesting).

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
cmake --build build-release
for e in fd_vcurve iv_roundtrip mc_convergence mc_coverage mc_efficiency qmc_convergence \
         greeks_vs_h greeks_vs_n greeks_matrix; do build-release/experiments/$e; done
python3 tools/make_figures.py
```

## License

MIT — see [LICENSE](LICENSE).
