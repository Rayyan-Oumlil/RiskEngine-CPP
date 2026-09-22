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

Early scaffold (Phase 0). See [`docs/riskengine_research.md`](docs/riskengine_research.md) §12 for the phased build plan, validation gates, and critical path.

## Building

```
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

Requires a C++20 compiler (GCC, Clang, or MSVC) and CMake 3.25+.

## License

MIT — see [LICENSE](LICENSE).
