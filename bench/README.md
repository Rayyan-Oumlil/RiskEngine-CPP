# Benchmarks

Micro-benchmarks behind report §9 (Google Benchmark, fetched by CMake). Off by default:

```
cmake -S . -B build-bench -DCMAKE_BUILD_TYPE=Release -DRISKENGINE_BUILD_BENCH=ON
cmake --build build-bench --target riskengine_bench
build-bench/bench/riskengine_bench --benchmark_repetitions=5 --benchmark_report_aggregates_only=true \
    --benchmark_out=data/results/bench.json --benchmark_out_format=json
```

The JSON records the machine (CPU, caches, load), and, through the experiment harness, the git commit,
compiler and flags. Timings belong to the machine that produced them: the committed
`data/results/bench.json` comes from a 4-vCPU KVM guest (see its context), not from bare metal.

What is measured:

- closed forms: Black-Scholes price and Greeks, implied vol (Brent on the log price);
- one Philox normal (uniform plus AS241 inversion);
- trees: CRR European and American at n = 100, 1,000, 10,000; Leisen-Reimer; BBS-Richardson;
- Monte Carlo on 2^20 paths of the canonical call: pseudo-random on 1 to 4 threads, antithetic,
  randomized QMC;
- a false-sharing study: per-thread counters adjacent in one cache line or padded to their own,
  atomic or plain, against the engine's local accumulators.
