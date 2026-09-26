// Micro-benchmarks behind report 9: the cost of each pricing method at the sizes the report uses,
// the scaling of the Monte Carlo engine with threads, and a false-sharing study of the kind of
// accumulator layout the engine avoids.
//
// Run from a Release build (see bench/README.md). Every Monte Carlo benchmark prices the canonical
// ATM call (S = K = 100, r = 5 %, q = 0, sigma = 20 %, T = 1).

#include <benchmark/benchmark.h>

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <new>
#include <thread>
#include <vector>

#include "riskengine/core/rng/random_stream.hpp"
#include "riskengine/methods/analytic/black_scholes.hpp"
#include "riskengine/methods/analytic/implied_vol.hpp"
#include "riskengine/methods/montecarlo/engine.hpp"
#include "riskengine/methods/montecarlo/longstaff_schwartz.hpp"
#include "riskengine/methods/tree/binomial.hpp"
#include "riskengine/models/gbm.hpp"
#include "riskengine/payoffs/vanilla.hpp"

#ifdef RISKENGINE_HAVE_HARNESS
#include "harness/experiment.hpp"
#endif

using namespace riskengine;

namespace {

const MarketState kMarket{Spot{100}, Rate{0.05}, Rate{0.0}, Vol{0.2}};
const VanillaOption kCall{Strike{100}, Maturity{1}, OptionType::Call};
const VanillaOption kPut{Strike{100}, Maturity{1}, OptionType::Put};

// --- Closed forms ---------------------------------------------------------------------------------

void BM_BlackScholesPrice(benchmark::State& state) {
    MarketState m = kMarket;
    for (auto _ : state) {
        benchmark::DoNotOptimize(m);
        benchmark::DoNotOptimize(black_scholes_price(kCall, m));
    }
}
BENCHMARK(BM_BlackScholesPrice);

void BM_BlackScholesGreeks(benchmark::State& state) {
    MarketState m = kMarket;
    for (auto _ : state) {
        benchmark::DoNotOptimize(m);
        benchmark::DoNotOptimize(black_scholes_greeks(kCall, m));
    }
}
BENCHMARK(BM_BlackScholesGreeks);

// Implied vol of an at-the-money and a 25 % out-of-the-money call (Brent on the log price). The
// quotes are priced at vols around 30 % (not 20 %, the solver's first bracket point, where a solve
// ends at once) and read through a volatile array, so the solve cannot be hoisted out of the loop.
// Not benchmark::DoNotOptimize on the double: with GCC 13, its "+m,r" constraint on a double lvalue
// made the loop read 100.0 instead of the quote, and every solve returned an error in 9 ns.
void BM_ImpliedVol(benchmark::State& state) {
    const VanillaOption o{Strike{static_cast<double>(state.range(0))}, Maturity{1}, OptionType::Call};
    volatile double quotes[8];
    for (int i = 0; i < 8; ++i) {
        MarketState quoted = kMarket;
        quoted.vol = Vol{0.28 + 0.005 * i};
        quotes[i] = black_scholes_price(o, quoted);
    }
    std::size_t i = 0, failures = 0;
    for (auto _ : state) {
        ImpliedVolResult r = implied_vol(o, quotes[i++ & 7], kMarket.spot, kMarket.rate, kMarket.div);
        failures += r.ok() ? 0 : 1;
        benchmark::DoNotOptimize(r);
    }
    if (failures > 0) state.SkipWithError("implied_vol failed on a benchmark quote");
}
BENCHMARK(BM_ImpliedVol)->Arg(100)->Arg(125);

// --- Random numbers -------------------------------------------------------------------------------

void BM_PhiloxNormals(benchmark::State& state) {
    RandomStream rng(SeedKey{1}, 0);
    for (auto _ : state) benchmark::DoNotOptimize(rng.normal());
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_PhiloxNormals);

// --- Trees ----------------------------------------------------------------------------------------

void BM_Tree(benchmark::State& state, TreeMethod method, Exercise exercise) {
    const auto n = static_cast<unsigned>(state.range(0));
    const VanillaOption& o = exercise == Exercise::American ? kPut : kCall;
    for (auto _ : state) benchmark::DoNotOptimize(binomial_price(method, o, exercise, kMarket, n));
    state.SetComplexityN(state.range(0));
}
BENCHMARK_CAPTURE(BM_Tree, crr_european, TreeMethod::Crr, Exercise::European)->Arg(100)->Arg(1000)->Arg(10000)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_CAPTURE(BM_Tree, crr_american, TreeMethod::Crr, Exercise::American)->Arg(100)->Arg(1000)->Arg(10000)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_CAPTURE(BM_Tree, leisen_reimer_european, TreeMethod::LeisenReimer, Exercise::European)->Arg(101)
    ->Arg(1001)->Unit(benchmark::kMicrosecond);
BENCHMARK_CAPTURE(BM_Tree, bbs_richardson_american, TreeMethod::BbsRichardson, Exercise::American)->Arg(1000)
    ->Unit(benchmark::kMicrosecond);

// --- Monte Carlo ----------------------------------------------------------------------------------

// 2^20 paths of the canonical call, on range(0) threads (the result is the same for every count).
void BM_MonteCarlo(benchmark::State& state, Sampling sampling, bool antithetic) {
    MonteCarloConfig config{.paths = 1u << 20, .threads = static_cast<unsigned>(state.range(0)),
                            .antithetic = antithetic, .sampling = sampling, .replications = 16};
    if (sampling == Sampling::RandomizedQmc) config.paths /= config.replications; // same total points
    const MonteCarlo<GBM, VanillaPayoff> pricer(VanillaPayoff{100.0, OptionType::Call}, Maturity{1}, config);
    for (auto _ : state) benchmark::DoNotOptimize(pricer.price(kMarket, SeedKey{2026}));
    state.SetItemsProcessed(state.iterations() * (1 << 20));
}
BENCHMARK_CAPTURE(BM_MonteCarlo, pseudo_random, Sampling::PseudoRandom, false)->DenseRange(1, 4)
    ->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK_CAPTURE(BM_MonteCarlo, antithetic, Sampling::PseudoRandom, true)->Arg(1)->Unit(benchmark::kMillisecond)
    ->UseRealTime();
BENCHMARK_CAPTURE(BM_MonteCarlo, rqmc, Sampling::RandomizedQmc, false)->Arg(1)->Unit(benchmark::kMillisecond)
    ->UseRealTime();

// --- False sharing --------------------------------------------------------------------------------

// Each of range(0) threads counts to 2^22 in its own slot; only where the slots live differs.
//   adjacent_atomic  atomic counters next to each other: one cache line shared by every thread, and
//                    every increment must own the line (a lock-prefixed read-modify-write);
//   padded_atomic    the same atomics, each on its own cache line;
//   adjacent_store   plain (volatile) doubles next to each other: the same sharing, but each thread
//                    reads its own last store back from its store buffer, which hides most of it;
//   local            a local sum written once at the end: what reduce_blocks does.
constexpr std::size_t kLine = 64; // std::hardware_destructive_interference_size is not in every library
struct alignas(kLine) PaddedCounter {
    std::atomic<std::uint64_t> count{0};
};
enum class Layout { AdjacentAtomic, PaddedAtomic, AdjacentStore, Local };

void BM_FalseSharing(benchmark::State& state, Layout layout) {
    const auto threads = static_cast<unsigned>(state.range(0));
    constexpr std::uint64_t kPerThread = 1u << 22;
    std::vector<std::atomic<std::uint64_t>> adjacent(threads);
    std::vector<PaddedCounter> padded(threads);
    std::vector<double> stores(threads);
    for (auto _ : state) {
        std::vector<std::jthread> pool;
        for (unsigned t = 0; t < threads; ++t)
            pool.emplace_back([&, t] {
                switch (layout) {
                    case Layout::AdjacentAtomic:
                    case Layout::PaddedAtomic: {
                        auto& c = layout == Layout::AdjacentAtomic ? adjacent[t] : padded[t].count;
                        for (std::uint64_t i = 0; i < kPerThread; ++i) c.fetch_add(1, std::memory_order_relaxed);
                        break;
                    }
                    case Layout::AdjacentStore: {
                        volatile double& slot = stores[t];
                        for (std::uint64_t i = 0; i < kPerThread; ++i) slot = slot + 1.0;
                        break;
                    }
                    case Layout::Local: {
                        volatile double sum = 0.0; // the same load and store per step, on this thread's stack
                        for (std::uint64_t i = 0; i < kPerThread; ++i) sum = sum + 1.0;
                        stores[t] = sum;
                        break;
                    }
                }
            });
    }
    state.SetItemsProcessed(state.iterations() * threads * kPerThread);
}
BENCHMARK_CAPTURE(BM_FalseSharing, adjacent_atomic, Layout::AdjacentAtomic)->DenseRange(1, 4)
    ->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK_CAPTURE(BM_FalseSharing, padded_atomic, Layout::PaddedAtomic)->DenseRange(1, 4)
    ->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK_CAPTURE(BM_FalseSharing, adjacent_store, Layout::AdjacentStore)->DenseRange(1, 4)
    ->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK_CAPTURE(BM_FalseSharing, local, Layout::Local)->DenseRange(1, 4)->Unit(benchmark::kMillisecond)
    ->UseRealTime();

// --- Path storage: vector-of-vectors against one contiguous allocation ---------------------------

// Longstaff-Schwartz (methods/montecarlo/longstaff_schwartz.hpp) stores every simulated path to
// regress on later. A std::vector<std::vector<double>> makes `rows` separate heap allocations,
// scattered with no locality guarantee between them, even though every consumer walks whole rows
// in order. PathMatrix is the same shape as one contiguous buffer. Both benchmarks do identical
// work (write `cols` doubles per row, `rows` rows) so the difference is purely the allocation and
// locality pattern, not the arithmetic.
void BM_PathStorageVectorOfVectors(benchmark::State& state) {
    const auto rows = static_cast<std::size_t>(state.range(0));
    constexpr std::size_t kCols = 50;
    for (auto _ : state) {
        std::vector<std::vector<double>> paths(rows, std::vector<double>(kCols));
        for (std::size_t i = 0; i < rows; ++i)
            for (std::size_t k = 0; k < kCols; ++k) paths[i][k] = static_cast<double>(i + k);
        benchmark::DoNotOptimize(paths);
    }
    state.SetItemsProcessed(state.iterations() * rows * kCols);
}
BENCHMARK(BM_PathStorageVectorOfVectors)->Arg(1 << 12)->Arg(1 << 16)->Arg(1 << 18)->Unit(benchmark::kMillisecond);

void BM_PathStorageMatrix(benchmark::State& state) {
    const auto rows = static_cast<std::size_t>(state.range(0));
    constexpr std::size_t kCols = 50;
    for (auto _ : state) {
        detail::PathMatrix paths(rows, kCols);
        for (std::size_t i = 0; i < rows; ++i) {
            const std::span<double> row = paths.row(i);
            for (std::size_t k = 0; k < kCols; ++k) row[k] = static_cast<double>(i + k);
        }
        benchmark::DoNotOptimize(paths);
    }
    state.SetItemsProcessed(state.iterations() * rows * kCols);
}
BENCHMARK(BM_PathStorageMatrix)->Arg(1 << 12)->Arg(1 << 16)->Arg(1 << 18)->Unit(benchmark::kMillisecond);

// End to end: one Longstaff-Schwartz price of the canonical American put, 50 exercise dates.
void BM_LongstaffSchwartz(benchmark::State& state) {
    const auto paths = static_cast<std::uint64_t>(state.range(0));
    const LongstaffSchwartz<GBM> lsm(kPut, LongstaffSchwartzConfig{.paths = paths, .steps = 50});
    for (auto _ : state) benchmark::DoNotOptimize(lsm.price(kMarket, SeedKey{2026}));
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(paths));
}
BENCHMARK(BM_LongstaffSchwartz)->Arg(1 << 16)->Arg(1 << 18)->Unit(benchmark::kMillisecond);

// --- Batched normals -------------------------------------------------------------------------------

// 4,096 normals per iteration from one stream, drawn one by one (normal()) or as one batch
// (normals()). The two produce identical bits; with RISKENGINE_ENABLE_AVX2 the batch inverts the
// central 85 % of uniforms four at a time and keeps only std::log scalar (core/rng/normal_icdf.hpp).
void BM_NormalsOneByOne(benchmark::State& state) {
    std::vector<double> z(4096);
    RandomStream rng(SeedKey{11}, 0);
    for (auto _ : state) {
        for (double& v : z) v = rng.normal();
        benchmark::DoNotOptimize(z.data());
    }
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(z.size()));
}
BENCHMARK(BM_NormalsOneByOne);

void BM_NormalsBatched(benchmark::State& state) {
    std::vector<double> z(4096);
    RandomStream rng(SeedKey{11}, 0);
    for (auto _ : state) {
        rng.normals(z);
        benchmark::DoNotOptimize(z.data());
    }
    state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(z.size()));
}
BENCHMARK(BM_NormalsBatched);

} // namespace

int main(int argc, char** argv) {
    benchmark::Initialize(&argc, argv);
    if (benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
#ifdef RISKENGINE_HAVE_HARNESS
    benchmark::AddCustomContext("git_commit", build_info::kGitCommit);
    benchmark::AddCustomContext("git_dirty", build_info::kGitDirty ? "true" : "false");
    benchmark::AddCustomContext("compiler", RISKENGINE_COMPILER);
    benchmark::AddCustomContext("build_type", RISKENGINE_BUILD_TYPE);
    benchmark::AddCustomContext("cxx_flags", RISKENGINE_CXX_FLAGS);
    benchmark::AddCustomContext("cpu", harness::detail::cpu_model());
#endif
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
