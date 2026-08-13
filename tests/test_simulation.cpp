#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "riskengine/core/simulation.hpp"

using namespace riskengine;

namespace {

double exp_normal(RandomStream& rng) { return std::exp(rng.normal()); }

} // namespace

TEST_CASE("Welford matches a two-pass computation and merges exactly as one stream", "[welford]") {
    std::vector<double> xs;
    RandomStream rng(SeedKey{5}, 0);
    for (int i = 0; i < 10'001; ++i) xs.push_back(1e6 + rng.normal()); // large offset: naive sum of squares cancels

    Welford all;
    for (double x : xs) all.add(x);
    // Accurate two-pass reference: x - 1e6 is exact here (Sterbenz), so the offset costs nothing.
    const double n = static_cast<double>(xs.size());
    double dev_sum = 0.0;
    for (double x : xs) dev_sum += x - 1e6;
    const double m = 1e6 + dev_sum / n;
    double ss = 0.0;
    for (double x : xs) ss += (x - m) * (x - m);
    CHECK(all.count() == xs.size());
    // Running-mean rounding grows like sqrt(n) ulps of the mean (~100 ulps = 2e-14 here).
    CHECK(std::abs(all.mean() - m) <= 2e-14 * m);
    CHECK(std::abs(all.variance() - ss / (n - 1)) <= 1e-9 * ss / (n - 1));

    Welford left, right;
    for (std::size_t i = 0; i < xs.size(); ++i) (i < 3'000 ? left : right).add(xs[i]);
    left.merge(right);
    CHECK(left.count() == all.count());
    CHECK(std::abs(left.mean() - all.mean()) <= 2e-14 * all.mean());
    CHECK(std::abs(left.variance() - all.variance()) <= 1e-9 * all.variance());

    Welford empty;
    empty.merge(all);
    CHECK(empty.mean() == all.mean());
    all.merge(Welford{});
    CHECK(all.count() == xs.size());
    CHECK(Welford{}.variance() == 0.0);
    CHECK(Welford{}.std_error() == 0.0);
}

TEST_CASE("Block plan splits samples exactly", "[simulation]") {
    const BlockPlan plan{1003, 10};
    std::uint64_t total = 0;
    for (std::uint32_t b = 0; b < plan.blocks; ++b) total += plan.block_samples(b);
    CHECK(total == 1003);
    CHECK(plan.block_samples(0) == 101);
    CHECK(plan.block_samples(9) == 100);
}

TEST_CASE("Simulation is bit-identical for 1, 2, 3, 8 and 64 threads", "[simulation][determinism]") {
    // Phase 2 gate. More threads than blocks and uneven splits included.
    const SeedKey key{2026};
    const Estimate reference = simulate(key, BlockPlan{100'000, 16, 1}, exp_normal);
    for (unsigned threads : {2u, 3u, 8u, 64u}) {
        INFO("threads=" << threads);
        const Estimate e = simulate(key, BlockPlan{100'000, 16, threads}, exp_normal);
        CHECK(e.value == reference.value);
        CHECK(e.std_error == reference.std_error);
        CHECK(e.samples == reference.samples);
    }
}

TEST_CASE("Simulation golden value is shared by GCC and Clang", "[simulation][golden]") {
    // Phase 2 gate: same bits on two compilers (exp and log come from the platform libm, so the
    // exact check is limited to platforms that share it).
    const Estimate e = simulate(SeedKey{2026}, BlockPlan{100'000, 16, 4}, exp_normal);
    CHECK(e.samples == 100'000);
    CHECK(e.discretization == 0.0);
    // E[e^Z] = e^{1/2}: statistical sanity, at 4 standard errors.
    CHECK(std::abs(e.value - std::exp(0.5)) <= 4.0 * e.std_error);
#if !defined(_MSC_VER)
    CHECK(e.value == 0x1.a18aa40e2060dp+0);
    CHECK(e.std_error == 0x1.b654fa0d7ec3ep-8);
#endif
}

TEST_CASE("Same key reproduces; a different stream is independent", "[simulation]") {
    const BlockPlan plan{10'000, 8, 4};
    const Estimate a = simulate(SeedKey{1}, plan, exp_normal);
    const Estimate b = simulate(SeedKey{1}, plan, exp_normal);
    const Estimate c = simulate(SeedKey{1, 1}, plan, exp_normal);
    CHECK(a.value == b.value);
    CHECK(a.value != c.value);
}

TEST_CASE("Standard errors are calibrated across independent replications", "[simulation][statistical]") {
    // 200 independent replications of the mean of a uniform (true value 1/2). The z-scores must
    // look standard normal: roughly 95 % inside +-1.96 and none beyond +-4.
    int inside_95 = 0, beyond_4 = 0;
    const int replications = 200;
    for (int s = 0; s < replications; ++s) {
        const Estimate e = simulate(SeedKey{static_cast<std::uint64_t>(s)}, BlockPlan{20'000, 8, 4},
                                    [](RandomStream& rng) { return rng.uniform(); });
        const double z = (e.value - 0.5) / e.std_error;
        inside_95 += std::abs(z) <= 1.96 ? 1 : 0;
        beyond_4 += std::abs(z) > 4.0 ? 1 : 0;
    }
    CHECK(beyond_4 == 0);
    // Binomial(200, 0.95): mean 190, sd ~3.1; the seeds are fixed, so this is deterministic.
    CHECK(inside_95 >= 180);
}

TEST_CASE("An exception in a block is propagated after all threads stop", "[simulation]") {
    const BlockPlan plan{1'000, 16, 4};
    auto failing = [](std::uint32_t b, std::uint64_t) -> Welford {
        if (b == 5) throw std::runtime_error("block 5 failed");
        return Welford{};
    };
    CHECK_THROWS_AS(reduce_blocks(plan, failing), std::runtime_error);
}
