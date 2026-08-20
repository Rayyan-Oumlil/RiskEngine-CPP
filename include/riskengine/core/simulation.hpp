#pragma once

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <exception>
#include <mutex>
#include <thread>
#include <type_traits>
#include <vector>

#include "riskengine/core/estimate.hpp"
#include "riskengine/core/rng/random_stream.hpp"
#include "riskengine/core/stats/welford.hpp"

namespace riskengine {

// How a simulation is split. The work is cut into a fixed number of blocks, independent of the
// thread count; each block owns its random stream and its accumulator, and the partial results
// are merged in block order. Floating-point addition is not associative, so this fixed order is
// what makes the result bit-identical with 1 or 64 threads. Changing `blocks` changes the result
// (within statistical error); changing `threads` never does.
struct BlockPlan {
    std::uint64_t samples;
    std::uint32_t blocks = 64;
    unsigned threads = 1;

    std::uint64_t block_samples(std::uint32_t b) const {
        const std::uint64_t base = samples / blocks, extra = samples % blocks;
        return base + (b < extra ? 1 : 0);
    }
};

// Runs block_fn(block, n) -> Acc for every block and merges the results in block order. Acc is any
// default-constructible accumulator with merge(const Acc&) (Welford, Covariance). block_fn is called
// concurrently from several threads and must not share mutable state between blocks. The first
// exception thrown by a block is rethrown once all threads have stopped.
template <class BlockFn>
auto reduce_blocks(const BlockPlan& plan, BlockFn&& block_fn) {
    using Acc = std::remove_cvref_t<std::invoke_result_t<BlockFn&, std::uint32_t, std::uint64_t>>;
    assert(plan.blocks > 0);
    std::vector<Acc> partial(plan.blocks);
    std::atomic<std::uint32_t> next{0};
    std::exception_ptr error;
    std::mutex error_mutex;

    auto worker = [&] {
        for (std::uint32_t b = next++; b < plan.blocks; b = next++) {
            try {
                partial[b] = block_fn(b, plan.block_samples(b));
            } catch (...) {
                const std::lock_guard lock(error_mutex);
                if (!error) error = std::current_exception();
                next = plan.blocks; // stop handing out work
            }
        }
    };

    const unsigned threads = std::clamp(plan.threads, 1u, std::max(plan.blocks, 1u));
    {
        std::vector<std::jthread> pool;
        pool.reserve(threads - 1);
        for (unsigned i = 1; i < threads; ++i) pool.emplace_back(worker);
        worker();
    } // jthreads join here
    if (error) std::rethrow_exception(error);

    Acc total;
    for (const Acc& p : partial) total.merge(p);
    return total;
}

// Monte Carlo mean of sample(rng) -> double, with its standard error. Block b draws from
// RandomStream(key, b), so the estimate is a pure function of (key, plan.samples, plan.blocks).
template <class SampleFn>
Estimate simulate(SeedKey key, const BlockPlan& plan, SampleFn&& sample) {
    const Welford w = reduce_blocks(plan, [&](std::uint32_t b, std::uint64_t n) {
        RandomStream rng(key, b);
        Welford acc;
        for (std::uint64_t i = 0; i < n; ++i) acc.add(sample(rng));
        return acc;
    });
    return to_estimate(w);
}

} // namespace riskengine
