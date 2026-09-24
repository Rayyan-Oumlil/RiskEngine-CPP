#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

namespace riskengine {

// Brownian-bridge construction of a path on the equally spaced grid t_k = k dt, k = 1..n.
//
// Maps n independent standard normals z, in order of importance, to the n standardized increments
// (W(t_k) - W(t_{k-1})) / sqrt(dt), which are again independent standard normals: z[0] sets W(T),
// z[1] the midpoint, then the quarter points, and so on (breadth first). The map is orthogonal,
// so it leaves pseudo-random simulation unchanged in distribution; with quasi-random points it
// puts the most important directions of the path on the best-distributed coordinates, which is
// what lowers the effective dimension of path-dependent payoffs (Glasserman 2003, 5.4).
class BrownianBridge {
public:
    explicit BrownianBridge(std::uint32_t steps) : steps_(steps) {
        assert(steps >= 1);
        // Work in units of dt, so W(k) has variance k. Point i of the construction is W(target[i])
        // given its already-built neighbours W(left[i]) and W(right[i]) (index 0 means W(0) = 0).
        target_.reserve(steps);
        target_.push_back(steps);
        left_.push_back(0);
        right_.push_back(0);
        left_weight_.push_back(0.0);
        right_weight_.push_back(0.0);
        sigma_.push_back(std::sqrt(static_cast<double>(steps)));
        std::vector<std::pair<std::uint32_t, std::uint32_t>> queue{{0, steps}};
        for (std::size_t q = 0; q < queue.size(); ++q) {
            const auto [l, r] = queue[q];
            if (r - l < 2) continue;
            const std::uint32_t m = l + (r - l) / 2;
            const double span = r - l;
            target_.push_back(m);
            left_.push_back(l);
            right_.push_back(r);
            left_weight_.push_back((r - m) / span);
            right_weight_.push_back((m - l) / span);
            sigma_.push_back(std::sqrt((m - l) * static_cast<double>(r - m) / span));
            queue.emplace_back(l, m);
            queue.emplace_back(m, r);
        }
        assert(target_.size() == steps);
    }

    std::uint32_t steps() const { return steps_; }

    // Writes the n standardized increments of the path built from z into `increments`.
    void transform(std::span<const double> z, std::span<double> increments) const {
        assert(z.size() == steps_ && increments.size() == steps_);
        // Build W(1..n) in place: W(k) is stored at increments[k - 1].
        const auto w = [&](std::uint32_t k) { return k == 0 ? 0.0 : increments[k - 1]; };
        increments[steps_ - 1] = sigma_[0] * z[0];
        for (std::uint32_t i = 1; i < steps_; ++i)
            increments[target_[i] - 1] =
                left_weight_[i] * w(left_[i]) + right_weight_[i] * w(right_[i]) + sigma_[i] * z[i];
        // Then difference, back to front: increment k = W(k) - W(k - 1), already in units of sqrt(dt).
        for (std::uint32_t k = steps_ - 1; k >= 1; --k) increments[k] -= increments[k - 1];
    }

private:
    std::uint32_t steps_;
    std::vector<std::uint32_t> target_, left_, right_;
    std::vector<double> left_weight_, right_weight_, sigma_;
};

} // namespace riskengine
