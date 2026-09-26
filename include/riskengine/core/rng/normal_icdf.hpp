#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>

#if defined(__AVX2__)
#include <immintrin.h>
#endif

namespace riskengine {

namespace detail {

// Horner evaluation of c[0] + c[1] x + ... + c[7] x^7, in a fixed operation order.
inline double poly7(const double (&c)[8], double x) {
    double r = c[7];
    for (int i = 6; i >= 0; --i) r = r * x + c[i];
    return r;
}

// AS241 central-region coefficients (|p - 0.5| <= 0.425), shared by the scalar and batched paths
// so the two can never drift apart.
inline constexpr double kIcdfCentralNum[8] = {
    3.3871328727963666080e0,  1.3314166789178437745e+2, 1.9715909503065514427e+3, 1.3731693765509461125e+4,
    4.5921953931549871457e+4, 6.7265770927008700853e+4, 3.3430575583588128105e+4, 2.5090809287301226727e+3};
inline constexpr double kIcdfCentralDen[8] = {
    1.0,                      4.2313330701600911252e+1, 6.8718700749205790830e+2, 5.3941960214247511077e+3,
    2.1213794301586595867e+4, 3.9307895800092710610e+4, 2.8729085735721942674e+4, 5.2264952788528545610e+3};

// AS241 tail coefficients: sqrt(-log(min(p, 1 - p))) <= 5 (Num/Den) and beyond (FarNum/FarDen).
inline constexpr double kIcdfTailNum[8] = {
    1.42343711074968357734e0, 4.63033784615654529590e0, 5.76949722146069140550e0, 3.64784832476320460504e0,
    1.27045825245236838258e0, 2.41780725177450611770e-1, 2.27238449892691845833e-2, 7.74545014278341407640e-4};
inline constexpr double kIcdfTailDen[8] = {
    1.0, 2.05319162663775882187e0, 1.67638483018380384940e0, 6.89767334985100004550e-1,
    1.48103976427480074590e-1, 1.51986665636164571966e-2, 5.47593808499534494600e-4, 1.05075007164441684324e-9};
inline constexpr double kIcdfFarTailNum[8] = {
    6.65790464350110377720e0, 5.46378491116411436990e0, 1.78482653991729133580e0, 2.96560571828504891230e-1,
    2.65321895265761230930e-2, 1.24266094738807843860e-3, 2.71155556874348757815e-5, 2.01033439929228813265e-7};
inline constexpr double kIcdfFarTailDen[8] = {
    1.0, 5.99832206555887937690e-1, 1.36929880922735805310e-1, 1.48753612908506148525e-2,
    7.86869131145613259100e-4, 1.84631831751005468180e-5, 1.42151175831644588870e-7, 2.04426310338993978564e-15};

} // namespace detail

// Inverse of the standard normal CDF: Wichura (1988), Algorithm AS 241 (PPND16), relative
// accuracy about 1e-16 over (0, 1).
//
// Used instead of std::normal_distribution, whose algorithm is implementation-defined (the same
// seed gives different numbers under libstdc++, libc++ and MSVC), and instead of Box-Muller,
// which destroys the structure of quasi-random points. One uniform in, one normal out, monotone:
// the mapping preserves stratification and antithetic symmetry.
inline double norm_icdf(double p) {

    if (!(p > 0.0 && p < 1.0)) {
        if (p == 0.0) return -std::numeric_limits<double>::infinity();
        if (p == 1.0) return std::numeric_limits<double>::infinity();
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double q = p - 0.5;
    if (std::abs(q) <= 0.425) {
        const double r = 0.180625 - q * q;
        return q * detail::poly7(detail::kIcdfCentralNum, r) / detail::poly7(detail::kIcdfCentralDen, r);
    }

    // Tails: work with the smaller of p and 1 - p, so the left tail keeps full relative precision.
    double r = std::sqrt(-std::log(q < 0.0 ? p : 1.0 - p));
    double x;
    if (r <= 5.0) {
        r -= 1.6;
        x = detail::poly7(detail::kIcdfTailNum, r) / detail::poly7(detail::kIcdfTailDen, r);
    } else {
        r -= 5.0;
        x = detail::poly7(detail::kIcdfFarTailNum, r) / detail::poly7(detail::kIcdfFarTailDen, r);
    }
    return q < 0.0 ? -x : x;
}

#if defined(__AVX2__)
namespace detail {

// poly7 on four lanes, in exactly the scalar operation order (multiply, then add, as two separately
// rounded IEEE operations). No FMA intrinsic: a fused multiply-add rounds once instead of twice and
// would give different bits from the scalar path, which the project builds with -ffp-contract=off.
inline __m256d poly7_x4(const double (&c)[8], __m256d x) {
    __m256d r = _mm256_set1_pd(c[7]);
    for (int i = 6; i >= 0; --i) r = _mm256_add_pd(_mm256_mul_pd(r, x), _mm256_set1_pd(c[i]));
    return r;
}

// The tail formula of norm_icdf on four lanes, given q = p - 0.5 and log(min(p, 1 - p)) computed by
// the caller with the scalar std::log (the one operation here that must stay scalar).
inline __m256d icdf_tail_x4(__m256d q, __m256d log_min) {
    const __m256d sign_bit = _mm256_set1_pd(-0.0);
    const __m256d r = _mm256_sqrt_pd(_mm256_xor_pd(log_min, sign_bit)); // sqrt(-log): both exact
    const __m256d near = _mm256_cmp_pd(r, _mm256_set1_pd(5.0), _CMP_LE_OQ);
    const __m256d r1 = _mm256_sub_pd(r, _mm256_set1_pd(1.6));
    const __m256d r2 = _mm256_sub_pd(r, _mm256_set1_pd(5.0));
    const __m256d x1 = _mm256_div_pd(poly7_x4(kIcdfTailNum, r1), poly7_x4(kIcdfTailDen, r1));
    const __m256d x2 = _mm256_div_pd(poly7_x4(kIcdfFarTailNum, r2), poly7_x4(kIcdfFarTailDen, r2));
    const __m256d x = _mm256_blendv_pd(x2, x1, near);
    const __m256d negative = _mm256_cmp_pd(q, _mm256_setzero_pd(), _CMP_LT_OQ);
    return _mm256_blendv_pd(x, _mm256_xor_pd(x, sign_bit), negative);
}

} // namespace detail
#endif

// norm_icdf over a span, bit-identical to calling the scalar norm_icdf on each element: the same
// values, not merely close ones, so results stay reproducible whether or not a build has AVX2.
//
// With AVX2, every operation that IEEE 754 rounds exactly (+, -, *, /, sqrt, sign flips) runs four
// lanes at a time, in the scalar code's order, so each lane computes the scalar bits. Only std::log
// stays scalar: vectorized logarithms round differently from libm. The work is split in two passes
// per block, so a chunk with one tail lane does not drag its three central lanes back to scalar:
//   1. central lanes (|p - 0.5| <= 0.425, about 85 %) are computed and written with a masked store;
//      tail lane indices are collected;
//   2. tail lanes are gathered four at a time: scalar log, then the rest of the formula vectorized.
// 0, 1, NaN and values outside (0, 1) take the scalar function. `p` and `out` may be the same buffer:
// pass 1 never writes a tail lane, so pass 2 still reads its original input.
inline void norm_icdf(std::span<const double> p, std::span<double> out) {
    assert(p.size() == out.size());
    std::size_t i = 0;
#if defined(__AVX2__)
    constexpr std::size_t kBlock = 256;
    const __m256d half = _mm256_set1_pd(0.5);
    const __m256d bound = _mm256_set1_pd(0.425);
    const __m256d shift = _mm256_set1_pd(0.180625);
    const __m256d sign_bit = _mm256_set1_pd(-0.0);
    while (p.size() - i >= 4) {
        const std::size_t end = i + std::min(kBlock, (p.size() - i) & ~std::size_t{3});
        std::uint32_t tail[kBlock];
        std::size_t tails = 0;

        for (std::size_t j = i; j < end; j += 4) { // pass 1: central lanes
            const __m256d q = _mm256_sub_pd(_mm256_loadu_pd(p.data() + j), half);
            const __m256d central = _mm256_cmp_pd(_mm256_andnot_pd(sign_bit, q), bound, _CMP_LE_OQ); // NaN: false
            const __m256d r = _mm256_sub_pd(shift, _mm256_mul_pd(q, q));
            const __m256d x = _mm256_div_pd(_mm256_mul_pd(q, detail::poly7_x4(detail::kIcdfCentralNum, r)),
                                            detail::poly7_x4(detail::kIcdfCentralDen, r));
            const int mask = _mm256_movemask_pd(central);
            if (mask == 0xF) {
                _mm256_storeu_pd(out.data() + j, x);
                continue;
            }
            _mm256_maskstore_pd(out.data() + j, _mm256_castpd_si256(central), x);
            // Branchless compaction: tail lanes fall at random positions, so an `if` here would
            // mispredict on ~15 % of lanes. Always write the index; advance the count by 0 or 1.
            for (int k = 0; k < 4; ++k) {
                tail[tails] = static_cast<std::uint32_t>(j - i + k);
                tails += ~mask >> k & 1;
            }
        }

        for (std::size_t t = 0; t < tails; t += 4) { // pass 2: tail lanes, four at a time
            const std::size_t n = std::min<std::size_t>(4, tails - t);
            alignas(32) double pin[4] = {0.25, 0.25, 0.25, 0.25}, log_min[4] = {-1.0, -1.0, -1.0, -1.0};
            bool scalar[4] = {false, false, false, false};
            for (std::size_t k = 0; k < n; ++k) {
                const double v = p[i + tail[t + k]];
                pin[k] = v;
                if (!(v > 0.0 && v < 1.0)) { scalar[k] = true; continue; } // 0, 1, NaN, out of range
                // Same value as the scalar (q < 0 ? p : 1 - p): for v < 0.5 the smaller is v, for v > 0.5
                // it is 1 - v, and v = 0.5 never reaches the tail. std::min compiles to a branchless
                // minsd, where the ternary is a coin-flip branch on random tails.
                log_min[k] = std::log(std::min(v, 1.0 - v));
            }
            alignas(32) double x[4];
            _mm256_store_pd(x, detail::icdf_tail_x4(_mm256_sub_pd(_mm256_load_pd(pin), half), _mm256_load_pd(log_min)));
            for (std::size_t k = 0; k < n; ++k) out[i + tail[t + k]] = scalar[k] ? norm_icdf(pin[k]) : x[k];
        }
        i = end;
    }
#endif
    for (; i < p.size(); ++i) out[i] = norm_icdf(p[i]);
}

} // namespace riskengine
