#pragma once

// Vector-width traits for the batched random-number kernels (core/rng/normal_icdf.hpp,
// core/rng/philox.hpp, core/rng/random_stream.hpp). Each kernel is written once, as a template over
// one of these, and runs 4 lanes wide with AVX2 or 8 lanes wide with AVX-512F.
//
// Every operation exposed here is either integer arithmetic or an IEEE 754 operation that is
// correctly rounded per lane (+, -, *, /, sqrt, sign changes, compares, selects). There is no fused
// multiply-add and no transcendental: that is what lets a kernel return exactly the bits of the
// scalar code, whatever the width. Only intrinsics from AVX2 and AVX-512F are used.

#include <cstdint>

#if defined(__AVX2__)
#include <immintrin.h>
#endif

namespace riskengine::simd {

#if defined(__AVX2__)
struct Avx2 {
    static constexpr int kLanes = 4;
    using Vec = __m256d;  // four doubles
    using Mask = __m256d; // lane mask, as a vector of all-ones / all-zeros lanes
    using Int = __m256i;  // four 64-bit integer lanes

    static Vec load(const double* p) { return _mm256_loadu_pd(p); }
    static void store(double* p, Vec v) { _mm256_storeu_pd(p, v); }
    static void store(double* p, Mask m, Vec v) { _mm256_maskstore_pd(p, _mm256_castpd_si256(m), v); }
    static Vec set1(double x) { return _mm256_set1_pd(x); }
    static Vec add(Vec a, Vec b) { return _mm256_add_pd(a, b); }
    static Vec sub(Vec a, Vec b) { return _mm256_sub_pd(a, b); }
    static Vec mul(Vec a, Vec b) { return _mm256_mul_pd(a, b); }
    static Vec div(Vec a, Vec b) { return _mm256_div_pd(a, b); }
    static Vec sqrt(Vec a) { return _mm256_sqrt_pd(a); }
    static Vec abs(Vec a) { return _mm256_andnot_pd(_mm256_set1_pd(-0.0), a); }
    static Vec neg(Vec a) { return _mm256_xor_pd(a, _mm256_set1_pd(-0.0)); } // exact sign flip, like -x
    static Mask le(Vec a, Vec b) { return _mm256_cmp_pd(a, b, _CMP_LE_OQ); } // false for NaN
    static Mask lt(Vec a, Vec b) { return _mm256_cmp_pd(a, b, _CMP_LT_OQ); }
    static unsigned bits(Mask m) { return static_cast<unsigned>(_mm256_movemask_pd(m)); }
    static Vec select(Mask m, Vec if_true, Vec if_false) { return _mm256_blendv_pd(if_false, if_true, m); }

    static Int set1(std::uint64_t x) { return _mm256_set1_epi64x(static_cast<long long>(x)); }
    static Int mul_lo32(Int a, Int b) { return _mm256_mul_epu32(a, b); } // low 32 x low 32 -> 64 bits
    static Int shr(Int a, int n) { return _mm256_srli_epi64(a, n); }
    static Int shl(Int a, int n) { return _mm256_slli_epi64(a, n); }
    static Int bit_and(Int a, Int b) { return _mm256_and_si256(a, b); }
    static Int bit_or(Int a, Int b) { return _mm256_or_si256(a, b); }
    static Int bit_xor(Int a, Int b) { return _mm256_xor_si256(a, b); }
    static Int from_lanes(const std::uint64_t (&v)[kLanes]) {
        return _mm256_loadu_si256(reinterpret_cast<const __m256i*>(v));
    }
    static Vec as_double(Int a) { return _mm256_castsi256_pd(a); }
    static Int as_int(Vec a) { return _mm256_castpd_si256(a); }

    // out[0 .. 2 * kLanes) = a0 b0 a1 b1 a2 b2 a3 b3.
    static void store_interleaved(double* out, Vec a, Vec b) {
        const Vec lo = _mm256_unpacklo_pd(a, b); // a0 b0 a2 b2
        const Vec hi = _mm256_unpackhi_pd(a, b); // a1 b1 a3 b3
        _mm256_storeu_pd(out, _mm256_permute2f128_pd(lo, hi, 0x20));
        _mm256_storeu_pd(out + 4, _mm256_permute2f128_pd(lo, hi, 0x31));
    }
};
#endif

#if defined(__AVX512F__)
struct Avx512 {
    static constexpr int kLanes = 8;
    using Vec = __m512d;
    using Mask = __mmask8; // one bit per lane
    using Int = __m512i;

    static Vec load(const double* p) { return _mm512_loadu_pd(p); }
    static void store(double* p, Vec v) { _mm512_storeu_pd(p, v); }
    static void store(double* p, Mask m, Vec v) { _mm512_mask_storeu_pd(p, m, v); } // untouched where m = 0
    static Vec set1(double x) { return _mm512_set1_pd(x); }
    static Vec add(Vec a, Vec b) { return _mm512_add_pd(a, b); }
    static Vec sub(Vec a, Vec b) { return _mm512_sub_pd(a, b); }
    static Vec mul(Vec a, Vec b) { return _mm512_mul_pd(a, b); }
    static Vec div(Vec a, Vec b) { return _mm512_div_pd(a, b); }
    static Vec sqrt(Vec a) { return _mm512_sqrt_pd(a); }
    static Vec abs(Vec a) { return _mm512_abs_pd(a); }
    // Sign flip through the integer unit: the double-typed xor needs AVX-512DQ, the integer one is F.
    static Vec neg(Vec a) {
        return _mm512_castsi512_pd(_mm512_xor_si512(_mm512_castpd_si512(a), _mm512_set1_epi64(INT64_MIN)));
    }
    static Mask le(Vec a, Vec b) { return _mm512_cmp_pd_mask(a, b, _CMP_LE_OQ); }
    static Mask lt(Vec a, Vec b) { return _mm512_cmp_pd_mask(a, b, _CMP_LT_OQ); }
    static unsigned bits(Mask m) { return static_cast<unsigned>(m); }
    static Vec select(Mask m, Vec if_true, Vec if_false) { return _mm512_mask_blend_pd(m, if_false, if_true); }

    static Int set1(std::uint64_t x) { return _mm512_set1_epi64(static_cast<long long>(x)); }
    static Int mul_lo32(Int a, Int b) { return _mm512_mul_epu32(a, b); }
    static Int shr(Int a, int n) { return _mm512_srli_epi64(a, static_cast<unsigned>(n)); }
    static Int shl(Int a, int n) { return _mm512_slli_epi64(a, static_cast<unsigned>(n)); }
    static Int bit_and(Int a, Int b) { return _mm512_and_si512(a, b); }
    static Int bit_or(Int a, Int b) { return _mm512_or_si512(a, b); }
    static Int bit_xor(Int a, Int b) { return _mm512_xor_si512(a, b); }
    static Int from_lanes(const std::uint64_t (&v)[kLanes]) { return _mm512_loadu_si512(v); }
    static Vec as_double(Int a) { return _mm512_castsi512_pd(a); }
    static Int as_int(Vec a) { return _mm512_castpd_si512(a); }

    // out[0 .. 2 * kLanes) = a0 b0 a1 b1 ... a7 b7.
    static void store_interleaved(double* out, Vec a, Vec b) {
        const __m512i lo = _mm512_set_epi64(11, 3, 10, 2, 9, 1, 8, 0);    // a0 b0 a1 b1 a2 b2 a3 b3
        const __m512i hi = _mm512_set_epi64(15, 7, 14, 6, 13, 5, 12, 4);  // a4 b4 ... a7 b7
        _mm512_storeu_pd(out, _mm512_permutex2var_pd(a, lo, b));
        _mm512_storeu_pd(out + 8, _mm512_permutex2var_pd(a, hi, b));
    }
};
#endif

// The widest width this build supports; kernels run at this width. Absent when neither is enabled.
#if defined(__AVX512F__)
using Widest = Avx512;
#define RISKENGINE_SIMD 1
#elif defined(__AVX2__)
using Widest = Avx2;
#define RISKENGINE_SIMD 1
#endif

} // namespace riskengine::simd
