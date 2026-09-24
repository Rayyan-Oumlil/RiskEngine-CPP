#pragma once

#include <cmath>
#include <limits>

namespace riskengine {

namespace detail {

// Horner evaluation of c[0] + c[1] x + ... + c[7] x^7, in a fixed operation order.
inline double poly7(const double (&c)[8], double x) {
    double r = c[7];
    for (int i = 6; i >= 0; --i) r = r * x + c[i];
    return r;
}

} // namespace detail

// Inverse of the standard normal CDF: Wichura (1988), Algorithm AS 241 (PPND16), relative
// accuracy about 1e-16 over (0, 1).
//
// Used instead of std::normal_distribution, whose algorithm is implementation-defined (the same
// seed gives different numbers under libstdc++, libc++ and MSVC), and instead of Box-Muller,
// which destroys the structure of quasi-random points. One uniform in, one normal out, monotone:
// the mapping preserves stratification and antithetic symmetry.
inline double norm_icdf(double p) {
    constexpr double a[8] = {3.3871328727963666080e0, 1.3314166789178437745e+2, 1.9715909503065514427e+3,
                             1.3731693765509461125e+4, 4.5921953931549871457e+4, 6.7265770927008700853e+4,
                             3.3430575583588128105e+4, 2.5090809287301226727e+3};
    constexpr double b[8] = {1.0, 4.2313330701600911252e+1, 6.8718700749205790830e+2, 5.3941960214247511077e+3,
                             2.1213794301586595867e+4, 3.9307895800092710610e+4, 2.8729085735721942674e+4,
                             5.2264952788528545610e+3};
    constexpr double c[8] = {1.42343711074968357734e0, 4.63033784615654529590e0, 5.76949722146069140550e0,
                             3.64784832476320460504e0, 1.27045825245236838258e0, 2.41780725177450611770e-1,
                             2.27238449892691845833e-2, 7.74545014278341407640e-4};
    constexpr double d[8] = {1.0, 2.05319162663775882187e0, 1.67638483018380384940e0, 6.89767334985100004550e-1,
                             1.48103976427480074590e-1, 1.51986665636164571966e-2, 5.47593808499534494600e-4,
                             1.05075007164441684324e-9};
    constexpr double e[8] = {6.65790464350110377720e0, 5.46378491116411436990e0, 1.78482653991729133580e0,
                             2.96560571828504891230e-1, 2.65321895265761230930e-2, 1.24266094738807843860e-3,
                             2.71155556874348757815e-5, 2.01033439929228813265e-7};
    constexpr double f[8] = {1.0, 5.99832206555887937690e-1, 1.36929880922735805310e-1, 1.48753612908506148525e-2,
                             7.86869131145613259100e-4, 1.84631831751005468180e-5, 1.42151175831644588870e-7,
                             2.04426310338993978564e-15};

    if (!(p > 0.0 && p < 1.0)) {
        if (p == 0.0) return -std::numeric_limits<double>::infinity();
        if (p == 1.0) return std::numeric_limits<double>::infinity();
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double q = p - 0.5;
    if (std::abs(q) <= 0.425) {
        const double r = 0.180625 - q * q;
        return q * detail::poly7(a, r) / detail::poly7(b, r);
    }

    // Tails: work with the smaller of p and 1 - p, so the left tail keeps full relative precision.
    double r = std::sqrt(-std::log(q < 0.0 ? p : 1.0 - p));
    double x;
    if (r <= 5.0) {
        r -= 1.6;
        x = detail::poly7(c, r) / detail::poly7(d, r);
    } else {
        r -= 5.0;
        x = detail::poly7(e, r) / detail::poly7(f, r);
    }
    return q < 0.0 ? -x : x;
}

} // namespace riskengine
