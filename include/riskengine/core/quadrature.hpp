#pragma once

#include <array>
#include <cmath>
#include <stdexcept>

namespace riskengine {

// Integral over [0, inf) of a smooth integrand that decays at infinity, by 16-point Gauss-Legendre
// on consecutive panels of width `h`, until `quiet_panels` successive panels each add less than
// `tol` times the running total (plus tol). Used for characteristic-function prices (Heston), whose
// integrands are smooth and decay at least exponentially. Throws if `max_panels` are not enough.
template <class F>
double integrate_half_line(F&& f, double h, double tol = 1e-15, int quiet_panels = 4, int max_panels = 20000) {
    // Nodes and weights on [-1, 1], symmetric: only the positive half is stored (numpy leggauss(16)).
    static constexpr std::array<std::array<double, 2>, 8> kGauss = {{
        {0.09501250983763744, 0.18945061045506864},
        {0.2816035507792589, 0.18260341504492364},
        {0.45801677765722737, 0.16915651939500265},
        {0.6178762444026438, 0.1495959888165767},
        {0.755404408355003, 0.12462897125553407},
        {0.8656312023878318, 0.0951585116824926},
        {0.9445750230732326, 0.062253523938647456},
        {0.9894009349916499, 0.027152459411754176},
    }};
    double total = 0.0;
    int quiet = 0;
    for (int k = 0; k < max_panels; ++k) {
        const double mid = (k + 0.5) * h, half = 0.5 * h;
        double panel = 0.0;
        for (const auto& [x, w] : kGauss) panel += w * (f(mid - half * x) + f(mid + half * x));
        panel *= half;
        total += panel;
        quiet = std::abs(panel) < tol * (std::abs(total) + 1.0) ? quiet + 1 : 0;
        if (quiet == quiet_panels) return total;
    }
    throw std::runtime_error("integrate_half_line: integrand has not decayed");
}

} // namespace riskengine
