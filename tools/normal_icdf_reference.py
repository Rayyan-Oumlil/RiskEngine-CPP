"""High-precision inverse normal CDF reference values for tests/test_rng.cpp.

Each p is the exact double given by its Python literal; x = N^{-1}(p) is found by bisection at
60 significant digits (250 halvings of [-40, 40], error < 1e-70), independently of AS241.

Usage: python3 tools/normal_icdf_reference.py   (requires: pip install mpmath)
"""
from mpmath import mp, mpf, erfc, sqrt, nstr

mp.dps = 60

# Includes the AS241 region boundaries: |p - 0.5| = 0.425 and sqrt(-log p) = 5 (p ~ 1.39e-11).
PROBABILITIES = [1e-300, 1e-100, 1e-20, 1e-11, 2e-11, 1e-10, 1e-5, 0.001, 0.02, 0.075,
                 0.07500000000000001, 0.1, 0.3, 0.45, 0.5, 0.55, 0.8, 0.925, 0.975, 0.999,
                 1 - 1e-10]


def ncdf(x):
    return erfc(-x / sqrt(2)) / 2


def icdf(p):
    p = mpf(p)
    lo, hi = mpf(-40), mpf(40)
    for _ in range(250):
        mid = (lo + hi) / 2
        if ncdf(mid) < p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


for p in PROBABILITIES:
    x = "0.0" if p == 0.5 else nstr(icdf(p), 20)  # the median is exactly 0
    print(f"    {{{p!r}, {x}}},")
