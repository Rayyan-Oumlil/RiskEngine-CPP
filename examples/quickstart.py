"""Quickstart from Python: the same option priced every way, through the C++ engine.

    cmake -S . -B build-py -DCMAKE_BUILD_TYPE=Release -DRISKENGINE_BUILD_PYTHON=ON
    cmake --build build-py --target riskengine_python
    PYTHONPATH=build-py/python python examples/quickstart.py
"""

import time

import numpy as np

import riskengine as re

market = dict(spot=100.0, strike=100.0, maturity=1.0, vol=0.2, rate=0.05)

bs = re.black_scholes(**market)
print(f"engine built with: {re.simd_width()}")
print(f"Black-Scholes call        {bs:.6f}")
print(f"  greeks                  {re.greeks(**market)}")
print(f"  implied vol round trip  {re.implied_vol(bs, 100.0, 100.0, 1.0, rate=0.05):.12f}")
print(f"Leisen-Reimer, n = 1001   {re.binomial(**market, steps=1001, method='leisen_reimer'):.6f}")
print(f"American put, tree        {re.binomial(**market, kind='put', exercise='american', steps=2000, method='bbs_richardson'):.6f}")

lsm = re.longstaff_schwartz(**market, paths=1 << 17, threads=4)
print(f"American put, LSM         {lsm.value:.6f} +/- {lsm.std_error:.6f} (regression bias {lsm.discretization:+.4f})")

# Monte Carlo convergence: the standard error falls as N^(-1/2).
print("\nMonte Carlo call, 4 threads:")
for p in (14, 16, 18, 20):
    t0 = time.perf_counter()
    e = re.monte_carlo(**market, paths=1 << p, seed=p, threads=4)
    ms = 1e3 * (time.perf_counter() - t0)
    print(f"  N = 2^{p:<3} {e.value:.5f} +/- {e.std_error:.5f}   error {e.value - bs:+.5f}   {ms:6.1f} ms")

# The generator as a NumPy array (zero-copy): reproducible from its seed on any machine.
z = re.normals(1_000_000, seed=2026)
print(f"\n1M normals: mean {z.mean():+.5f}, sd {z.std():.5f}, P(|z| > 3) = {np.mean(np.abs(z) > 3):.5f} (theory 0.00270)")

print(f"\nHeston  {re.heston(100, 100, 1, 0.05, 0.0, v0=0.04, kappa=2.0, theta=0.04, xi=0.3, rho=-0.7):.6f}")
print(f"Merton  {re.merton(100, 100, 1, 0.2, 0.05, 0.0, jump_intensity=0.25, jump_mean=-0.2, jump_sd=0.1):.6f}")
