"""Independent high-precision Black-Scholes reference values for tests/test_black_scholes.cpp.

Prices are evaluated with mpmath at 50 significant digits. Greeks are *numerical* derivatives
of that price (mpmath.diff), not the closed-form Greek formulas, so they independently check
the C++ closed forms, sign conventions and dividend handling.

Usage: python3 tools/bs_reference.py   (requires: pip install mpmath)
"""
from mpmath import mp, mpf, erfc, exp, log, sqrt, diff

mp.dps = 50


def ncdf(x):
    return erfc(-x / sqrt(2)) / 2


def price(kind, S, K, r, q, sigma, T):
    d1 = (log(S / K) + (r - q + sigma**2 / 2) * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    if kind == "call":
        return S * exp(-q * T) * ncdf(d1) - K * exp(-r * T) * ncdf(d2)
    return K * exp(-r * T) * ncdf(-d2) - S * exp(-q * T) * ncdf(-d1)


CASES = {
    # name: (S, K, r, q, sigma, T)
    "canonical": ("100", "100", "0.05", "0", "0.2", "1"),
    "dividend": ("100", "95", "0.03", "0.02", "0.25", "0.5"),
    "hull_15_6": ("42", "40", "0.1", "0", "0.2", "0.5"),
    # Tails: the out-of-the-money leg is ~1e-12 of the spot; tests its relative precision.
    "wing_k200": ("100", "200", "0.05", "0", "0.2", "0.25"),
    "wing_k50": ("100", "50", "0.05", "0", "0.2", "0.25"),
}

for name, args in CASES.items():
    S, K, r, q, sigma, T = map(mpf, args)
    print(f"// {name}: S={args[0]} K={args[1]} r={args[2]} q={args[3]} sigma={args[4]} T={args[5]}")
    for kind in ("call", "put"):
        f = lambda s=S, k=K, rr=r, qq=q, v=sigma, t=T: price(kind, s, k, rr, qq, v, t)
        vals = {
            "price": f(),
            "delta": diff(lambda x: f(s=x), S),
            "gamma": diff(lambda x: f(s=x), S, 2),
            "vega": diff(lambda x: f(v=x), sigma),
            "theta": -diff(lambda x: f(t=x), T),
            "rho": diff(lambda x: f(rr=x), r),
        }
        print(f"//   {kind}: " + ", ".join(f"{k}={mp.nstr(v, 17)}" for k, v in vals.items()))
