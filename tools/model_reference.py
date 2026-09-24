"""Reference prices for the Heston and Merton models (tests/test_models.cpp), independent of the C++.

Heston: the Gil-Pelaez integrals of the characteristic function in the Albrecher et al. (2007)
"little trap" form, integrated by mpmath's adaptive tanh-sinh quadrature at 30 digits (the C++ uses
fixed Gauss-Legendre panels in double precision). Merton: the Poisson series of Black-Scholes prices,
summed in 30-digit arithmetic until the remaining Poisson mass is below 1e-30.

Usage: python3 tools/model_reference.py
"""
import mpmath as mp

mp.mp.dps = 30


def heston_cf(u, log_fwd, t, v0, kappa, theta, xi, rho):
    i = mp.mpc(0, 1)
    beta = kappa - rho * xi * i * u
    d = mp.sqrt(beta**2 + xi**2 * (i * u + u**2))
    g = (beta - d) / (beta + d)
    e = mp.exp(-d * t)
    big_d = (beta - d) / xi**2 * (1 - e) / (1 - g * e)
    big_c = kappa * theta / xi**2 * ((beta - d) * t - 2 * mp.log((1 - g * e) / (1 - g)))
    return mp.exp(big_c + big_d * v0 + i * u * log_fwd)


def heston_call(s, k, r, q, t, v0, kappa, theta, xi, rho):
    s, k, r, q, t = map(mp.mpf, (s, k, r, q, t))
    i = mp.mpc(0, 1)
    log_fwd, log_k = mp.log(s) + (r - q) * t, mp.log(k)
    args = (t, v0, kappa, theta, xi, rho)
    phi_mi = heston_cf(-i, log_fwd, *args)
    p1 = mp.mpf(0.5) + mp.quad(lambda u: mp.re(mp.exp(-i * u * log_k) * heston_cf(u - i, log_fwd, *args)
                                               / (i * u * phi_mi)), [0, mp.inf]) / mp.pi
    p2 = mp.mpf(0.5) + mp.quad(lambda u: mp.re(mp.exp(-i * u * log_k) * heston_cf(u, log_fwd, *args)
                                               / (i * u)), [0, mp.inf]) / mp.pi
    return s * mp.exp(-q * t) * p1 - k * mp.exp(-r * t) * p2, mp.exp(-r * t) * p2


def bs_call(s, k, r, q, sigma, t):
    d1 = (mp.log(s / k) + (r - q + sigma**2 / 2) * t) / (sigma * mp.sqrt(t))
    d2 = d1 - sigma * mp.sqrt(t)
    n = lambda x: mp.erfc(-x / mp.sqrt(2)) / 2
    return s * mp.exp(-q * t) * n(d1) - k * mp.exp(-r * t) * n(d2)


def merton_call(s, k, r, q, sigma, t, lam, mu, delta):
    s, k, r, q, sigma, t, lam, mu, delta = map(mp.mpf, (s, k, r, q, sigma, t, lam, mu, delta))
    kbar = mp.exp(mu + delta**2 / 2) - 1
    intensity = lam * (1 + kbar) * t
    price, mass, n = mp.mpf(0), mp.mpf(0), 0
    while True:
        w = mp.exp(-intensity) * intensity**n / mp.factorial(n)
        sigma_n = mp.sqrt(sigma**2 + n * delta**2 / t)
        r_n = r - lam * kbar + n * mp.log(1 + kbar) / t
        price += w * bs_call(s, k, r_n, q, sigma_n, t)
        mass += w
        n += 1
        if n > intensity and 1 - mass < mp.mpf(10) ** -30:
            return price


HESTON = [  # s, k, r, q, t, v0, kappa, theta, xi, rho
    ("feller_atm", 100, 100, 0.05, 0.0, 1.0, 0.04, 2.0, 0.04, 0.3, -0.7),
    ("feller_k80", 100, 80, 0.05, 0.0, 1.0, 0.04, 2.0, 0.04, 0.3, -0.7),
    ("feller_k120", 100, 120, 0.05, 0.0, 1.0, 0.04, 2.0, 0.04, 0.3, -0.7),
    ("andersen_case1", 100, 100, 0.0, 0.0, 10.0, 0.04, 0.5, 0.04, 1.0, -0.9),
    ("short_dividend", 100, 95, 0.03, 0.02, 0.1, 0.09, 1.0, 0.06, 0.8, -0.5),
]
MERTON = [  # s, k, r, q, sigma, t, lambda, jump_mean, jump_sd
    ("crash_atm", 100, 100, 0.05, 0.0, 0.2, 1.0, 0.5, -0.1, 0.15),
    ("crash_k80", 100, 80, 0.05, 0.0, 0.2, 1.0, 0.5, -0.1, 0.15),
    ("frequent_short", 100, 110, 0.03, 0.01, 0.15, 0.25, 2.0, 0.05, 0.05),
]

if __name__ == "__main__":
    for name, *a in HESTON:
        call, digital = heston_call(*a)
        print(f'{{"{name}", {", ".join(repr(x) for x in a)}, {mp.nstr(call, 17)}, {mp.nstr(digital, 17)}}},')
    for name, *a in MERTON:
        print(f'{{"{name}", {", ".join(repr(x) for x in a)}, {mp.nstr(merton_call(*a), 17)}}},')
