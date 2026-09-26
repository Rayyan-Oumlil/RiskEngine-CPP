"""Tests of the Python bindings against the engine's own reference values (docs/model_risk_report.md,
Appendix C). The bindings call the C++ engine directly, so these are the C++ numbers."""

import math

import numpy as np
import pytest

import riskengine as re

CANONICAL = dict(spot=100.0, strike=100.0, maturity=1.0, vol=0.2, rate=0.05)


def test_black_scholes_reference_values():
    assert re.black_scholes(**CANONICAL) == pytest.approx(10.450584, abs=1e-6)
    assert re.black_scholes(**CANONICAL, kind="put") == pytest.approx(5.573526, abs=1e-6)


def test_put_call_parity():
    call = re.black_scholes(**CANONICAL, dividend=0.02)
    put = re.black_scholes(**CANONICAL, dividend=0.02, kind="put")
    forward = 100.0 * math.exp(-0.02) - 100.0 * math.exp(-0.05)
    assert call - put == pytest.approx(forward, abs=1e-12)


def test_greeks_reference_values():
    g = re.greeks(**CANONICAL)
    assert g["delta"] == pytest.approx(0.636831, abs=1e-6)
    assert g["gamma"] == pytest.approx(0.018762, abs=1e-6)
    assert g["vega"] == pytest.approx(37.524035, abs=1e-6)
    assert g["rho"] == pytest.approx(53.232482, abs=1e-6)


def test_implied_vol_round_trip_and_arbitrage():
    price = re.black_scholes(**CANONICAL)
    assert re.implied_vol(price, 100.0, 100.0, 1.0, rate=0.05) == pytest.approx(0.2, abs=1e-12)
    with pytest.raises(ValueError):  # below the discounted intrinsic value: no vol exists
        re.implied_vol(1.0, 150.0, 100.0, 1.0, rate=0.05)


def test_trees():
    assert re.binomial(**CANONICAL, steps=1001, method="leisen_reimer") == pytest.approx(10.450584, abs=1e-6)
    american = re.binomial(**CANONICAL, kind="put", exercise="american", steps=20000)
    assert american == pytest.approx(6.090333, abs=5e-7)
    with pytest.raises(ValueError):  # Leisen-Reimer needs an odd number of steps
        re.binomial(**CANONICAL, steps=100, method="leisen_reimer")


def test_monte_carlo_is_accurate_and_thread_invariant():
    exact = re.black_scholes(**CANONICAL)
    one = re.monte_carlo(**CANONICAL, paths=1 << 18, seed=7, threads=1)
    assert abs(one.value - exact) < 4 * one.std_error
    for threads in (2, 3, 8):
        many = re.monte_carlo(**CANONICAL, paths=1 << 18, seed=7, threads=threads)
        assert (many.value, many.std_error) == (one.value, one.std_error)  # same bits, not approx


def test_longstaff_schwartz():
    one = re.longstaff_schwartz(**CANONICAL, paths=1 << 16, threads=1)
    assert abs(one.value - 6.090333) < 4 * one.std_error
    many = re.longstaff_schwartz(**CANONICAL, paths=1 << 16, threads=4)
    assert (many.value, many.std_error, many.discretization) == (one.value, one.std_error, one.discretization)


def test_model_risk_prices():
    heston = re.heston(100.0, 100.0, 1.0, 0.05, 0.0, v0=0.04, kappa=2.0, theta=0.04, xi=0.3, rho=-0.7)
    merton = re.merton(100.0, 100.0, 1.0, 0.2, 0.05, 0.0, jump_intensity=0.25, jump_mean=-0.2, jump_sd=0.1)
    assert heston == pytest.approx(10.394219, abs=1e-6)
    assert merton == pytest.approx(11.393999, abs=1e-6)
    with pytest.raises(ValueError):  # |rho| > 1
        re.heston(100.0, 100.0, 1.0, 0.05, 0.0, v0=0.04, kappa=2.0, theta=0.04, xi=0.3, rho=-1.5)


def test_normal_var_es():
    r = re.normal_var_es(0.0, 1.0, 0.99)
    assert r["var"] == pytest.approx(2.3263, abs=1e-4)
    assert r["es"] == pytest.approx(2.6652, abs=1e-4)


def test_normals_are_a_reproducible_numpy_array():
    z = re.normals(1_000_000, seed=2026)
    assert isinstance(z, np.ndarray) and z.dtype == np.float64 and z.shape == (1_000_000,)
    assert np.array_equal(z, re.normals(1_000_000, seed=2026))  # same seed, same numbers
    assert not np.array_equal(z, re.normals(1_000_000, seed=2027))
    assert abs(z.mean()) < 5e-3 and abs(z.std() - 1.0) < 5e-3
    assert np.array_equal(re.normals(10, seed=1, stream=3)[:4], re.normals(4, seed=1, stream=3))  # prefix of the stream


@pytest.mark.parametrize(
    "call",
    [
        lambda: re.black_scholes(**CANONICAL, kind="straddle"),
        lambda: re.black_scholes(spot=-1.0, strike=100.0, maturity=1.0, vol=0.2),
        lambda: re.monte_carlo(**CANONICAL, paths=1),
        lambda: re.monte_carlo(**CANONICAL, paths=1001, antithetic=True),
        lambda: re.binomial(**CANONICAL, exercise="bermudan"),
        lambda: re.normal_var_es(0.0, 1.0, confidence=1.5),
        lambda: re.normals(10, seed=1, stream=1 << 31),
    ],
)
def test_bad_inputs_raise_value_error(call):
    with pytest.raises(ValueError):
        call()
