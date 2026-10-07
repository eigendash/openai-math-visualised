"""Zeta zeros and the fixed zero-free half-plane."""

from __future__ import annotations

import mpmath as mp
import numpy as np
import pytest

from openmath.quasi_riemann import (
    CLASSICAL_C,
    ZERO_FREE_BOUNDARY,
    classical_zero_free_region,
    densities_right_of,
    zero_count_check,
    zero_real_parts,
    zeros_respect_boundary,
    zeta_zeros,
)

# The first six ordinates of the nontrivial zeros of zeta.
FIRST_ORDINATES = [14.134725142, 21.022039639, 25.010857580, 30.424876126, 32.935061588, 37.586178159]


def test_zeta_zeros_match_the_published_ordinates():
    got = zeta_zeros(6)
    assert got.tolist() == pytest.approx(FIRST_ORDINATES, abs=1e-6)


def test_zeta_zeros_are_increasing_and_validate_count():
    ts = zeta_zeros(20)
    assert np.all(np.diff(ts) > 0)
    with pytest.raises(ValueError):
        zeta_zeros(0)


def test_zero_real_parts_are_one_half():
    re_parts = zero_real_parts(12)
    assert np.all(re_parts == pytest.approx(0.5, abs=1e-12))


def test_zeros_respect_the_seven_eighths_boundary():
    summary = zeros_respect_boundary(40)
    assert summary["violations"] == 0
    assert summary["max_real_part"] <= ZERO_FREE_BOUNDARY
    assert summary["slack"] == pytest.approx(ZERO_FREE_BOUNDARY - 0.5)


def test_boundary_constant_is_seven_eighths():
    assert ZERO_FREE_BOUNDARY == pytest.approx(0.875)


def test_zero_count_check_is_close_to_the_smooth_count():
    t_max = float(zeta_zeros(80)[-1])
    result = zero_count_check(t_max, 80)
    assert result["zeros_below_t_max"] == 80
    assert abs(result["deviation"]) <= 3.0


def test_siegel_count_grows_with_height():
    small = zero_count_check(float(zeta_zeros(5)[-1]), 5)
    large = zero_count_check(float(zeta_zeros(50)[-1]), 50)
    assert large["expected_smooth"] > small["expected_smooth"]


def test_classical_region_is_below_one_and_increases_with_t():
    ts = np.array([20.0, 100.0, 1000.0])
    sigma = classical_zero_free_region(ts, CLASSICAL_C)
    assert np.all(sigma < 1.0)
    assert np.all(np.diff(sigma) > 0)
    # at very large height the classical boundary approaches 1
    assert float(classical_zero_free_region(1e9, CLASSICAL_C)) > 0.99


def test_classical_region_validates_height():
    with pytest.raises(ValueError):
        classical_zero_free_region(np.array([1.0, 2.0]))


def test_densities_right_of_are_a_step_at_one_half():
    t_max = float(zeta_zeros(20)[-1])
    sigmas = np.array([0.3, 0.5, 0.7, 0.9])
    dens = densities_right_of(t_max, sigmas)
    assert dens.tolist() == pytest.approx([1.0, 0.0, 0.0, 0.0])
    with pytest.raises(ValueError):
        densities_right_of(t_max, np.array([1.5]))


def test_mpmath_cross_check_of_one_zero():
    # an independent check that mpmath's zero is a zero of zeta
    t = float(mp.im(mp.zetazero(1)))
    value = complex(mp.zeta(mp.mpf("0.5") + 1j * mp.mpf(str(t))))
    assert abs(value) < 1e-6
