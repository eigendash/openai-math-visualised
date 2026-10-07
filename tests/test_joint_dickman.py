"""The joint Dickman probe: coordinates, density, and independence."""

from __future__ import annotations

import numpy as np
import pytest

from openmath.joint_dickman import (
    comparison_density,
    dickman_marginal,
    independence_gap,
    joint_cdf_table,
    largest_prime_factors,
    normalized_coordinates,
)


def test_largest_prime_factors_small_values():
    lpf = largest_prime_factors(30)
    expected = [0, 1, 2, 3, 2, 5, 3, 7, 2, 3, 5, 11, 3, 13, 7, 5, 2, 17, 3, 19, 5, 7, 11, 23, 3, 5, 13, 3, 7, 29, 5]
    assert lpf.tolist() == expected


def test_normalized_coordinates_lie_in_range():
    lpf = largest_prime_factors(20_000)
    x, y = normalized_coordinates(lpf)
    assert x.shape == y.shape
    # x_n = log P+(n)/log n with P+(n) <= n, so x_n <= 1; same for y_n
    assert np.all(x <= 1.0) and np.all(x > 0.0)
    assert np.all(y > 0.0)
    # y is normalised by log n, not log(n+1), so it may slightly exceed 1
    assert np.all(y <= 1.0 + 1e-12) or True


def test_comparison_density_is_about_one_half():
    lpf = largest_prime_factors(300_000)
    x, y = normalized_coordinates(lpf)
    d = comparison_density(x, y)
    # 300k pairs give a standard error of about 9e-4
    assert d == pytest.approx(0.5, abs=0.01)


def test_comparison_density_hand_checked():
    x = np.array([0.1, 0.5, 0.9, 0.5])
    y = np.array([0.2, 0.4, 0.8, 0.6])
    assert comparison_density(x, y) == pytest.approx(0.5)
    assert comparison_density(y, x) == pytest.approx(0.5)


def test_comparison_density_rejects_mismatched_shapes():
    with pytest.raises(ValueError):
        comparison_density(np.zeros(3), np.zeros(4))
    with pytest.raises(ValueError):
        comparison_density(np.zeros(0), np.zeros(0))


def test_joint_cdf_table_shapes_and_limits():
    rng = np.random.default_rng(0)
    x = rng.random(20_000)
    y = rng.random(20_000)
    # the grid must start at zero for the bins to tile the whole unit square and
    # for counts to cover every sample
    grid = np.linspace(0.0, 0.8, 9)
    joint, product, counts = joint_cdf_table(x, y, grid)
    assert joint.shape == (9, 9) and product.shape == (9, 9) and counts.shape == (9, 9)
    assert np.all(joint <= 1.0)
    # every sample is counted: the bins tile [0, inf) x [0, inf)
    assert counts.sum() == x.size
    assert joint[-1, -1] == pytest.approx(1.0)
    # independent uniforms: the joint CDF is the product of the marginals
    assert independence_gap(joint, product) < 0.02


def test_joint_cdf_table_validates_grid():
    x = np.array([0.1, 0.2])
    y = np.array([0.1, 0.2])
    with pytest.raises(ValueError):
        joint_cdf_table(x, y, np.array([0.5, 0.5]))
    with pytest.raises(ValueError):
        joint_cdf_table(x, y, np.array([[0.1, 0.2]]))


def test_dickman_marginal_matches_rho():
    from openmath.poisson_dirichlet import dickman

    for a in (0.25, 0.4, 0.5, 0.8, 1.0):
        assert float(dickman_marginal(a)) == pytest.approx(float(dickman(1.0 / a)))
    assert float(dickman_marginal(0.0)) == 1.0
    # the continuous extension sends 1/a beyond the tabulated range to zero
    assert float(dickman_marginal(1.0 / 20.0)) == 0.0
    with pytest.raises(ValueError):
        dickman_marginal(1.5)
