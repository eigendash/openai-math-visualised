"""The Dickman function and the PD(1) sampler, against published values."""

from __future__ import annotations

import math

import numpy as np
import pytest

from openmath.poisson_dirichlet import (
    MAX_DICKMAN_U,
    dickman,
    normalized_log_factors,
    pd1_order_statistic_cdf,
    pd1_order_statistic_mean,
    predecessor_factors,
    stick_breaking_pd1,
)

# Published values of the Dickman--de Bruijn function (Dickman 1930).
RHO = {
    1.0: 1.0,
    1.5: 1.0 - math.log(1.5),
    2.0: 0.3068528194400547,
    2.5: 0.1303195618322514,
    3.0: 0.04860838829113157,
}


@pytest.mark.parametrize("u,expected", sorted(RHO.items()))
def test_dickman_matches_published_values(u, expected):
    assert dickman(u) == pytest.approx(expected, abs=1e-9)


def test_dickman_elementary_pieces_are_exact():
    assert dickman(0.0) == 1.0
    assert dickman(1.0) == 1.0
    for u in (1.1, 1.25, 1.5, 1.75, 1.999):
        assert dickman(u) == pytest.approx(1.0 - math.log(u), rel=0, abs=1e-15)


def test_dickman_is_decreasing_after_one():
    grid = np.linspace(1.0, 3.0, 400)
    values = dickman(grid)
    assert np.all(np.diff(values) < 0)


def test_dickman_satisfies_the_integrated_delay_relation():
    """Check the delay relation in the form that holds for the standard rho.

    For ``rho = 1`` on ``[0, 1]`` and ``u rho'(u) = -rho(u-1)``, multiplying by
    ``du`` and integrating gives ``d(u rho) = rho(u-1) du``, hence

        u rho(u) = int_{u-1}^{u} rho(s) ds      for u >= 2,

    since ``rho`` vanishes below ``u = 1``.  The solver builds its table from this
    relation, so the check validates the quadrature rather than independently
    confirming the function; the published-value test does the latter.
    """
    for u in (2.2, 2.5, 2.7, 3.0, 3.5):
        grid = np.linspace(u - 1.0, u, 200_001)
        integral = np.trapezoid(dickman(grid), grid)
        assert u * float(np.atleast_1d(dickman(u))[0]) == pytest.approx(integral, abs=5e-5)


def test_dickman_matches_an_independent_euler_integration():
    """A bare-bones forward Euler solve of the delay ODE, written independently.

    This shares no code with :func:`dickman`: it marches ``rho`` from ``u = 1`` with
    a small fixed step and reads the lag by linear interpolation, which is the
    textbook way to solve a delay equation numerically.
    """
    h = 1e-3
    u = np.arange(0.0, 3.0 + h, h)
    r = np.ones_like(u)
    start = int(round(1.0 / h))
    for i in range(start, len(u) - 1):
        t = u[i]
        lag = 1.0 if t - 1.0 <= 1.0 else float(np.interp(t - 1.0, u[: i + 1], r[: i + 1]))
        r[i + 1] = r[i] - h * lag / t
    for target in (2.0, 2.5, 3.0):
        euler = r[int(round(target / h))]
        # forward Euler at step 1e-3 carries an O(h) error of a few 1e-4
        assert float(np.atleast_1d(dickman(target))[0]) == pytest.approx(euler, abs=2e-3)


def test_dickman_right_hand_limit_is_below_the_left():
    """rho is decreasing on [1, inf), so its right limit at each integer is smaller."""
    for n in (1.0, 2.0, 3.0):
        left = float(np.atleast_1d(dickman(n - 1e-9))[0])
        right = float(np.atleast_1d(dickman(n + 1e-9))[0])
        assert right < left


def test_dickman_vectorises_and_rejects_out_of_range():
    out = dickman([0.5, 1.0, 2.0, 3.0])
    assert out.shape == (4,)
    assert np.isscalar(dickman(2.0)) or np.ndim(dickman(2.0)) == 0
    with pytest.raises(ValueError):
        dickman(-0.1)
    with pytest.raises(ValueError):
        dickman(MAX_DICKMAN_U + 1.0)


def test_stick_breaking_components_are_ordered_and_sum_below_one():
    draws = stick_breaking_pd1(4_000, top=8, seed=7)
    assert draws.shape == (4_000, 8)
    assert np.all(np.diff(draws, axis=1) <= 0)  # decreasing across columns
    assert np.all(draws >= 0.0)
    total = draws.sum(axis=1)
    assert np.all(total <= 1.0)
    # the discarded tail is prod_{i<=top} U_i, whose mean is 2^-top
    assert float(np.mean(1.0 - total)) == pytest.approx(2.0**-8, rel=0.05)


def test_pd1_largest_component_has_the_dickman_law():
    draws = stick_breaking_pd1(200_000, top=8, seed=3)
    for x in (0.2, 0.35, 0.5, 0.7, 0.9):
        empirical = float(np.mean(draws[:, 0] <= x))
        assert empirical == pytest.approx(float(dickman(1.0 / x)), abs=0.004), x


def test_pd1_cdf_agrees_with_the_sampler_for_all_components():
    draws = stick_breaking_pd1(100_000, top=8, seed=11)
    for j in (1, 2, 3):
        for x in (0.1, 0.25, 0.5):
            sampled = float(np.mean(draws[:, j - 1] <= x))
            exact = float(pd1_order_statistic_cdf(x, j, samples=100_000, seed=11))
            assert sampled == pytest.approx(exact, abs=0.01), (j, x)


def test_pd1_means_decrease_and_are_consistent_with_the_sampler():
    draws = stick_breaking_pd1(50_000, top=6, seed=5)
    means = draws.mean(axis=0)
    assert np.all(np.diff(means) < 0)
    # the helper draws with the standard depth; the leading columns agree with any
    # shallower draw made from the same seed
    deep = stick_breaking_pd1(50_000, top=12, seed=5)
    assert float(pd1_order_statistic_mean(1, samples=50_000, seed=5)) == pytest.approx(
        float(deep[:, 0].mean()), abs=1e-12
    )
    # a shallower draw uses the same leading uniforms but a different number of
    # columns, so only the ordering and the trend are comparable, not the values
    assert means[0] > means[1] > means[2]


def test_predecessor_factors_uses_multiplicity_and_sorts_descending():
    ps, facs = predecessor_factors(100)
    for p, f in zip(ps, facs):
        prod = 1
        for q in f:
            prod *= q
        assert prod == p - 1
        assert f == sorted(f, reverse=True)


def test_normalized_log_factors_sum_to_one_and_pad_with_zero():
    v = normalized_log_factors([2, 2, 2, 2, 3, 7], 1009, 4)
    assert v.shape == (4,)
    assert np.all(v > 0)
    # padding beyond the factor list contributes zeros and the total is exactly 1
    full = normalized_log_factors([2, 3], 7, 5)
    assert full[2:].tolist() == [0.0, 0.0, 0.0]
    assert float(full.sum()) == pytest.approx(1.0)
    # V_1(p) = log q_1 / log(p - 1), so it is 1 exactly when p - 1 is a prime power
    # whose largest prime factor equals p - 1: true for p = 3 (p - 1 = 2) and false
    # for p = 17 (p - 1 = 2^4, largest factor 2).
    assert normalized_log_factors([2], 3, 1)[0] == pytest.approx(1.0)
    assert normalized_log_factors([2, 2, 2, 2], 17, 1)[0] == pytest.approx(1.0 / 4.0)
