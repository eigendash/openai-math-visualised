"""Prime gaps: exceedance density, stability, and the Prachar comparison."""

from __future__ import annotations

import numpy as np
import pytest

from openmath.prime_gaps import (
    cramer_reference,
    density_stability,
    exceedance_curve,
    exceedance_density,
    prachar_density,
    prime_gaps,
)


def test_prime_gaps_small():
    ps, gaps = prime_gaps(30)
    assert ps.tolist() == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    assert gaps.tolist() == [1, 2, 2, 4, 2, 4, 2, 4, 6]


def test_prime_gaps_alignment_invariant():
    ps, gaps = prime_gaps(10_000)
    assert ps.size == gaps.size + 1
    assert np.all(gaps > 0)


def test_exceedance_density_hand_checked():
    ps = np.array([2, 3, 5, 7, 11, 13])
    gaps = np.array([1, 2, 2, 4, 2])
    # log p_n for the five gaps: log 2, log 3, log 5, log 7, log 11
    # d_n > 1 * log p_n holds only for the gap 4 at p = 7 (log 7 = 1.946)
    density, count = exceedance_density(ps, gaps, 1.0)
    assert count == int(np.sum(gaps > np.log(ps[:-1])))
    assert density == pytest.approx(count / gaps.size)


def test_exceedance_density_validates():
    ps = np.array([2, 3, 5])
    gaps = np.array([1, 2])
    with pytest.raises(ValueError):
        exceedance_density(ps, np.array([1, 2, 3]), 1.0)
    with pytest.raises(ValueError):
        exceedance_density(ps, gaps, 0.0)


def test_exceedance_is_positive_and_decreasing_in_c():
    ps, gaps = prime_gaps(2_000_000)
    cs = np.array([0.25, 0.5, 1.0, 2.0, 3.0])
    dens = exceedance_curve(ps, gaps, cs)
    assert np.all(dens > 0.0)
    assert np.all(np.diff(dens) < 0.0)


def test_cramer_reference_shape():
    cs = np.array([0.0, 1.0, 2.0])
    assert cramer_reference(cs).tolist() == pytest.approx([1.0, np.exp(-1), np.exp(-2)])
    with pytest.raises(ValueError):
        cramer_reference(np.array([-1.0]))


def test_exceedance_density_is_stable_across_the_range():
    ps, gaps = prime_gaps(20_000_000)
    for c in (0.5, 1.0):
        blocks = density_stability(ps, gaps, c, splits=5)
        assert blocks.size == 5
        # positive lower density with no strong decay towards the end
        assert blocks.min() > 0.0
        assert blocks[-1] > 0.7 * blocks[0]


def test_prachar_density_is_positive():
    ps, gaps = prime_gaps(2_000_000)
    d = prachar_density(ps, gaps)
    assert 0.2 < d < 0.6
    with pytest.raises(ValueError):
        prachar_density(ps, gaps[:-1])


def test_density_stability_validates_splits():
    ps, gaps = prime_gaps(1000)
    with pytest.raises(ValueError):
        density_stability(ps, gaps, 1.0, splits=1)
