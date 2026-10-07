"""Jacobsthal's function against the published primorial values."""

from __future__ import annotations

import numpy as np
import pytest

from openmath.jacobsthal import (
    MAX_PERIOD,
    h_lower_bound,
    iwaniec_envelope,
    jacobsthal_exact,
    largest_reachable_k,
    longest_covered_run,
    primorial_jacobsthal,
    quadratic_envelope,
)

# j(P_k) for the first ten primorials (OEIS A048670).
PRIMORIAL_J = {1: 2, 2: 4, 3: 6, 4: 10, 5: 14, 6: 22, 7: 26, 8: 34, 9: 40, 10: 46}


@pytest.mark.parametrize("k,expected", sorted(PRIMORIAL_J.items()))
def test_primorial_jacobsthal_matches_literature(k, expected):
    value, _ = primorial_jacobsthal(k)
    assert value == expected


def test_longest_covered_run_counts_runs_and_wrap():
    assert longest_covered_run(np.array([], dtype=bool)) == 0
    assert longest_covered_run(np.array([False, False])) == 0
    assert longest_covered_run(np.array([True, True])) == 2
    assert longest_covered_run(np.array([True, False, True, True, False])) == 2
    # wrap-around: the trailing and leading runs join
    assert longest_covered_run(np.array([True, False, True])) == 2
    assert longest_covered_run(np.array([True, True, False, True])) == 3


def test_jacobsthal_exact_is_the_run_plus_one():
    # For n = 2 the longest run of non-units is 1 (the single integer 2), so j = 2
    assert jacobsthal_exact(2, [2]) == 2
    # For n = 6 the period is 6 and the run [2,4] has length 3, so j = 4
    assert jacobsthal_exact(6, [2, 3]) == 4
    # For n = 30 the run [2, 6] has length 5, so j = 6
    assert jacobsthal_exact(30, [2, 3, 5]) == 6


def test_jacobsthal_exact_rejects_a_period_that_does_not_fit():
    with pytest.raises(ValueError):
        jacobsthal_exact(5, [2, 3])  # period 6 > 5


def test_segmented_and_direct_paths_agree():
    primes = [2, 3, 5, 7, 11]
    period = 2 * 3 * 5 * 7 * 11
    direct = jacobsthal_exact(period, primes, block=1 << 30)
    segmented = jacobsthal_exact(period, primes, block=997)
    assert direct == segmented == PRIMORIAL_J[5]


def test_h_lower_bound_is_at_least_the_primorial_value():
    for k in range(1, 8):
        assert h_lower_bound(k, extra_primes=0) >= PRIMORIAL_J[k]


def test_envelopes_and_reachability():
    k = np.array([3.0, 5.0, 10.0])
    assert np.all(quadratic_envelope(k) > 0)
    assert np.all(iwaniec_envelope(k) > 0)
    # The two envelopes are not ordered at every k: k^2 log^2 k exceeds
    # k^2/(log log 3k)^2 only once log log 3k > 1, i.e. k >~ 4.
    assert iwaniec_envelope(np.array([10.0]))[0] > quadratic_envelope(np.array([10.0]))[0]
    with pytest.raises(ValueError):
        quadratic_envelope(np.array([2.0]))
    assert largest_reachable_k() >= 10
    assert MAX_PERIOD >= 10**10
