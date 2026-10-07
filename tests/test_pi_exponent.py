"""Rational approximation to pi: exact convergents and empirical exponent."""

from __future__ import annotations

import mpmath as mp
import numpy as np
import pytest

from openmath.pi_exponent import (
    approximation_table,
    best_approximation_quality,
    continued_fraction,
    convergents,
    convergents_exact,
    empirical_exponent,
    semiconvergents,
    to_float,
)


def test_continued_fraction_of_pi_prefix():
    cf = continued_fraction(mp.pi, 10, dps=60)
    assert cf == [3, 7, 15, 1, 292, 1, 1, 1, 2, 1]


def test_continued_fraction_of_rational_stops():
    # mpmath keeps exact rationals, so a rational's expansion terminates
    cf = continued_fraction(mp.mpf(8) / 5, 10, dps=50)
    assert cf[:2] == [1, 1]
    assert len(cf) <= 10


def test_continued_fraction_validates_terms():
    with pytest.raises(ValueError):
        continued_fraction(mp.pi, 0)


def test_convergents_exact_matches_float_version():
    cf = [3, 7, 15, 1, 292, 1, 1, 1, 2, 1]
    p_int, q_int = convergents_exact(cf)
    p_f, q_f = convergents(cf)
    assert p_int == [int(v) for v in p_f]
    assert q_int == [int(v) for v in q_f]
    assert p_int[:4] == [3, 22, 333, 355]
    assert q_int[:4] == [1, 7, 106, 113]


def test_convergents_exact_is_the_continued_fraction_recurrence():
    cf = [1, 2, 2, 2, 2]
    p, q = convergents_exact(cf)
    # approximants of sqrt(2)
    assert q == [1, 2, 5, 12, 29]
    assert p == [1, 3, 7, 17, 41]


def test_convergents_validate_input():
    with pytest.raises(ValueError):
        convergents([])
    with pytest.raises(ValueError):
        convergents_exact([])


def test_pi_convergents_are_the_best_approximations():
    table = approximation_table(mp.pi, terms=12, dps=100)
    qs = table["q"]
    errs = table["err"]
    for qi, ei in zip(qs, errs):
        if qi < 2:
            continue
        # no rational with denominator below q_i gets closer to pi
        limit = int(qi) - 1
        for q in range(2, min(limit, 300) + 1):
            p = int(mp.nint(mp.pi * q))
            with mp.workdps(100):
                e = abs(mp.pi - mp.mpf(p) / q)
            assert e >= ei, (qi, q)


def test_squared_quality_is_small_and_bounded():
    table = approximation_table(mp.pi, terms=25, dps=150)
    sq = table["squared_quality"]
    for i, (q, s) in enumerate(zip(table["q"], sq)):
        if q < 2:
            continue
        assert s < 1.0, i
        assert s > 0.0


def test_empirical_exponent_approaches_two_for_pi():
    table = approximation_table(mp.pi, terms=30, dps=150)
    keep = [i for i, q in enumerate(table["q"]) if q >= 2]
    exps = empirical_exponent([table["err"][i] for i in keep], [table["q"][i] for i in keep])
    # every convergent satisfies the exponent-2 bound, and the tail is close to 2
    assert np.all(exps > 1.9)
    assert float(np.median(exps[-10:])) == pytest.approx(2.0, abs=0.06)


def test_empirical_exponent_hand_checked():
    errors = np.array([1e-3, 1e-6])
    dens = np.array([10.0, 1000.0])
    got = empirical_exponent(errors, dens)
    assert got[0] == pytest.approx(3.0)
    assert got[1] == pytest.approx(2.0)


def test_empirical_exponent_validates():
    with pytest.raises(ValueError):
        empirical_exponent(np.array([0.0]), np.array([10.0]))
    with pytest.raises(ValueError):
        empirical_exponent(np.array([1e-3]), np.array([1.0]))
    with pytest.raises(ValueError):
        empirical_exponent(np.array([1e-3, 1e-4]), np.array([10.0]))


def test_semiconvergents_are_intermediants():
    cf = [3, 7, 15, 1]
    p_exact, q_exact = convergents_exact(cf)
    semis = semiconvergents(cf, 0)
    # index 0 is 3/1 with a_1 = 7, so the semiconvergents are t*3/1 + 1/0 = 3t
    # for t = 1, ..., a_1 - 1
    assert semis[0] == (3, 1)
    assert semis[-1] == (18, 6)
    assert len(semis) == 6
    # at index 1 (22/7) the semiconvergents run towards the next convergent 333/106
    semis1 = semiconvergents(cf, 1)
    # at index 1 the first semiconvergent is (p_1 + p_0)/(q_1 + q_0) = 25/8
    assert semis1[0] == (p_exact[1] + p_exact[0], q_exact[1] + q_exact[0])
    assert len(semis1) == cf[2] - 1
    # all of them lie strictly between the two convergents
    for p, q in semis1:
        assert q_exact[1] < q < q_exact[2]
    with pytest.raises(ValueError):
        semiconvergents(cf, 9)


def test_to_float_handles_large_integers_and_mpf():
    big = 10**30 + 7
    assert float(to_float([big])[0]) == float(big)
    assert float(to_float([mp.mpf("1e-30")])[0]) == pytest.approx(1e-30)


def test_best_approximation_quality_filters_q():
    q, sq = best_approximation_quality(mp.pi, q_max=1000)
    assert np.all(q <= 1000)
    assert np.all(q >= 2)
    assert sq.shape == q.shape
