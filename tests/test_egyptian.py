"""Egyptian fraction search: validity, minimality, and the greedy comparison."""

from __future__ import annotations

from fractions import Fraction

import numpy as np
import pytest

from openmath.egyptian import (
    double_log_envelope,
    erdos_graham_envelope,
    greedy_expansion,
    greedy_expansion_exact,
    greedy_length_table,
    max_min_length,
    max_min_length_table,
    shortest_expansion,
    vose_envelope,
)


def total(terms):
    return sum(Fraction(1, n) for n in terms)


def test_greedy_expansion_is_valid_for_many_fractions():
    for b in range(2, 120):
        for a in range(1, b):
            terms = greedy_expansion(a, b)
            assert total(terms) == Fraction(a, b), (a, b, terms)
            assert len(set(terms)) == len(terms)
            assert terms == sorted(terms)


def test_greedy_matches_exact_rational_implementation():
    for b in range(2, 60):
        for a in range(1, b):
            assert greedy_expansion(a, b) == greedy_expansion_exact(a, b)


def test_shortest_expansion_is_valid_and_no_shorter_exists():
    cases = {1: 1, 2: 2, 3: 3, 4: 3, 5: 3, 6: 3, 7: 4, 8: 4}
    for b, known in cases.items():
        for a in range(1, b):
            got = shortest_expansion(a, b, max_terms=6)
            assert got is not None
            assert total(got) == Fraction(a, b)
            assert len(set(got)) == len(got)
            # no expansion with fewer terms exists
            if len(got) > 1:
                assert shortest_expansion(a, b, max_terms=len(got) - 1) is None


def test_shortest_expansion_known_examples():
    # 5/121 needs three terms; the greedy expansion needs five
    assert len(shortest_expansion(5, 121, max_terms=6)) == 3
    assert len(greedy_expansion(5, 121)) == 5
    assert len(shortest_expansion(1, 7, max_terms=4)) == 1


def test_shortest_expansion_returns_none_within_a_tight_cap():
    assert shortest_expansion(6, 7, max_terms=1) is None


def test_max_min_length_is_the_maximum_over_numerators():
    b = 13
    best = max(len(shortest_expansion(a, b, max_terms=6)) for a in range(1, b))
    value, arg = max_min_length(b, max_terms=6)
    assert value == best
    assert 1 <= arg < b
    assert len(shortest_expansion(arg, b, max_terms=6)) == value


def test_max_min_length_table_is_consistent():
    bs, table = max_min_length_table(30, max_terms=5)
    assert bs.tolist() == list(range(2, 31))
    values = table[0]
    assert np.all(values >= 1)
    # N(b) grows slowly: at most five terms for b <= 30
    assert int(values.max()) <= 5


def test_greedy_is_never_shorter_than_the_minimum():
    bs, greedy = greedy_length_table(60, step=7)
    for b, g in zip(bs, greedy):
        exact, _ = max_min_length(int(b), max_terms=6)
        assert g >= exact


def test_envelopes_ordering_at_moderate_b():
    b = np.array([50.0, 100.0, 200.0])
    dlog = double_log_envelope(b)
    eg = erdos_graham_envelope(b)
    vo = vose_envelope(b)
    # The conjectured order dominates neither of the earlier bounds pointwise at
    # small b; only the double-logarithmic envelope is the smallest of the three
    # throughout this range, which is the fact the figure relies on.
    assert np.all(dlog < vo) and np.all(dlog < eg)
    assert np.all(vo < eg)


def test_input_validation():
    with pytest.raises(ValueError):
        greedy_expansion(0, 5)
    with pytest.raises(ValueError):
        greedy_expansion(5, 5)
    with pytest.raises(ValueError):
        shortest_expansion(5, 5)
    with pytest.raises(ValueError):
        max_min_length(1)
    with pytest.raises(ValueError):
        double_log_envelope(np.array([1.0]))
