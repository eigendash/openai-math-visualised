"""Sieve primitives: correctness against brute force and against each other."""

from __future__ import annotations

import numpy as np
import pytest

from openmath.sieves import (
    factorize_from_lpf,
    lpf_sieve,
    prime_sieve,
    primes_upto,
    spf_sieve,
)


def brute_primes(limit: int) -> list[int]:
    return [n for n in range(2, limit + 1) if all(n % d for d in range(2, int(n**0.5) + 1))]


def test_prime_sieve_matches_brute_force():
    for limit in (0, 1, 2, 10, 97, 1000):
        got = np.flatnonzero(prime_sieve(limit)).tolist()
        assert got == brute_primes(limit), limit


def test_prime_sieve_known_prefixes():
    assert primes_upto(30).tolist() == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    assert primes_upto(1).tolist() == []
    assert int(primes_upto(10_000).size) == 1229  # pi(10^4)


def test_smallest_prime_factor():
    spf = spf_sieve(1000)
    assert spf[0] == 1 and spf[1] == 1
    for n in range(2, 1001):
        expected = next(p for p in brute_primes(n) if n % p == 0)
        assert int(spf[n]) == expected, n


def test_largest_prime_factor_known_values():
    lpf = lpf_sieve(100)
    assert lpf[0] == 0 and lpf[1] == 1
    assert int(lpf[2]) == 2
    assert int(lpf[97]) == 97
    assert int(lpf[100]) == 5
    assert int(lpf[96]) == 3
    # P+(n) <= n always, with equality exactly at the primes
    for n in range(2, 101):
        assert int(lpf[n]) <= n
        assert (int(lpf[n]) == n) == (n in brute_primes(n))


def test_factorize_from_lpf_is_non_increasing_and_correct():
    lpf = lpf_sieve(5000)
    for n in range(2, 5000):
        f = factorize_from_lpf(n, lpf)
        assert f == sorted(f, reverse=True)
        prod = 1
        for p in f:
            prod *= p
        assert prod == n
        for p in f:
            assert p in brute_primes(p)


def test_factorize_rejects_out_of_range():
    lpf = lpf_sieve(100)
    with pytest.raises(IndexError):
        factorize_from_lpf(101, lpf)
    assert factorize_from_lpf(1, lpf) == []


def test_sieves_validate_inputs():
    with pytest.raises(ValueError):
        prime_sieve(-1)
    with pytest.raises(ValueError):
        spf_sieve(0)
    with pytest.raises(ValueError):
        lpf_sieve(0)
