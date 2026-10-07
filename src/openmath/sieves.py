"""Array-based sieves used by several of the probes.

The three primitives here are

* :func:`prime_sieve` -- a boolean sieve of Eratosthenes,
* :func:`spf_sieve`   -- the smallest prime factor of every integer in a range,
* :func:`lpf_sieve`   -- the largest prime factor of every integer in a range.

``spf`` and ``lpf`` are the two functions needed to factor a whole interval at
once rather than one integer at a time. ``lpf`` is what makes the prime-gap and
joint-Dickman probes possible at their stated ranges.
"""

from __future__ import annotations

import numpy as np

__all__ = ["prime_sieve", "primes_upto", "spf_sieve", "lpf_sieve", "factorize_from_lpf"]


def prime_sieve(limit: int) -> np.ndarray:
    """Return a boolean array ``is_prime`` of length ``limit + 1``.

    ``is_prime[n]`` is True exactly when ``n`` is prime. Sizes 0 and 1 are
    handled by the array length rather than by special cases.
    """
    if limit < 0:
        raise ValueError("limit must be non-negative")
    is_prime = np.ones(limit + 1, dtype=bool)
    if limit >= 0:
        is_prime[0] = False
    if limit >= 1:
        is_prime[1] = False
    for p in range(2, int(limit**0.5) + 1):
        if is_prime[p]:
            is_prime[p * p :: p] = False
    return is_prime


def primes_upto(limit: int) -> np.ndarray:
    """Return the primes at most ``limit`` as an ``int64`` array."""
    return np.flatnonzero(prime_sieve(limit)).astype(np.int64)


def spf_sieve(limit: int) -> np.ndarray:
    """Smallest prime factor of every integer in ``[0, limit]``.

    ``spf[0] = spf[1] = 1`` by convention, and ``spf[p] = p`` for prime ``p``.
    The array is ``int64`` so that it also serves as a table of the *distinct*
    prime factors of each index.
    """
    if limit < 1:
        raise ValueError("limit must be at least 1")
    spf = np.zeros(limit + 1, dtype=np.int64)
    spf[0] = 1
    if limit >= 1:
        spf[1] = 1
    for p in range(2, int(limit**0.5) + 1):
        if spf[p] == 0:  # p is prime
            block = spf[p * p :: p]
            block[block == 0] = p
    # Anything still zero is a prime larger than sqrt(limit).
    remaining = spf == 0
    remaining[0] = False
    spf[remaining] = np.flatnonzero(remaining)
    return spf


def lpf_sieve(limit: int) -> np.ndarray:
    """Largest prime factor of every integer in ``[0, limit]``.

    ``lpf[0] = 0`` and ``lpf[1] = 1`` by convention, so that
    ``lpf[n] = P^+(n)`` in the notation of the joint-Dickman manuscript.
    """
    if limit < 1:
        raise ValueError("limit must be at least 1")
    lpf = np.ones(limit + 1, dtype=np.int64)
    lpf[0] = 0
    for p in range(2, limit + 1):
        if lpf[p] == 1:  # p is prime, so it has not been written yet
            lpf[p::p] = p
    return lpf


def factorize_from_lpf(n: int, lpf: np.ndarray) -> list[int]:
    """Factor ``n`` into primes with multiplicity using a supplied ``lpf`` table.

    Returns the factors in non-increasing order, matching the ordering in the
    Poisson--Dirichlet manuscript, where the prime factors of ``p - 1`` are
    listed as ``q_1 >= q_2 >= ...``.
    """
    if n < 2:
        return []
    if n >= lpf.shape[0]:
        raise IndexError("n is outside the range covered by the lpf table")
    factors: list[int] = []
    while n > 1:
        p = int(lpf[n])
        while n % p == 0:
            factors.append(p)
            n //= p
    return factors
