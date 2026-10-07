"""Positive lower density of large prime gaps (family 026).

Let ``p_n`` be the ``n``-th prime and ``d_n = p_{n+1} - p_n``.  The manuscript
*Positive lower density of large prime gaps* proves that for every fixed real
``C > 0`` there are constants ``c(C) > 0`` and ``N_0(C)`` with

    #{1 <= n <= N : d_n > C log p_n} >= c(C) N      for all N >= N_0(C).

So each threshold ``C`` is exceeded by a positive proportion of gaps, with the
density measured by counting prime indices rather than by counting integers.  A
corollary is that ``{n : p_n / n < p_{n+1} / (n+1)}`` has positive lower density,
because that inequality is equivalent to ``d_n > p_n / n`` and ``p_n / n ~ log p_n``.

The probes here compute the gap sequence by sieving and measure the empirical
proportion for a range of thresholds, together with the density of the Erdős--
Prachar comparison.
"""

from __future__ import annotations

import numpy as np

from .sieves import prime_sieve

__all__ = [
    "prime_gaps",
    "exceedance_density",
    "exceedance_curve",
    "prachar_density",
    "cramer_reference",
]


def prime_gaps(limit: int) -> tuple[np.ndarray, np.ndarray]:
    """Primes at most ``limit`` and the consecutive gaps between them.

    Returns ``(primes, gaps)`` where ``gaps[i] = primes[i + 1] - primes[i]``, so
    ``gaps`` has length ``len(primes) - 1``.  With ``primes`` written ``p_n`` for
    ``n`` starting at 1, ``gaps[n - 1] = d_n``.
    """
    if limit < 3:
        raise ValueError("limit must be at least 3")
    ps = np.flatnonzero(prime_sieve(limit)).astype(np.int64)
    return ps, np.diff(ps)


def exceedance_density(
    primes: np.ndarray, gaps: np.ndarray, c: float
) -> tuple[float, int]:
    """Empirical proportion of gaps with ``d_n > c log p_n``.

    Returns ``(density, count)``.  ``primes`` must align with ``gaps`` so that
    ``primes[i]`` is the prime preceding ``gaps[i]``.
    """
    if primes.size != gaps.size + 1:
        raise ValueError("primes must have exactly one more entry than gaps")
    if c <= 0:
        raise ValueError("the threshold C must be positive")
    logp = np.log(primes[:-1].astype(float))
    mask = gaps > c * logp
    return float(mask.mean()), int(mask.sum())


def exceedance_curve(
    primes: np.ndarray, gaps: np.ndarray, thresholds: np.ndarray
) -> np.ndarray:
    """Empirical densities for a vector of thresholds ``C``."""
    return np.array([exceedance_density(primes, gaps, float(c))[0] for c in thresholds])


def prachar_density(primes: np.ndarray, gaps: np.ndarray) -> float:
    """Density of ``p_n / n < p_{n+1} / (n + 1)``, the Erdős--Prachar set.

    The comparison is equivalent to ``d_n > p_n / n``, so the density is measured
    directly from that inequality rather than from the ratios.
    """
    if primes.size != gaps.size + 1:
        raise ValueError("primes must have exactly one more entry than gaps")
    n = np.arange(1, gaps.size + 1, dtype=float)
    mask = gaps.astype(float) > primes[:-1].astype(float) / n
    return float(mask.mean())


def cramer_reference(thresholds: np.ndarray) -> np.ndarray:
    """The Cramér model's prediction ``exp(-C)`` for the exceedance density.

    Under the heuristic that prime gaps near ``p`` are approximately geometric with
    mean ``log p``, the probability that a gap exceeds ``C log p`` is ``exp(-C)``.
    This is a heuristic reference curve, not a theorem, and it is not a substitute
    for the positive-density statement: that statement asserts a positive
    *proportion* for every fixed ``C`` without identifying the value.
    """
    c = np.asarray(thresholds, dtype=float)
    if np.any(c < 0):
        raise ValueError("thresholds must be non-negative")
    return np.exp(-c)


def density_stability(
    primes: np.ndarray, gaps: np.ndarray, c: float, splits: int = 8
) -> np.ndarray:
    """Exceedance density evaluated on successive contiguous blocks of the data.

    The theorem asserts a positive proportion for every sufficiently large initial
    segment, so a genuine positive limiting density should be roughly stable across
    blocks, whereas a density driven by finitely many record gaps would decay.
    """
    if splits < 2:
        raise ValueError("splits must be at least 2")
    logp = np.log(primes[:-1].astype(float))
    mask = (gaps > c * logp).astype(float)
    blocks = np.array_split(mask, splits)
    return np.array([float(b.mean()) for b in blocks])
