"""Jacobsthal's function and the quadratic bound ``h(k) << k^2`` (family 021).

For a positive integer ``n``, ``j(n)`` is the least ``m`` such that every interval
of ``m`` consecutive integers contains an integer coprime to ``n``.  Since
replacing ``n`` by its radical does not change ``j(n)``, the interesting object is

    h(k) = sup { j(n) : omega(n) <= k },

where ``omega(n)`` is the number of distinct prime divisors.  ``h(k) - 1`` is the
greatest length of a run of consecutive integers that can be covered by the
divisibility classes of at most ``k`` primes.

The manuscript *A quadratic bound for Jacobsthal's function* proves

    h(k) <= C k^2 / (log log (3k))^2      for every k >= 1,

with ``C`` absolute, improving Iwaniec's ``h(k) << k^2 log^2 k`` and answering
Jacobsthal's question whether ``h(k) << k^2``.

Only lower bounds on ``h`` are computable.  For a set of primes the value ``j(n)``
is found by sieving the multiples of those primes over one period of their
product and taking the longest covered run, and ``h(k)`` is then bounded below by
maximising over the prime sets that are searched.
"""

from __future__ import annotations

import numpy as np

from .sieves import primes_upto

#: Largest period, in positions, that :func:`jacobsthal_exact` will sieve directly.
#: A period is the product of the chosen primes, which grows like ``e^{k log k}``,
#: so this caps the reachable ``k`` at about 10 for the primorial case.
MAX_PERIOD = 1 << 34

__all__ = [
    "jacobsthal_exact",
    "largest_reachable_k",
    "MAX_PERIOD",
    "longest_covered_run",
    "primorial_jacobsthal",
    "h_lower_bound",
    "quadratic_envelope",
    "scan_shift_range",
]


def longest_covered_run(covered: np.ndarray) -> int:
    """Length of the longest run of True values in ``covered``.

    A run may wrap around the end of the array.  For a period of a product of
    primes this cannot happen: position ``0`` is covered by every prime and
    position ``P - 1`` is covered by none, so the wrap case is handled only for
    generality.
    """
    if covered.size == 0:
        return 0
    if covered.all():
        return int(covered.size)
    if not covered.any():
        return 0
    padded = np.concatenate(([False], covered, [False]))
    diff = np.diff(padded.astype(np.int8))
    starts = np.flatnonzero(diff == 1)
    ends = np.flatnonzero(diff == -1)
    best = int((ends - starts).max())
    if covered[0] and covered[-1]:
        head = int(np.argmin(covered))
        tail = int(covered.size - 1 - np.argmin(covered[::-1]))
        best = max(best, head + (covered.size - 1 - tail))
    return best


def _max_run_from_segments(primes: list[int], product: int, block: int) -> int:
    """Longest covered run over ``[0, product)`` using bounded memory.

    The range is walked in blocks of ``block`` positions.  ``position`` is always
    a multiple of the block size, so the multiples of each prime inside a block
    start at a computable offset; only the current run length has to be carried
    between blocks.
    """
    best = 0
    current = 0
    for start in range(0, product, block):
        length = min(block, product - start)
        covered = np.zeros(length, dtype=bool)
        for p in primes:
            p = int(p)
            first = (-start) % p
            if first < length:
                covered[first::p] = True
        if covered.all():
            current += length
            best = max(best, current)
            continue
        idx = np.flatnonzero(~covered)
        # run of Trues before the first gap
        current += int(idx[0])
        best = max(best, current)
        # runs strictly inside the block
        padded = np.concatenate(([False], covered, [False]))
        diff = np.diff(padded.astype(np.int8))
        s_in = np.flatnonzero(diff == 1)
        e_in = np.flatnonzero(diff == -1)
        if s_in.size:
            best = max(best, int((e_in - s_in).max()))
        # run after the last gap continues into the next block
        current = length - 1 - int(idx[-1])
        best = max(best, current)
    return best


def jacobsthal_exact(
    bound: int, primes: list[int], block: int = 1 << 26, max_period: int = MAX_PERIOD
) -> int:
    """Compute ``j(prod primes)`` for the given prime set.

    The multiples of a set of primes are periodic with period equal to their
    product, so it suffices to sieve one period.  ``bound`` must be at least that
    product.  Periods beyond ``block`` positions are walked in blocks by
    :func:`_max_run_from_segments`; periods beyond ``max_period`` raise, because
    the sieving work grows linearly in the period and is not worth attempting.

    Reachability is the real limit here: the period of the first ``k`` primes is
    ``e^{(1+o(1)) k log k}``, so ``k = 10`` already needs ``6.5e9`` positions and
    ``k = 11`` about ``2e11``.  Values of ``h(k)`` for larger ``k`` are taken from
    the literature instead of recomputed.
    """
    if not primes:
        return 1
    product = 1
    for p in primes:
        product *= int(p)
    if product > bound:
        raise ValueError(f"bound {bound} is smaller than the period {product}")
    if product > max_period:
        raise ValueError(
            f"period {product:.3e} exceeds max_period {max_period:.3e}; "
            "the run cannot be computed at this k"
        )
    if product > block:
        run = _max_run_from_segments(primes, product, block)
    else:
        covered = np.zeros(product, dtype=bool)
        for p in primes:
            covered[0 :: int(p)] = True
        run = longest_covered_run(covered)
    # j(n) is the least m for which every m consecutive integers contain a unit,
    # i.e. one more than the longest run of non-units.  The period P contains the
    # run [1, run], and every window of length run + 1 meets a unit.
    return run + 1


def primorial_jacobsthal(k: int, max_period: int = MAX_PERIOD) -> tuple[int, list[int]]:
    """``j(P_k)`` for ``P_k`` the product of the first ``k`` primes.

    Returns ``(value, primes)``.  This is the case Jacobsthal conjectured to be
    extremal, so it is a lower bound for ``h(k)``.  It is known to be strict at
    ``k = 24``, where ``j(P_24) = 234`` while ``h(24) = 236``, so the primorial
    value is never claimed here to equal ``h(k)``.
    """
    if k < 1:
        raise ValueError("k must be at least 1")
    ps = [int(p) for p in primes_upto(_kth_prime_bound(k))[:k]]
    if len(ps) < k:  # pragma: no cover - defensive
        raise RuntimeError("prime bound too small")
    return jacobsthal_exact(_period_of(ps), ps, max_period=max_period), ps


def largest_reachable_k(max_period: int = MAX_PERIOD) -> int:
    """The largest ``k`` whose primorial period fits under ``max_period``."""
    product = 1
    k = 0
    for p in primes_upto(1000):
        if product * int(p) > max_period:
            break
        product *= int(p)
        k += 1
    return k


def h_lower_bound(k: int, extra_primes: int = 0, max_period: int = MAX_PERIOD) -> int:
    """Lower bound for ``h(k)`` from the first ``k`` primes and nearby variants.

    Adding primes can only lengthen the longest covered run, and ``h(k)`` is a
    supremum over sets of at most ``k`` primes.  The family searched is the first
    ``k`` primes together with the sets obtained by replacing one of the largest
    members with one of the next ``extra_primes`` primes.

    The reported value is therefore a genuine lower bound for ``h(k)`` and, for the
    ``k`` reachable here, is expected to be attained by the primorial set.
    """
    if k < 1:
        raise ValueError("k must be at least 1")
    pool = [
        int(p)
        for p in primes_upto(_kth_prime_bound(k + extra_primes + 6))[: k + extra_primes]
    ]
    if len(pool) < k:
        raise RuntimeError("prime bound too small")
    best = jacobsthal_exact(_period_of(pool[:k]), pool[:k], max_period=max_period)
    if extra_primes:
        for i in range(k - 1, max(k - extra_primes - 1, 0) - 1, -1):
            for j in range(k, len(pool)):
                cand = list(pool[:k])
                cand[i] = pool[j]
                if len(set(cand)) != k:
                    continue
                try:
                    best = max(
                        best, jacobsthal_exact(_period_of(cand), cand, max_period=max_period)
                    )
                except ValueError:
                    continue
    return best


def scan_shift_range(primes: list[int], span: int) -> np.ndarray:
    """For each offset in ``[0, span)``, the covered run length starting there.

    This is diagnostic: it shows how the extremal run sits inside the period of
    the prime product.
    """
    product = _period_of(primes)
    covered = np.zeros(product, dtype=bool)
    for p in primes:
        covered[0::int(p)] = True
    span = min(span, product)
    out = np.empty(span, dtype=np.int32)
    n = covered.size
    for start in range(span):
        length = 0
        while covered[(start + length) % n]:
            length += 1
            if length == n:
                break
        out[start] = length
    return out


def quadratic_envelope(k: np.ndarray | float) -> np.ndarray:
    """The manuscript's bound shape ``k^2 / (log log (3k))^2``, up to the constant.

    All logarithms are natural.  The value is defined for ``k >= 1``; at ``k = 1``
    the double logarithm is negative, so the envelope is only used for ``k >= 3``
    in the probes.
    """
    k_arr = np.atleast_1d(np.asarray(k, dtype=float))
    if np.any(k_arr < 3.0):
        raise ValueError("the envelope is only meaningful for k >= 3")
    return k_arr**2 / np.log(np.log(3.0 * k_arr)) ** 2


def iwaniec_envelope(k: np.ndarray | float) -> np.ndarray:
    """Iwaniec's previous bound shape ``k^2 log^2 k``, for comparison."""
    k_arr = np.atleast_1d(np.asarray(k, dtype=float))
    return k_arr**2 * np.log(k_arr) ** 2


def _period_of(primes: list[int]) -> int:
    product = 1
    for p in primes:
        product *= int(p)
    return product


def _kth_prime_bound(k: int) -> int:
    """A bound for the ``k``-th prime, valid for ``k >= 6`` (Rosser's bound)."""
    if k < 6:
        return 16
    return int(k * (np.log(k) + np.log(np.log(k)))) + 8
