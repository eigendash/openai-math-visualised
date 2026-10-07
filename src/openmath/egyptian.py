"""Short Egyptian fractions and the double-logarithmic order of ``N(b)`` (family 025).

For integers ``1 <= a < b`` let ``N(a, b)`` be the least ``k`` for which

    a/b = 1/n_1 + ... + 1/n_k,   2 <= n_1 < ... < n_k,

with no upper bound on the denominators, and put ``N(b) = max_{1 <= a < b} N(a, b)``.
The manuscript *Short Egyptian fractions* proves that there are absolute constants
``c_1, c_2 > 0`` with

    c_1 log log b <= N(b) <= c_2 log log b      for all large b,

resolving Erdős's conjecture and improving Vose's ``N(b) << sqrt(log b)`` and the
classical ``N(b) << log b / log log b``.  The individual numerator can be far
cheaper: ``N(1, b) = 1`` for every ``b``.

``N(a, b)`` is computed here exactly for small ``b`` by depth-limited search: a
representation of length ``k`` with current remainder ``r`` and ``t`` terms still to
place uses a next denominator ``n >= max(n_prev + 1, ceil(1/r))``, and feasibility
is pruned by ``t / n >= r`` and by the requirement that the remaining reciprocals
can still sum to ``r``.  The search is exhaustive below the depth cap, so a
failure to find a representation of length ``k`` certifies ``N(a, b) > k``.
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np

__all__ = [
    "greedy_expansion",
    "shortest_expansion",
    "max_min_length",
    "max_min_length_table",
    "greedy_length_table",
    "double_log_envelope",
]


def greedy_expansion(a: int, b: int) -> list[int]:
    """Fibonacci--Sylvester greedy expansion of ``a / b`` into distinct unit fractions.

    Returns the denominators in increasing order.  The greedy expansion exists for
    every ``0 < a < b`` and its length is ``O(log b / log log b)`` in the worst case,
    so it is a valid but not minimal expansion.
    """
    if not (0 < a < b):
        raise ValueError("need 0 < a < b")
    out: list[int] = []
    num, den = a, b
    while num > 0:
        n = -(-den // num)  # ceil(den / num)
        if out and n <= out[-1]:
            n = out[-1] + 1
        out.append(n)
        num, den = num * n - den, den * n
    return out


def greedy_expansion_exact(a: int, b: int) -> list[int]:
    """Greedy expansion using exact rational arithmetic, as a cross-check."""
    if not (0 < a < b):
        raise ValueError("need 0 < a < b")
    out: list[int] = []
    rest = Fraction(a, b)
    while rest > 0:
        n = -(-rest.denominator // rest.numerator)
        if out and n <= out[-1]:
            n = out[-1] + 1
        out.append(n)
        rest -= Fraction(1, n)
    return out


def shortest_expansion(a: int, b: int, max_terms: int = 10) -> list[int] | None:
    """A shortest distinct unit-fraction expansion of ``a / b``, or None.

    Iterative deepening on the number of terms, with exhaustive search at each
    depth using exact rational arithmetic.  Two bounds prune the tree at a node
    with remainder ``r``, ``t`` terms still to place and smallest usable
    denominator ``low``:

    * ``t / low >= r``, since the remaining terms are at most ``1/low`` each;
    * ``1 / (low + t - 1) <= r``, since the remaining terms are at least
      ``1/(low + t - 1)`` each.

    A return of None certifies that no expansion of length at most ``max_terms``
    exists, so the true minimum length exceeds ``max_terms``.
    """
    if not (0 < a < b):
        raise ValueError("need 0 < a < b")
    for k in range(1, max_terms + 1):
        found = _dfs(Fraction(a, b), k, 2, [])
        if found is not None:
            return found
    return None


def _dfs(rem: Fraction, t: int, low: int, acc: list[int]) -> list[int] | None:
    """Exhaustive search for ``t`` distinct unit fractions, each at least ``low``.

    ``rem`` must end up as exactly zero.  Two bounds drive the pruning of the next
    denominator ``n``:

    * ``1 / n <= rem``, since the next term alone may not overshoot;
    * ``sum_{i=n}^{n+t-1} 1/i >= rem``, because the ``t`` smallest usable denominators
      are ``n, n+1, ..., n+t-1`` and no choice of ``t`` distinct denominators at
      least ``n`` can sum to less than that.  This sum decreases in ``n``, so it cuts
      off the tail of the range exactly.
    """
    if t == 0:
        return list(acc) if rem == 0 else None
    if rem <= 0:
        return None
    if t == 1:
        # A single term must equal 1/rem exactly.
        if rem.numerator == 1 and rem.denominator >= low:
            return acc + [rem.denominator]
        return None
    n_min = max(low, -(-rem.denominator // rem.numerator))
    # Upper bound: the largest n with harmonic(n, t) >= rem.  harmonic is strictly
    # decreasing in n for fixed t >= 2, so binary search on a monotone predicate.
    if _harmonic_sum(n_min, t) < rem:
        return None
    lo, hi = n_min, n_min
    while _harmonic_sum(hi, t) >= rem:
        hi *= 2
        if hi > 10**150:  # pragma: no cover - unreachable for t >= 2
            break
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if _harmonic_sum(mid, t) >= rem:
            lo = mid
        else:
            hi = mid
    n_max = lo
    for n in range(n_min, n_max + 1):
        rest = rem - Fraction(1, n)
        if rest < 0:
            break
        got = _dfs(rest, t - 1, n + 1, acc + [n])
        if got is not None:
            return got
    return None


def _harmonic_sum(n: int, t: int) -> Fraction:
    """``1/n + 1/(n+1) + ... + 1/(n+t-1)`` as an exact rational."""
    return sum((Fraction(1, n + i) for i in range(t)), Fraction(0))


def max_min_length(b: int, max_terms: int = 8) -> tuple[int, int]:
    """``(N(b), argmax)`` for the given denominator.

    ``N(b) = max_{1 <= a < b} N(a, b)``, and the second entry is the numerator
    attaining it.  Numerators are scanned in decreasing order because ``a = b - 1``
    is the classical extremal case, so the maximum is usually found early and the
    remaining numerators can be skipped once the cap is reached.
    """
    if b < 2:
        raise ValueError("b must be at least 2")
    best, arg = 0, 1
    for a in range(b - 1, 0, -1):
        exp = shortest_expansion(a, b, max_terms=max_terms)
        length = 0 if exp is None else len(exp)
        if exp is None:
            # No expansion within the cap; the true value is larger than max_terms.
            length = max_terms + 1
        if length > best:
            best, arg = length, a
        if best >= max_terms:
            break
    return best, arg


def max_min_length_table(b_max: int, max_terms: int = 8) -> tuple[np.ndarray, np.ndarray]:
    """``(N(b), argmax)`` for every ``2 <= b <= b_max``."""
    if b_max < 2:
        raise ValueError("b_max must be at least 2")
    bs = np.arange(2, b_max + 1)
    vals = np.empty(bs.size, dtype=np.int64)
    args = np.empty(bs.size, dtype=np.int64)
    for i, b in enumerate(bs):
        vals[i], args[i] = max_min_length(int(b), max_terms=max_terms)
    return bs, np.stack([vals, args])


def greedy_length_table(b_max: int, step: int = 1) -> tuple[np.ndarray, np.ndarray]:
    """Worst-case greedy expansion length over ``a`` for each ``b <= b_max``."""
    if b_max < 2:
        raise ValueError("b_max must be at least 2")
    bs = np.arange(2, b_max + 1, step)
    out = np.empty(bs.size, dtype=np.int64)
    for i, b in enumerate(bs):
        out[i] = max(len(greedy_expansion(a, int(b))) for a in range(1, int(b)))
    return bs, out


def double_log_envelope(b: np.ndarray | float, c: float = 1.0) -> np.ndarray:
    """The conjectured order ``c log log b``."""
    b_arr = np.atleast_1d(np.asarray(b, dtype=float))
    if np.any(b_arr <= 1):
        raise ValueError("log log b needs b > 1")
    return c * np.log(np.log(b_arr))


def erdos_graham_envelope(b: np.ndarray | float) -> np.ndarray:
    """Erdős's earlier bound shape ``log b / log log b``."""
    b_arr = np.atleast_1d(np.asarray(b, dtype=float))
    return np.log(b_arr) / np.log(np.log(b_arr))


def vose_envelope(b: np.ndarray | float) -> np.ndarray:
    """Vose's bound shape ``sqrt(log b)``."""
    b_arr = np.atleast_1d(np.asarray(b, dtype=float))
    return np.sqrt(np.log(b_arr))
