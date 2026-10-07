"""The joint Dickman law for consecutive integers (family 012).

The manuscript *The joint Dickman law for consecutive integers* proves that for
every fixed ``a, b`` in ``(0, 1)``,

    (1/X) #{2 <= n <= X : P+(n) <= n^a, P+(n+1) <= n^b} -> rho(1/a) rho(1/b),

where ``P+(n)`` is the largest prime factor of ``n`` and ``rho`` is the
Dickman--de Bruijn function.  Equivalently, the two coordinates

    x_n = log P+(n) / log n,     y_n = log P+(n+1) / log n

are asymptotically independent with common distribution function ``a -> rho(1/a)``.
Two consequences are probed here:

* the comparison ``P+(n) < P+(n+1)`` has natural density ``1/2``;
* the joint law factorises, so the joint CDF is the product of the marginals.

Both ``x_n`` and ``y_n`` are restricted to ``[0, 1]`` together with
``P+(1) = 1``, and ``x_n = 0`` occurs only for ``n = 1``.
"""

from __future__ import annotations

import numpy as np

from .poisson_dirichlet import MAX_DICKMAN_U, dickman

__all__ = [
    "largest_prime_factors",
    "normalized_coordinates",
    "comparison_density",
    "joint_cdf_table",
    "independence_gap",
]


def largest_prime_factors(limit: int) -> np.ndarray:
    """``P+(n)`` for ``0 <= n <= limit``, with ``P+(0) = 0`` and ``P+(1) = 1``.

    Computed by sieving the primes at most ``limit`` in increasing order and
    stamping each prime into every multiple by strided assignment; the last stamp
    is the largest prime factor.
    """
    if limit < 2:
        raise ValueError("limit must be at least 2")
    from .sieves import primes_upto

    out = np.ones(limit + 1, dtype=np.int64)
    out[0] = 0
    for p in primes_upto(limit):
        p = int(p)
        out[p::p] = p
    return out


def normalized_coordinates(
    lpf: np.ndarray, start: int = 2
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(x_n, y_n)`` for ``start <= n <= len(lpf) - 2``.

    Both coordinates are ``log P+(m) / log m``, so each lies in ``[0, 1]`` because
    ``P+(m) <= m`` and ``log m > 0`` for ``m >= 2``.
    """
    n = np.arange(start, len(lpf) - 1, dtype=np.int64)
    logn = np.log(n.astype(float))
    x = np.log(lpf[n].astype(float)) / logn
    y = np.log(lpf[n + 1].astype(float)) / logn
    return x, y


def comparison_density(x: np.ndarray, y: np.ndarray) -> float:
    """Empirical density of ``P+(n) < P+(n+1)``, i.e. of ``x_n < y_n``."""
    if x.shape != y.shape:
        raise ValueError("x and y must have the same shape")
    if x.size == 0:
        raise ValueError("need at least one sample")
    return float(np.mean(x < y))


def joint_cdf_table(
    x: np.ndarray, y: np.ndarray, grid: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Binned joint and product-of-marginals CDFs on ``grid x grid``.

    Returns ``(joint, product, counts)``.  ``joint[i, j]`` estimates
    ``Pr[x <= grid[i], y <= grid[j]]``, ``product[i, j]`` is the product of the two
    one-dimensional empirical CDFs, and ``counts[i, j]`` is the number of samples in
    the cell.  The grid must be strictly increasing.

    The bins tile ``[grid[0], inf) x [grid[0], inf)``, so ``counts.sum()`` is the
    number of samples with both coordinates at least ``grid[0]``.  Samples in
    ``[0, grid[0])`` are excluded, which is why a caller that needs the CDF to reach
    one should start the grid at zero.
    """
    g = np.asarray(grid, dtype=float)
    if g.ndim != 1 or g.size < 1 or np.any(np.diff(g) <= 0):
        raise ValueError("grid must be a strictly increasing 1-D array")
    if x.shape != y.shape:
        raise ValueError("x and y must have the same shape")
    # The last bin must catch everything above grid[-1].  A finite edge is used
    # rather than inf: passing both an infinite edge and a range makes numpy clip
    # the edges to the data maximum and drop every sample above grid[-1].
    hi = float(max(x.max(), y.max())) + 1.0
    edges = np.append(g, hi)
    hist, _, _ = np.histogram2d(x, y, bins=[edges, edges])
    joint = np.cumsum(np.cumsum(hist, axis=0), axis=1)[: g.size, : g.size]
    marginal_x = np.cumsum(hist.sum(axis=1))[: g.size]
    marginal_y = np.cumsum(hist.sum(axis=0))[: g.size]
    total = x.size
    joint = joint / total
    product = np.outer(marginal_x, marginal_y) / (total * total)
    return joint, product, hist


def independence_gap(joint: np.ndarray, product: np.ndarray) -> float:
    """Supremum distance between the joint CDF and the product of marginals."""
    return float(np.max(np.abs(joint - product)))


def dickman_marginal(a: np.ndarray | float) -> np.ndarray | float:
    """The limiting marginal ``a -> rho(1/a)``, with value 1 at ``a = 0``.

    The manuscript continuously extends the law to ``[0, 1]``; ``rho(t) -> 0``
    super-exponentially as ``t -> inf``, so ``a -> 0`` is filled in with the limit
    zero beyond the tabulated range.
    """
    scalar = np.ndim(a) == 0
    a_arr = np.atleast_1d(np.asarray(a, dtype=float))
    if np.any((a_arr < 0) | (a_arr > 1)):
        raise ValueError("a must lie in [0, 1]")
    out = np.ones_like(a_arr)
    pos = a_arr > 0
    if np.any(pos):
        inv = 1.0 / a_arr[pos]
        cap = MAX_DICKMAN_U
        out[pos] = np.where(inv <= cap, np.atleast_1d(dickman(np.minimum(inv, cap))), 0.0)
    return float(out[0]) if scalar else out
