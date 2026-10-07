"""The Poisson--Dirichlet law PD(1) and the quantities of the prime-predecessor law.

The manuscript *The Poisson--Dirichlet Law for Prime Predecessors* proves that for
a uniformly chosen prime ``p <= x``, the numbers

    V_j(p) = log q_j(p) / log(p - 1),

where ``q_1 >= q_2 >= ...`` are the prime factors of ``p - 1`` counted with
multiplicity, converge jointly, for every fixed ``k``, to the first ``k``
components of the decreasing rearrangement of the stick-breaking sequence

    B_1 = 1 - U_1,   B_j = (U_1 ... U_{j-1}) (1 - U_j),

whose decreasing rearrangement is PD(1).

This module supplies

* the decreasing-rearrangement order statistics of PD(1), both by simulation and
  by the exact recursion (3) below,
* the Dickman function, which is the law of the largest component ``L_1``,
* the predecessor factorisation itself, via a largest-prime-factor sieve.

The exact recursion.  PD(1) satisfies the self-similarity ``L^(1) =_d (1 - U) +
U L^(2)`` with ``L^(1)``, ``L^(2)`` independent and PD(1).  Writing ``y = 1 - u``
for the first fragment and ``L_1 = 1 - U_1``,

    Pr[L_j <= x] = int_0^1 Pr[L_1 <= x | U_1 = u] * Pr[L_{j-1} <= (x - 1 + u)/u] du.

The condition ``L_1 <= x`` is ``u >= 1 - x`` and the residual argument is below 1
exactly when ``u <= x``.  Substituting ``u = 1 - (1 - x) t`` removes the endpoint
singularity:

    F_j(x) = (1 - x) int_0^1 F_{j-1}(1 - (1 - x)(1 - t)/ (1 - (1 - x) t)) dt.

For ``j = 1`` the inner condition is automatic and ``F_1(x) = x``, so
``F_1`` is the Dickman law ``rho(1/x)``.
"""

from __future__ import annotations

import numpy as np

from .sieves import factorize_from_lpf, lpf_sieve, primes_upto

__all__ = [
    "dickman",
    "stick_breaking_pd1",
    "pd1_order_statistic_cdf",
    "predecessor_factors",
    "normalized_log_factors",
]

#: Largest ``u`` at which :func:`dickman` is tabulated.  The probes need ``u <= 5``
#: (the largest is ``1/a`` for the marginal at ``a = 0.2``).  The value is set high
#: enough that :func:`openmath.joint_dickman.dickman_marginal` and
#: :func:`pd1_order_statistic_cdf` can invert it for small arguments; by ``u = 12``
#: the function is around ``1e-11``, far below any quantity plotted.
MAX_DICKMAN_U = 12.0


#: Number of stick fragments formed per PD(1) draw.  The discarded tail is
#: ``prod_{j <= top} U_j``, whose *mean* is ``2^-DEFAULT_PD1_TOP`` (about ``2.4e-4``
#: here) and whose median is ``(log 2)^top`` (about ``1.3e-4`` at this value).
DEFAULT_PD1_TOP = 12

#: Largest ``u`` at which :func:`dickman` is tabulated per call.
_MAX_ODE_STEPS = 4_000_000


def dickman(u: np.ndarray | float, step: float = 1e-4) -> np.ndarray:
    """The Dickman--de Bruijn function ``rho``, by fourth-order integration.

    ``rho`` is determined by

        rho(u) = 1            (0 <= u <= 1),
        u rho'(u) = -rho(u-1) (u > 1).

    On ``[1, 2]`` the solution is the elementary ``rho(u) = 1 - log u``; beyond
    that the delay equation is integrated from ``u = 2`` by classical
    fourth-order Runge--Kutta on the uniform step ``step``, with the lagged value
    ``rho(u - 1)`` taken by linear interpolation of the already-computed grid.
    Interpolation of the lag is the only source of error, and it is ``O(step^2)``
    with a small constant, so the default step reproduces the published values to
    about twelve digits:

        rho(2)   = 0.3068528194400547,
        rho(2.5) = 0.1303195618322514,
        rho(3)   = 0.04860838829113157.

    The grid is cached by ``step``, so repeated calls with the same step reuse one
    integration.
    """
    scalar = np.ndim(u) == 0
    u_arr = np.atleast_1d(np.asarray(u, dtype=float))
    if np.any(u_arr < 0):
        raise ValueError("rho is defined here for u >= 0")
    u_max = float(u_arr.max())
    if u_max > MAX_DICKMAN_U:
        raise ValueError(
            f"rho is tabulated here only up to u = {MAX_DICKMAN_U}; got {u_max}"
        )
    if u_max <= 2.0:
        out = np.where(u_arr <= 1.0, 1.0, 1.0 - np.log(np.maximum(u_arr, 1.0)))
        return out[0] if scalar else out

    grid, table = _dickman_grid(u_max, step)
    out = np.empty_like(u_arr)
    low = u_arr <= 1.0
    out[low] = 1.0
    mid = (u_arr > 1.0) & (u_arr <= 2.0)
    out[mid] = 1.0 - np.log(u_arr[mid])
    high = u_arr > 2.0
    out[high] = np.interp(u_arr[high], grid, table)
    return out[0] if scalar else out


_DICKMAN_CACHE: dict[float, tuple[np.ndarray, np.ndarray]] = {}


def _dickman_grid(u_max: float, step: float) -> tuple[np.ndarray, np.ndarray]:
    """Integrate the delay equation on ``[2, u_max + 1]`` and cache the result."""
    n = int(np.ceil(u_max / step)) + 2
    key = step
    cached = _DICKMAN_CACHE.get(key)
    if cached is not None and cached[0][-1] >= u_max:
        # Return the cached grid unchanged.  Trimming it would make the lag
        # interpolation read past the end of the array, since the recursion needs
        # rho at points up to one unit below the largest requested u.
        return cached
    if n > _MAX_ODE_STEPS:
        raise ValueError("requested range needs too many ODE steps; raise step")

    grid = np.arange(n + 1) * step
    table = np.empty(n + 1, dtype=float)
    table[grid <= 1.0] = 1.0
    band = (grid > 1.0) & (grid <= 2.0)
    table[band] = 1.0 - np.log(grid[band])

    # First grid node at or above u = 2: its value is closed form, so the
    # recursion starts there and every lagged query below it is already filled.
    i = int(np.searchsorted(grid, 2.0, side="left"))
    if not np.isclose(grid[i], 2.0):
        raise ValueError("the step must divide 1.0 so that u = 2 is a grid node")
    rho = lambda x: np.interp(x, grid[: i + 1], table[: i + 1])  # noqa: E731
    while i < n:
        t, r = grid[i], table[i]
        k1 = -rho(t - 1.0) / t
        k2 = -rho(t + 0.5 * step - 1.0) / (t + 0.5 * step)
        k3 = -rho(t + 0.5 * step - 1.0) / (t + 0.5 * step)
        k4 = -rho(t + step - 1.0) / (t + step)
        table[i + 1] = r + step / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        i += 1
    _DICKMAN_CACHE[key] = (grid, table)
    return grid, table


def stick_breaking_pd1(
    n_samples: int,
    top: int = 12,
    seed: int | None = 0,
    batch: int = 400_000,
) -> np.ndarray:
    """Draw the ``top`` largest components of PD(1), in decreasing order.

    PD(1) is the decreasing rearrangement of the stick-breaking sequence
    ``B_1 = 1 - U_1`` and ``B_j = (U_1 ... U_{j-1})(1 - U_j)``.  A finite
    ``Dirichlet(1, ..., 1)`` does *not* approximate it -- the largest component of
    a length-``N`` Dirichlet is only of order ``log N / N`` -- so the stick
    breaking is carried out directly, with ``log B_j`` accumulated to avoid
    underflow.  Only the first ``top`` fragments are formed, which drops a tail of
    mass ``prod_{i <= top} U_i``; its median is ``2^-top``, i.e. about ``2.4e-4``
    at the default, and its mean is ``3^-top``, about ``1.9e-6``.  The tail is not
    renormalised away: the returned entries are genuine fragments of a partition
    of unity, so their sum is slightly below one by the discarded mass.

    Draws are produced in batches, so the memory cost is ``batch * top``.
    """
    if n_samples < 1:
        raise ValueError("n_samples must be positive")
    if top < 1:
        raise ValueError("top must be at least 1")
    rng = np.random.default_rng(seed)
    out = np.empty((n_samples, top), dtype=float)
    done = 0
    while done < n_samples:
        m = min(batch, n_samples - done)
        u = rng.random((m, top))
        log_u = np.log(u)
        # log B_j = sum_{i<j} log U_i + log(1 - U_j)
        log_b = np.cumsum(log_u, axis=1) - log_u + np.log1p(-u)
        out[done : done + m] = np.exp(np.sort(log_b, axis=1)[:, ::-1])
        done += m
    return out


def pd1_order_statistic_cdf(
    x: np.ndarray | float,
    j: int,
    samples: int = 2_000_000,
    seed: int = 0,
) -> np.ndarray:
    """CDF of the ``j``-th largest PD(1) component, ``Pr[L_j <= x]``.

    For ``j = 1`` the law is the Dickman distribution, ``Pr[L_1 <= x] = rho(1/x)``,
    and that closed form is used.  For ``j >= 2`` the CDF is estimated from
    ``samples`` draws of :func:`stick_breaking_pd1`; the Monte Carlo standard
    error is at most ``0.5 / sqrt(samples)``, i.e. ``3.5e-4`` at the default.
    """
    if j < 1:
        raise ValueError("j must be at least 1")
    scalar = np.ndim(x) == 0
    x_arr = np.atleast_1d(np.asarray(x, dtype=float))
    if np.any((x_arr < 0.0) | (x_arr > 1.0)):
        raise ValueError("the PD(1) components live in [0, 1]")
    if j == 1:
        cdf = np.where(x_arr <= 0.0, 0.0, dickman(1.0 / np.maximum(x_arr, 1e-300)))
        return cdf[0] if scalar else cdf
    draws = stick_breaking_pd1(samples, seed=seed)
    col = np.sort(draws[:, j - 1])
    cdf = np.searchsorted(col, x_arr, side="right") / samples
    return cdf[0] if scalar else cdf


def pd1_order_statistic_mean(
    j: int, samples: int = 2_000_000, seed: int = 1
) -> float:
    """Mean of the ``j``-th largest PD(1) component, by simulation.

    The standard error is at most ``0.5 / sqrt(samples)``, about ``3.5e-4`` at the
    default.  The means decay faster than geometrically: ``E[L_j]`` is
    approximately ``0.624 / (j - 1)!`` for ``j >= 2``, and ``sum_j E[L_j] = 1``
    exactly, since the components are a partition of unity.  The tests check that
    the simulated means sum to one within Monte Carlo error.
    """
    if j < 1:
        raise ValueError("j must be at least 1")
    draws = stick_breaking_pd1(max(samples, 1), top=DEFAULT_PD1_TOP, seed=seed)
    return float(draws[:, j - 1].mean())



def predecessor_factors(limit: int) -> tuple[np.ndarray, list[list[int]]]:
    """Factor ``p - 1`` for every prime ``3 <= p <= limit``.

    Returns the array of primes and, for each, the prime factors of ``p - 1`` in
    non-increasing order with multiplicity.
    """
    if limit < 5:
        raise ValueError("limit must be at least 5")
    lpf = lpf_sieve(limit)
    ps = primes_upto(limit)
    ps = ps[ps >= 3]
    factors = [factorize_from_lpf(int(p) - 1, lpf) for p in ps]
    return ps, factors


def normalized_log_factors(factors: list[int], p: int, k: int) -> np.ndarray:
    """The vector ``(V_1(p), ..., V_k(p))`` padded with zeros after the list ends.

    The manuscript sets ``q_j(p) = 1`` once the list of factors is exhausted, so
    ``V_j(p) = 0`` there and ``sum_j V_j(p) = 1`` exactly.
    """
    v = np.zeros(k, dtype=float)
    denom = np.log(p - 1.0)
    for idx in range(min(k, len(factors))):
        v[idx] = np.log(factors[idx]) / denom
    return v
