"""The irrationality exponent of ``pi`` (family 017).

The irrationality exponent of an irrational real ``x`` is

    mu(x) = sup { nu > 0 : 0 < |x - p/q| < q^{-nu} for infinitely many p, q },

with ``p/q`` in lowest terms and ``q >= 2``.  Pigeonhole gives ``mu(x) >= 2``, and
the classical convergent theory of continued fractions identifies the extremal
approximations: ``mu(x) = 2`` exactly when the partial quotients are unbounded and
the convergents satisfy ``q |q x - p| -> 0``.

The manuscript *The irrationality exponent of ``pi`` is 2* proves

    mu(pi) = 2,   and for every nu > 2 there is Q(nu) with
    |pi - p/q| >= q^{-nu}  for all p in Z, q >= Q(nu),

using all integer numerators and denominators, reduced or not.  The nonexistence of
a uniform positive lower bound near ``q = 1`` is why the threshold ``Q(nu)`` is
needed.  As a consequence the Flint--Hills series ``sum 1/(n^3 sin^2 n)`` converges.

The computation here is a numerical study of the convergents rather than a proof:
it reports the continued fraction, the quantity ``q^2 |pi - p/q|`` for each
convergent, and the empirical exponent ``-log|pi - p/q| / log q``.
"""

from __future__ import annotations

import mpmath as mp
import numpy as np

__all__ = [
    "continued_fraction",
    "convergents",
    "convergents_exact",
    "to_float",
    "approximation_table",
    "empirical_exponent",
    "semiconvergents",
    "best_approximation_quality",
]


def continued_fraction(x: mp.mpf, terms: int, dps: int = 80) -> list[int]:
    """Partial quotients ``a_0, a_1, ...`` of ``x`` to ``terms`` terms.

    Computed by the standard recurrence on the fractional part, at ``dps`` decimal
    digits of working precision.
    """
    if terms < 1:
        raise ValueError("terms must be at least 1")
    with mp.workdps(dps):
        v = mp.mpf(x)
        out: list[int] = []
        for _ in range(terms):
            a = int(mp.floor(v))
            out.append(a)
            frac = v - a
            if frac == 0:
                break
            v = 1 / frac
        return out


def convergents(cf: list[int]) -> tuple[np.ndarray, np.ndarray]:
    """Numerators and denominators of the convergents of ``cf``.

    ``p_k / q_k`` with ``p_{-2}=0, p_{-1}=1`` and ``q_{-2}=1, q_{-1}=0``.  The
    arrays have the same length as ``cf`` and are float64, which limits the study to
    convergents whose denominator is below ``2^53``; that is plenty, since the
    interesting regime is the approach of ``q^2 |x - p/q|`` to a small value.
    """
    if not cf:
        raise ValueError("cf must be non-empty")
    p_m2, p_m1 = 0, 1
    q_m2, q_m1 = 1, 0
    ps = np.empty(len(cf), dtype=float)
    qs = np.empty(len(cf), dtype=float)
    for i, a in enumerate(cf):
        p = a * p_m1 + p_m2
        q = a * q_m1 + q_m2
        ps[i], qs[i] = p, q
        p_m2, p_m1 = p_m1, p
        q_m2, q_m1 = q_m1, q
    return ps, qs


def convergents_exact(cf: list[int]) -> tuple[list[int], list[int]]:
    """Numerators and denominators of the convergents of ``cf``, as exact integers.

    ``p_{-2} = 0``, ``p_{-1} = 1``, ``q_{-2} = 1``, ``q_{-1} = 0``, and
    ``p_k = a_k p_{k-1} + p_{k-2}``, ``q_k = a_k q_{k-1} + q_{k-2}``.  Python integers
    are unbounded, so no precision is lost however large the denominators become.
    """
    if not cf:
        raise ValueError("cf must be non-empty")
    p_m2, p_m1 = 0, 1
    q_m2, q_m1 = 1, 0
    ps: list[int] = []
    qs: list[int] = []
    for a in cf:
        p = a * p_m1 + p_m2
        q = a * q_m1 + q_m2
        ps.append(p)
        qs.append(q)
        p_m2, p_m1 = p_m1, p
        q_m2, q_m1 = q_m1, q
    return ps, qs


def approximation_table(
    x: mp.mpf, terms: int = 40, dps: int = 120
) -> dict[str, object]:
    """Continued fraction data for ``x``, with exact integer convergents.

    Returns a dict with

    * ``a``: the partial quotients;
    * ``p``, ``q``: lists of Python integers, the numerators and denominators of the
      convergents, exact and unbounded in size;
    * ``err``: ``|x - p/q|`` as arbitrary-precision ``mpf`` values; and
    * ``quality``, ``squared_quality``: ``q*err`` and ``q**2*err`` as ``mpf``.

    The point of keeping integers and ``mpf`` values rather than converting to
    float64 is that the convergents of ``pi`` reach denominators above ``10**20``
    within forty terms, where float64 no longer represents ``q`` exactly and the
    squared quality ``q**2 |x - p/q|`` is destroyed by rounding.  Callers that want
    arrays for plotting should convert with :func:`to_float`.
    """
    cf = continued_fraction(x, terms, dps=dps)
    p_int, q_int = convergents_exact(cf)
    ps = [int(v) for v in p_int]
    qs = [int(v) for v in q_int]
    with mp.workdps(dps):
        errs = [abs(mp.mpf(x) - mp.mpf(p) / mp.mpf(q)) for p, q in zip(ps, qs)]
        quality = [mp.mpf(q) * e for q, e in zip(qs, errs)]
        squared = [mp.mpf(q) ** 2 * e for q, e in zip(qs, errs)]
    return {
        "a": np.array(cf, dtype=np.int64),
        "p": ps,
        "q": qs,
        "err": errs,
        "quality": quality,
        "squared_quality": squared,
    }


def to_float(values) -> np.ndarray:
    """Convert a list of Python ints or ``mpf`` values to a float64 array."""
    return np.array([float(v) for v in values], dtype=float)


def empirical_exponent(errors, denominators) -> np.ndarray:
    """The empirical exponent ``-log|err| / log q``, a lower estimate for ``mu``.

    For a convergent of an irrational with bounded partial quotients this tends to
    the value ``1 + liminf log q_{k+1} / log q_k``, which for bounded partial
    quotients equals 2.  Being an infimum-like statistic over finitely many
    approximants, it is a lower bound for ``mu``, never an upper bound.
    """
    err = np.asarray([float(e) for e in errors], dtype=float)
    den = np.asarray([float(d) for d in denominators], dtype=float)
    if err.shape != den.shape:
        raise ValueError("errors and denominators must have the same shape")
    if np.any(err <= 0):
        raise ValueError("errors must be positive")
    if np.any(den <= 1):
        raise ValueError("denominators must exceed 1")
    return -np.log(err) / np.log(den)


def semiconvergents(cf: list[int], index: int) -> list[tuple[int, int]]:
    """Intermediants between the convergents of ``index`` and ``index + 1``.

    Returns ``(t p_index + p_{index-1}) / (t q_index + q_{index-1})`` for
    ``1 <= t < a_{index+1}``, which are the best approximations of the second kind
    that are not themselves convergents.  These matter for the irrationality
    exponent when a partial quotient is unusually large, since they can approach
    ``x`` much more closely than the size of the partial quotient alone suggests.
    The recurrence seeds are ``p_{-1} = 0`` and ``q_{-1} = 1``; taking ``t = 1``
    returns the convergent at ``index`` itself.
    """
    if index < 0 or index + 1 >= len(cf):
        raise ValueError("index out of range for semiconvergents")
    ps, qs = convergents_exact(cf)
    p_prev = ps[index - 1] if index >= 1 else 0
    q_prev = qs[index - 1] if index >= 1 else 0
    return [
        (t * ps[index] + p_prev, t * qs[index] + q_prev)
        for t in range(1, int(cf[index + 1]))
    ]


def best_approximation_quality(
    x: mp.mpf, q_max: int, dps: int = 120
) -> tuple[np.ndarray, np.ndarray]:
    """For each convergent with ``q <= q_max``, the squared quality ``q^2 |x - p/q|``.

    Returns ``(q, squared_quality)`` restricted to ``q <= q_max``.
    """
    table = approximation_table(x, terms=60, dps=dps)
    keep = [i for i, q in enumerate(table["q"]) if 2 <= q <= q_max]
    return to_float([table["q"][i] for i in keep]), to_float(
        [table["squared_quality"][i] for i in keep]
    )
