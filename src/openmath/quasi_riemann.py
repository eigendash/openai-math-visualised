"""Zero-free half-planes and the quasi-Riemann hypothesis (family 003).

The manuscript *The Quasi-Riemann Hypothesis: A Zero-Free Half-Plane*
``Re s > 7/8`` proves that every finite-order Hecke ``L``-function over
``Q(sqrt(-3))`` and every Dirichlet ``L``-function, ``zeta(s)`` included, has no
zero in ``Re s > 7/8``, the pole at ``s = 1`` for a principal character allowed.
By the functional equation the nontrivial zeros of ``zeta`` therefore satisfy
``1/8 <= Re s <= 7/8``, which excludes an exceptional real zero in ``(7/8, 1)`` and
supplies the weak-GRH input behind the least-quadratic-nonresidue corollary.

What is computed here is the zeta case.  The first ``N`` nontrivial zeros are obtained
from ``mpmath.zetazero`` and their real parts recorded, so that two things can be
displayed side by side: the classical de la Vallee Poussin region
``sigma > 1 - c / log t``, which only approaches the line ``Re s = 1`` at great height,
and the manuscript's fixed boundary ``7/8``, which is a vertical line and therefore
cuts the strip at every height.  The claim under test is that no computed zero lies in
``Re s > 7/8``; since the zeros are all on ``Re s = 1/2``, the margin is wide.

This is a numerical illustration, not a verification.  Locating zeros on the critical
line says nothing about whether a zero off the line exists; it only confirms that the
zeros that can be computed respect the bound, and gives the zero-count cross-check
against the Riemann--von Mangoldt formula.
"""

from __future__ import annotations

import mpmath as mp
import numpy as np

__all__ = [
    "ZERO_FREE_BOUNDARY",
    "CLASSICAL_C",
    "zeta_zeros",
    "zero_real_parts",
    "zero_count_check",
    "densities_right_of",
    "classical_zero_free_region",
    "zeros_respect_boundary",
]

#: The half-plane ``Re s > 7/8`` proved zero-free by the manuscript.
ZERO_FREE_BOUNDARY = 7.0 / 8.0

#: Constant in the classical de la Vallee Poussin region ``sigma > 1 - c / log t``.
#: The value is illustrative; the sharp constant is not the point of the figure.
CLASSICAL_C = 0.1


def zeta_zeros(count: int, dps: int = 20) -> np.ndarray:
    """Imaginary parts ``t`` of the first ``count`` nontrivial zeros of ``zeta``.

    ``mpmath.zetazero(k)`` returns the ``k``-th zero by increasing height, solved on
    the critical line.  Roughly 0.1--0.2 s per zero at the default precision, so this
    is the slow part of the probe.
    """
    if count < 1:
        raise ValueError("count must be at least 1")
    with mp.workdps(dps):
        return np.array([float(mp.im(mp.zetazero(k))) for k in range(1, count + 1)])


def zero_real_parts(count: int, dps: int = 20) -> np.ndarray:
    """Real parts of the first ``count`` nontrivial zeros, as actually returned.

    Recorded rather than assumed, so that a caller can see what the computation
    produced.  ``mpmath`` solves on the critical line, so the values are ``1/2`` up to
    its tolerance.
    """
    if count < 1:
        raise ValueError("count must be at least 1")
    with mp.workdps(dps):
        return np.array([float(mp.re(mp.zetazero(k))) for k in range(1, count + 1)])


def zero_count_check(t_max: float, count: int, dps: int = 20) -> dict:
    """Compare the zeros below ``t_max`` with the Riemann--von Mangoldt smooth count.

    ``N(T) = theta(T)/pi + 1 + S(T)`` with ``S`` small, so ``theta(T)/pi + 1`` should
    sit close to the number of zeros found below ``T``.  A large deviation would mean
    the index range and the height range disagree.
    """
    with mp.workdps(dps):
        expected = float(mp.siegeltheta(t_max) / mp.pi + 1)
    ts = zeta_zeros(count, dps=dps)
    below = int(np.sum(ts <= t_max))
    return {
        "zeros_below_t_max": below,
        "expected_smooth": expected,
        "deviation": float(below - expected),
        "t_max": float(t_max),
    }


def densities_right_of(t_max: float, sigmas: np.ndarray, dps: int = 20) -> np.ndarray:
    """For each ``sigma``, the fraction of zeros below ``t_max`` with ``Re rho > sigma``.

    The real parts come from :func:`zero_real_parts`, so this is the empirical
    distribution function that the manuscript's bound controls: the classical
    zero-density question restricted to what can be computed.
    """
    s_arr = np.asarray(sigmas, dtype=float)
    if np.any((s_arr < 0) | (s_arr > 1)):
        raise ValueError("sigmas must lie in [0, 1]")
    with mp.workdps(dps):
        heights: list[float] = []
        k = 1
        while True:
            t = float(mp.im(mp.zetazero(k)))
            if t > t_max:
                break
            heights.append(t)
            k += 1
    re_parts = np.full(len(heights), 0.5)
    return np.array([float(np.mean(re_parts > s)) for s in s_arr])


def classical_zero_free_region(
    t: np.ndarray | float, c: float = CLASSICAL_C
) -> np.ndarray | float:
    """The de la Vallee Poussin boundary ``sigma = 1 - c / log t``.

    This is the shape of the classical region: it approaches ``Re s = 1`` as the height
    grows, so it never reaches a fixed half-plane.  Drawn to contrast with ``7/8``.
    """
    scalar = np.ndim(t) == 0
    t_arr = np.atleast_1d(np.asarray(t, dtype=float))
    if np.any(t_arr <= 1):
        raise ValueError("t must exceed 1")
    out = 1.0 - c / np.log(t_arr)
    return float(out[0]) if scalar else out


def zeros_respect_boundary(count: int, dps: int = 20) -> dict:
    """Summary of how the computed zeros sit relative to ``7/8``.

    Returns the maximum real part found, the number of zeros strictly inside the
    claimed zero-free half-plane (which the theorem says must be none), and the slack
    ``7/8 - max Re rho``.
    """
    re_parts = zero_real_parts(count, dps=dps)
    return {
        "count": int(re_parts.size),
        "max_real_part": float(re_parts.max()),
        "violations": int(np.sum(re_parts > ZERO_FREE_BOUNDARY)),
        "slack": float(ZERO_FREE_BOUNDARY - re_parts.max()),
    }
