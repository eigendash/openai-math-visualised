#!/usr/bin/env python
"""Compute the data series consumed by the animated GitHub Pages site.

Writes ``docs/data.js``, which assigns one object to ``window.OM_DATA``.  It is a
``.js`` file rather than ``.json`` so that the page works when opened from the
filesystem as well as from a server: a ``fetch`` of a local JSON file is blocked by
the browser's origin rules, while a ``<script>`` is not.

Every series carries the raw measurement together with the value the corresponding
manuscript predicts, so the animation can show the two against each other rather
than merely displaying numbers.  Nothing here is a proof; see the repository README.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from openmath import config as cfgmod  # noqa: E402

DOCS = ROOT / "docs"
DOCS.mkdir(exist_ok=True)


def r(x, n=5):
    """Round a float for compact embedding; JSON has no float32 distinction."""
    if isinstance(x, (list, tuple, np.ndarray)):
        return [r(v, n) for v in x]
    v = float(x)
    return round(v, n) if np.isfinite(v) else None


# --------------------------------------------------------------------------- 011
def pd_series(cfg) -> dict:
    """Prime predecessor: the empirical law against PD(1), as primes accumulate."""
    from openmath.poisson_dirichlet import (
        dickman,
        normalized_log_factors,
        predecessor_factors,
        stick_breaking_pd1,
    )

    k = 6
    ps, facs = predecessor_factors(cfg.pd_prime_limit)
    V = np.array([normalized_log_factors(f, int(p), k) for f, p in zip(facs, ps)])
    print(f"  011: {ps.size:,} primes factored")

    # The primes are already in increasing order, so prefix windows are nested
    # samples of the same law, which is exactly what "as x grows" means.
    n_windows = 40
    edges = np.linspace(0, ps.size, n_windows + 1).astype(int)
    means = [[float(V[: e, j].mean()) for j in range(k)] for e in edges[1:]]
    xs = [float(ps[e - 1]) for e in edges[1:]]

    # Density of V_1 on a fixed grid, cumulative over the same windows.
    grid = np.linspace(0.0, 1.0, 51)
    hist = np.zeros(grid.size - 1)
    hist_series = []
    for e in edges[1:]:
        h, _ = np.histogram(V[:e, 0], bins=grid, density=True)
        hist_series.append(r(h))

    sim = stick_breaking_pd1(400_000, top=12, seed=cfg.seeds["pd"])
    sim_mean = [float(sim[:, j].mean()) for j in range(k)]

    # The exact law of L_1 is Pr[L_1 <= x] = rho(1/x).  Differentiating that gives
    # the density; the cap keeps the inverse in the tabulated range, where the
    # survival function is already indistinguishable from zero.
    from openmath.poisson_dirichlet import MAX_DICKMAN_U

    # Pr[L_1 <= x] = rho(1/x), so the density is d/dx (1 - rho(1/x)).  Taking the
    # derivative numerically fails, because rho(1/x) is within 1e-9 of 1 across most
    # of the range and the gradient is lost to cancellation.  Instead use the delay
    # equation u rho'(u) = -rho(u-1) with u = 1/x:
    #     f(x) = rho(u-1) / (u x^2) = x rho(1/x - 1).
    centres = 0.5 * (grid[:-1] + grid[1:])
    xpos = np.maximum(centres, 1e-9)
    u = 1.0 / xpos
    lag = np.maximum(u - 1.0, 0.0)
    # For u <= 2 the lag is at most 1, where rho is exactly 1, giving f(x) = 1.
    lagged = np.where(lag <= 1.0, 1.0, 0.0)
    inside = (lag > 1.0) & (lag <= MAX_DICKMAN_U)
    if np.any(inside):
        lagged[inside] = np.atleast_1d(dickman(lag[inside]))
    density = lagged / u / xpos**2

    assert len(xs) == len(means), "the x axis and the means must have one entry each"
    return {
        "x": r(xs, 1),
        "means": r(means),
        "simMean": r(sim_mean),
        "grid": r(centres),
        "histSeries": hist_series,
        "theoryDensity": r(np.clip(density, 0, None)),
        "nPrimes": int(ps.size),
        "limit": cfg.pd_prime_limit,
    }


# --------------------------------------------------------------------------- 012
def dickman_series(cfg) -> dict:
    """Consecutive integers: independence and the density-1/2 comparison."""
    from openmath.joint_dickman import (
        comparison_density,
        joint_cdf_table,
        largest_prime_factors,
        normalized_coordinates,
    )

    lpf = largest_prime_factors(cfg.dickman_limit)
    x, y = normalized_coordinates(lpf)
    print(f"  012: {x.size:,} pairs")

    n_windows = 40
    edges = np.linspace(0, x.size, n_windows + 1).astype(int)
    dens = [comparison_density(x[: e], y[: e]) for e in edges[1:]]
    xs = [float(e) for e in edges[1:]]

    # The independence gap, measured cumulatively, against the binning bias.
    grid = np.linspace(0.0, 0.8, 17)
    gaps = []
    for e in edges[1:]:
        joint, product, _ = joint_cdf_table(x[:e], y[:e], grid)
        gaps.append(float(np.max(np.abs(joint - product))))

    # Marginal of x_n against the predicted rho(1/a).
    qs = np.linspace(0.05, 0.95, 19)
    emp = [float(np.mean(x <= q)) for q in qs]
    from openmath.joint_dickman import dickman_marginal

    pred = [float(np.atleast_1d(dickman_marginal(q))[0]) for q in qs]

    return {
        "x": r(xs, 0),
        "density": r(dens),
        "gap": r(gaps),
        "q": r(qs, 3),
        "marginalEmpirical": r(emp),
        "marginalPredicted": r(pred),
        "nPairs": int(x.size),
        "limit": cfg.dickman_limit,
    }


# --------------------------------------------------------------------------- 021
def jacobsthal_series(cfg) -> dict:
    """Jacobsthal: primorial values against the quadratic envelope."""
    from openmath.jacobsthal import (
        MAX_PERIOD,
        h_lower_bound,
        iwaniec_envelope,
        primorial_jacobsthal,
        quadratic_envelope,
    )

    ks, vals = [], []
    for k in range(1, cfg.jacobsthal_k_max + 1):
        v, _ = primorial_jacobsthal(k)
        ks.append(k)
        vals.append(v)
        print(f"  021: j(P_{k}) = {v}")

    ks_arr = np.array(ks, dtype=float)
    vals_arr = np.array(vals, dtype=float)
    mask = ks_arr >= 3
    return {
        "k": ks,
        "j": r(vals, 0),
        "envelope": r(quadratic_envelope(ks_arr[mask])),
        "envelopeK": ks[2:],
        "iwaniec": r(iwaniec_envelope(ks_arr[mask])),
        "ratio": r(vals_arr[mask] / quadratic_envelope(ks_arr[mask])),
        "ratioIwaniec": r(vals_arr[mask] / iwaniec_envelope(ks_arr[mask])),
        "kMax": cfg.jacobsthal_k_max,
        "maxPeriod": float(MAX_PERIOD),
    }


# --------------------------------------------------------------------------- 026
def gaps_series(cfg) -> dict:
    """Prime gaps: how often the gap exceeds C log p, and how that settles."""
    from openmath.prime_gaps import cramer_reference, prime_gaps

    ps, gaps = prime_gaps(cfg.gap_limit)
    print(f"  026: {ps.size:,} primes")
    logp = np.log(ps[:-1].astype(float))

    n_blocks = 60
    edges = np.linspace(0, gaps.size, n_blocks + 1).astype(int)
    cs = np.linspace(0.25, 3.0, 24)
    curves = []
    prachar = []
    records = []
    for e in edges[1:]:
        g, lp = gaps[:e], logp[:e]
        curves.append(r([float(np.mean(g > c * lp)) for c in cs]))
        prachar.append(float(np.mean(g.astype(float) > ps[:e].astype(float) / np.arange(1, e + 1))))
        records.append(float(np.max(g)))
    xs = [float(e) for e in edges[1:]]

    return {
        "thresholds": r(cs, 3),
        "curves": curves,
        "prachar": r(prachar),
        "recordGap": r(records, 0),
        "cramer": r(cramer_reference(cs)),
        "x": r(xs, 0),
        "nPrimes": int(ps.size),
        "limit": cfg.gap_limit,
        "maxGap": int(gaps.max()),
    }


# --------------------------------------------------------------------------- 025
def egyptian_series(cfg) -> dict:
    """Egyptian fractions: N(b), the greedy comparison, and one worked example."""
    from openmath.egyptian import (
        double_log_envelope,
        greedy_expansion,
        max_min_length_table,
    )

    bs, table = max_min_length_table(cfg.egyptian_b_max, max_terms=cfg.egyptian_max_terms)
    N = table[0]
    print(f"  025: N(b) to b={cfg.egyptian_b_max}, max {N.max()}")

    # A worked greedy expansion: 5/121 needs three terms but greedy takes five.
    from fractions import Fraction

    a, b = 5, 121
    rest = Fraction(a, b)
    steps = []
    while rest > 0:
        n = -(-rest.denominator // rest.numerator)
        rest -= Fraction(1, n)
        steps.append({"den": n, "rest": float(rest), "sum": float(Fraction(a, b) - rest)})

    # The optimal expansion of the same fraction, for contrast.
    from openmath.egyptian import shortest_expansion

    best = shortest_expansion(a, b, max_terms=8)

    return {
        "b": bs.tolist(),
        "N": N.tolist(),
        "ratio": r([float(v) / float(np.log(np.log(b))) for v, b in zip(N, bs) if b >= 10]),
        "ratioB": [int(b) for b in bs if b >= 10],
        "envelope": r(double_log_envelope(bs[bs >= 10])),
        "worked": {"a": a, "b": b, "greedy": steps, "optimal": best},
        "maxN": int(N.max()),
        "bMax": cfg.egyptian_b_max,
    }


# --------------------------------------------------------------------------- 017
def pi_series(cfg) -> dict:
    """pi: the irrationality exponent as convergents accumulate."""
    import mpmath as mp

    from openmath.pi_exponent import (
        approximation_table,
        empirical_exponent,
        semiconvergents,
        to_float,
    )

    with mp.workdps(cfg.pi_dps):
        tab = approximation_table(mp.pi, terms=cfg.pi_terms, dps=cfg.pi_dps)
        s2 = approximation_table(mp.sqrt(2), terms=cfg.pi_terms, dps=cfg.pi_dps)

    keep = [i for i, q in enumerate(tab["q"]) if q >= 2]
    q = to_float([tab["q"][i] for i in keep])
    err = [tab["err"][i] for i in keep]
    quality = to_float([tab["squared_quality"][i] for i in keep])
    expo = empirical_exponent(err, [tab["q"][i] for i in keep])

    keep2 = [i for i, qq in enumerate(s2["q"]) if qq >= 2]
    q2 = to_float([s2["q"][i] for i in keep2])
    quality2 = to_float([s2["squared_quality"][i] for i in keep2])

    # The exponents dip below 2 at individual convergents, since the definition uses
    # the limsup.  A running maximum is a valid lower bound but is pinned forever by
    # the first small-q convergent, so it says nothing as q grows; a running median
    # over a short window shows the trend settling instead.
    w = 5
    running = np.array([
        float(np.median(expo[max(0, i - w + 1): i + 1])) for i in range(len(expo))
    ])

    # The first few semiconvergents of pi, to show that the best approximations of
    # the second kind also approach the bound.
    semis = []
    cf = [int(v) for v in tab["a"]]
    for idx in range(min(6, len(cf) - 1)):
        for p, qs in semiconvergents(cf, idx):
            with mp.workdps(cfg.pi_dps):
                e = abs(mp.pi - mp.mpf(p) / mp.mpf(qs))
            semis.append({"q": int(qs), "quality": float(mp.mpf(qs) ** 2 * e)})

    return {
        "q": r(np.log10(q)),
        "quality": r(quality),
        "qualitySqrt2": r(quality2),
        "qSqrt2": r(np.log10(q2)),
        "exponent": r(expo),
        "smoothedExponent": r(running),
        "partialQuotients": cf,
        "semiconvergents": [{"lq": r(np.log10(s["q"]), 3), "quality": r(s["quality"])} for s in semis],
        "maxQ": float(q.max()),
        "minQuality": float(quality.min()),
    }


# --------------------------------------------------------------------------- 003
def quasi_series(cfg) -> dict:
    """Zeta: the zero-free region and the zeros that fall outside it."""
    from openmath.quasi_riemann import (
        CLASSICAL_C,
        ZERO_FREE_BOUNDARY,
        classical_zero_free_region,
        zeta_zeros,
    )

    ts = zeta_zeros(cfg.zeta_zeros, dps=cfg.zeta_dps)
    print(f"  003: {ts.size} zeros, max t = {ts.max():.1f}")

    t_grid = np.linspace(10.0, float(ts.max()), 120)
    classical = classical_zero_free_region(t_grid, CLASSICAL_C)

    # The counting residual S(T) = k - (theta(t_k)/pi + 1), sampled along the zeros.
    import mpmath as mp

    with mp.workdps(cfg.zeta_dps):
        smooth = np.array(
            [float(mp.siegeltheta(float(t)) / mp.pi + 1) for t in ts]
        )
    index = np.arange(1, ts.size + 1, dtype=float)
    residual = index - smooth

    return {
        "t": r(ts),
        "critical": [0.5] * int(ts.size),
        "tGrid": r(t_grid, 2),
        "classical": r(classical),
        "boundary": ZERO_FREE_BOUNDARY,
        "residual": r(residual, 3),
        "count": int(ts.size),
        "violations": int(np.sum(np.full(ts.size, 0.5) > ZERO_FREE_BOUNDARY)),
        "tMax": float(ts.max()),
    }


BUILDERS = {
    "pd": pd_series,
    "dickman": dickman_series,
    "jacobsthal": jacobsthal_series,
    "gaps": gaps_series,
    "egyptian": egyptian_series,
    "pi": pi_series,
    "quasi": quasi_series,
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", default="", help="comma-separated keys, e.g. pi,quasi")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = cfgmod.load(Path(args.config) if args.config else None)
    keys = [s.strip() for s in args.only.split(",") if s.strip()] or list(BUILDERS)
    print(f"building {len(keys)} series: {', '.join(keys)}")

    # Reuse existing output whenever it is present, so that a partial rebuild
    # (``--only``) replaces just the named series and leaves the others in place.
    # Without this, running ``--only pi`` would silently drop every other series.
    out_path = DOCS / "data.js"
    payload: dict = {}
    if out_path.exists():
        text = out_path.read_text()
        try:
            payload = json.loads(text[text.index("{"): text.rindex("}") + 1])
        except (ValueError, IndexError):
            payload = {}

    for key in keys:
        if key not in BUILDERS:
            print(f"unknown key {key}", file=sys.stderr)
            continue
        t0 = time.time()
        print(f"[{key}]")
        payload[key] = BUILDERS[key](cfg)
        print(f"  {time.time() - t0:.1f}s")

    payload["meta"] = {
        "generated": time.strftime("%Y-%m-%d"),
        "source": "openai/math",
        "series": sorted(k for k in payload if k != "meta"),
    }
    out_path.write_text("window.OM_DATA = " + json.dumps(payload, separators=(",", ":")) + ";\n")
    print(f"wrote docs/data.js  ({out_path.stat().st_size / 1024:.0f} KiB)")


if __name__ == "__main__":
    main()
