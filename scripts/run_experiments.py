#!/usr/bin/env python
"""Run the seven numerical probes and write the figures and results tables.

Usage::

    python scripts/run_experiments.py [--only 011,026] [--config config.json]

Every probe writes a JSON record of its summary statistics under ``results/`` and one
or two figures under ``results/figures/``.  Nothing is asserted about the truth of the
manuscripts; the script reports what the computation found.
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
from openmath.plotting import apply_style, colour, save  # noqa: E402

RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)


def _write(name: str, payload: dict) -> None:
    (RESULTS / f"{name}.json").write_text(json.dumps(payload, indent=2, default=str) + "\n")
    print(f"  wrote results/{name}.json")


# --------------------------------------------------------------------------- 011
def probe_011(cfg) -> None:
    """Poisson--Dirichlet law for the prime predecessor."""
    import matplotlib.pyplot as plt

    from openmath.poisson_dirichlet import (
        dickman,
        normalized_log_factors,
        pd1_order_statistic_cdf,
        predecessor_factors,
        stick_breaking_pd1,
    )

    t0 = time.time()
    k = cfg.pd_components
    ps, facs = predecessor_factors(cfg.pd_prime_limit)
    V = np.array([normalized_log_factors(f, int(p), k) for f, p in zip(facs, ps)])
    print(f"  011: {ps.size} primes <= {cfg.pd_prime_limit:,} factored")

    # The empirical distribution of the j-th largest normalised log factor.
    scores = V[:, :k]
    # Moments against the simulated law.
    sim = stick_breaking_pd1(min(cfg.pd_simulations, 400_000), top=k, seed=cfg.seeds["pd"])
    emp_mean = scores.mean(axis=0)
    sim_mean = sim.mean(axis=0)
    emp_sum = float(scores.sum(axis=1).mean())
    sim_sum = float(sim.sum(axis=1).mean())

    # Kolmogorov--Smirnov distance per component, using the simulation as reference.
    ks = np.array(
        [
            float(
                np.max(
                    np.abs(
                        np.searchsorted(np.sort(sim[:, j]), np.linspace(0, 1, 2000)) / sim.shape[0]
                        - np.searchsorted(np.sort(scores[:, j]), np.linspace(0, 1, 2000)) / scores.shape[0]
                    )
                )
            )
            for j in range(k)
        ]
    )

    # Convergence in x: the mean of the largest component over growing windows.
    windows = np.array_split(np.arange(ps.size), cfg.pd_windows)
    wins_x = np.array([float(ps[w[-1]]) for w in windows])
    wins_mean = np.array([float(scores[w, 0].mean()) for w in windows])

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.3))

    ax = axes[0]
    grid = np.linspace(0.005, 0.999, 400)
    for j in (1, 2, 3, 4):
        ax.plot(
            grid,
            pd1_order_statistic_cdf(grid, j, samples=400_000, seed=cfg.seeds["pd"]),
            color=colour(j),
            lw=1.2,
            label=rf"$j={j}$",
        )
        ax.plot(
            np.sort(scores[:, j - 1]),
            np.arange(1, scores.shape[0] + 1) / scores.shape[0],
            color=colour(j),
            lw=1.0,
            ls="--",
            alpha=0.9,
        )
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$\Pr[V_j \leq x]$")
    ax.set_title("Order statistics against PD(1)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(title="solid: PD(1)\ndashed: primes", loc="lower right")

    ax = axes[1]
    for j in range(4):
        ax.hist(
            scores[:, j],
            bins=90,
            range=(0, 1),
            density=True,
            histtype="step",
            lw=1.0,
            color=colour(j),
            label=rf"$V_{{{j + 1}}}$",
        )
    ax.set_xlabel("normalised log factor")
    ax.set_ylabel("density")
    ax.set_title(r"Prime predecessors, $p \leq 3\times10^{7}$")
    ax.legend()

    ax = axes[2]
    ax.plot(wins_x, wins_mean, "o-", color=colour(0), lw=1.2, ms=4, label=r"$\overline{V_1}$")
    ax.axhline(sim_mean[0], color=colour(1), lw=1.2, ls="--", label=r"PD(1) mean")
    ax.axhline(float(dickman(1.0)), color=colour(2), lw=1.0, ls=":", label=r"$\rho(1)=1$")
    ax.set_xscale("log")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"mean of $V_1$")
    ax.set_title(r"Convergence of the largest component")
    ax.legend()

    fig.tight_layout()
    save(fig, "011-poisson-dirichlet")

    _write(
        "011-poisson-dirichlet",
        {
            "prime_limit": cfg.pd_prime_limit,
            "primes": int(ps.size),
            "components": k,
            "empirical_mean": emp_mean.tolist(),
            "pd1_simulated_mean": sim_mean.tolist(),
            "empirical_sum_top_k": emp_sum,
            "pd1_sum_top_k": sim_sum,
            "ks_distance": ks.tolist(),
            "window_x": wins_x.tolist(),
            "window_mean_V1": wins_mean.tolist(),
            "seconds": round(time.time() - t0, 1),
        },
    )


# --------------------------------------------------------------------------- 012
def probe_012(cfg) -> None:
    """Joint Dickman law for consecutive integers."""
    import matplotlib.pyplot as plt

    from openmath.joint_dickman import (
        comparison_density,
        dickman_marginal,
        independence_gap,
        joint_cdf_table,
        largest_prime_factors,
        normalized_coordinates,
    )

    t0 = time.time()
    lpf = largest_prime_factors(cfg.dickman_limit)
    x, y = normalized_coordinates(lpf)
    dens = comparison_density(x, y)
    print(f"  012: {x.size:,} pairs, density(P+(n)<P+(n+1)) = {dens:.5f}")

    grid = np.linspace(0.0, 1.0, cfg.dickman_grid + 1)[1:]
    joint, product, counts = joint_cdf_table(x, y, grid)
    gap = independence_gap(joint, product)

    # Marginal against the predicted rho(1/a).
    qs = np.linspace(0.05, 0.95, 25)
    emp_marginal = np.array([float(np.mean(x <= q)) for q in qs])
    predicted = dickman_marginal(qs)

    # Density in growing windows, to display the approach to 1/2.
    windows = np.array_split(np.arange(x.size), cfg.dickman_slices)
    wx = np.array([float(w[-1] + 2) for w in windows])
    wd = np.array([comparison_density(x[w], y[w]) for w in windows])

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.3))

    ax = axes[0]
    im = ax.pcolormesh(grid, grid, np.diff(np.diff(joint, axis=0, prepend=0), axis=1, prepend=0), cmap="magma_r", shading="auto")
    ax.set_xlabel(r"$x_n = \log P^+(n)/\log n$")
    ax.set_ylabel(r"$y_n = \log P^+(n+1)/\log n$")
    ax.set_title("Joint density of the two coordinates")
    fig.colorbar(im, ax=ax, fraction=0.046, label="density")

    ax = axes[1]
    ax.plot([0, 1], [0, 1], color="#888888", lw=0.8, ls=":")
    ax.plot(grid, grid - np.array([joint[i, i] for i in range(grid.size)]), color=colour(0), lw=1.3, label="empirical joint CDF")
    ax.plot(grid, grid - np.array([product[i, i] for i in range(grid.size)]), color=colour(1), lw=1.3, ls="--", label="product of marginals")
    ax.set_xlabel(r"$a = b$")
    ax.set_ylabel("CDF along the diagonal")
    ax.set_title(f"Independence: sup gap {gap:.3f}")
    ax.legend()

    ax = axes[2]
    ax.plot(qs, emp_marginal, "o", color=colour(0), ms=3.5, label=r"empirical $\Pr[x_n \leq a]$")
    ax.plot(qs, predicted, color=colour(1), lw=1.4, label=r"$\rho(1/a)$")
    ax.set_xlabel(r"$a$")
    ax.set_ylabel("marginal CDF")
    ax.set_title("Marginal against the Dickman law")
    ax.legend()

    fig.tight_layout()
    save(fig, "012-joint-dickman")

    _write(
        "012-joint-dickman",
        {
            "limit": cfg.dickman_limit,
            "pairs": int(x.size),
            "density_P_less": dens,
            "independence_gap": gap,
            "window_x": wx.tolist(),
            "window_density": wd.tolist(),
            "marginal_q": qs.tolist(),
            "marginal_empirical": emp_marginal.tolist(),
            "marginal_predicted": predicted.tolist(),
            "seconds": round(time.time() - t0, 1),
        },
    )


# --------------------------------------------------------------------------- 021
def probe_021(cfg) -> None:
    """Jacobsthal's function and the quadratic bound."""
    import matplotlib.pyplot as plt

    from openmath.jacobsthal import (
        MAX_PERIOD,
        h_lower_bound,
        iwaniec_envelope,
        primorial_jacobsthal,
        quadratic_envelope,
    )

    t0 = time.time()
    ks, vals = [], []
    for k in range(1, cfg.jacobsthal_k_max + 1):
        v, _ = primorial_jacobsthal(k)
        ks.append(k)
        vals.append(v)
        print(f"  021: j(P_{k}) = {v}")
    ks = np.array(ks, dtype=float)
    vals = np.array(vals, dtype=float)

    # Implied constant C in h(k) <= C k^2 / (log log 3k)^2, over k >= 3.
    mask = ks >= 3
    ratio = vals[mask] / quadratic_envelope(ks[mask])
    iwaniec_ratio = vals[mask] / iwaniec_envelope(ks[mask])

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.4))

    ax = axes[0]
    kk = np.linspace(3, cfg.jacobsthal_k_max, 200)
    ax.plot(ks, vals, "o-", color=colour(0), lw=1.3, ms=4.5, label=r"$j(P_k)$, computed")
    ax.plot(kk, quadratic_envelope(kk), color=colour(1), lw=1.2, label=r"$k^2/(\log\log 3k)^2$")
    ax.plot(kk, iwaniec_envelope(kk), color=colour(2), lw=1.2, ls="--", label=r"$k^2\log^2 k$ (Iwaniec)")
    ax.set_xlabel(r"$k$")
    ax.set_ylabel(r"$h(k)$ lower bound")
    ax.set_yscale("log")
    ax.set_title("Jacobsthal values against the two bounds")
    ax.legend(loc="upper left", fontsize=7.5)

    ax = axes[1]
    ax.plot(ks[mask], ratio, "o-", color=colour(0), lw=1.3, ms=4.5, label="vs new bound")
    ax.plot(ks[mask], iwaniec_ratio, "s--", color=colour(2), lw=1.3, ms=4, label="vs Iwaniec bound")
    ax.set_xlabel(r"$k$")
    ax.set_ylabel("value / envelope")
    ax.set_title("Implied constant in each bound")
    ax.legend()

    fig.tight_layout()
    save(fig, "021-jacobsthal")

    _write(
        "021-jacobsthal",
        {
            "k": ks.tolist(),
            "primorial_j": vals.tolist(),
            "quadratic_envelope": quadratic_envelope(ks[mask]).tolist(),
            "implied_C": ratio.tolist(),
            "implied_C_iwaniec": iwaniec_ratio.tolist(),
            "max_period": MAX_PERIOD,
            "seconds": round(time.time() - t0, 1),
        },
    )


# --------------------------------------------------------------------------- 026
def probe_026(cfg) -> None:
    """Positive lower density of large prime gaps."""
    import matplotlib.pyplot as plt

    from openmath.prime_gaps import (
        cramer_reference,
        density_stability,
        exceedance_curve,
        prachar_density,
        prime_gaps,
    )

    t0 = time.time()
    ps, gaps = prime_gaps(cfg.gap_limit)
    print(f"  026: {ps.size:,} primes <= {cfg.gap_limit:,}")

    cs = np.linspace(0.25, 3.0, cfg.gap_thresholds)
    dens = exceedance_curve(ps, gaps, cs)
    prachar = prachar_density(ps, gaps)
    print(f"  026: Prachar density = {prachar:.5f}")

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.3))

    ax = axes[0]
    ax.plot(cs, dens, "-", color=colour(0), lw=1.4, label="measured")
    ax.plot(cs, cramer_reference(cs), color=colour(2), lw=1.1, ls=":", label=r"Cramér $e^{-C}$")
    ax.set_xlabel(r"$C$")
    ax.set_ylabel(r"proportion with $d_n > C\log p_n$")
    ax.set_title("Exceedance density")
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 1.0)
    ax.legend()

    ax = axes[1]
    logp = np.log(ps[:-1].astype(float))
    r = gaps / logp
    ax.hist(r[r < 6], bins=120, range=(0, 6), color=colour(0), histtype="step", lw=1.1)
    for c, lab in [(1.0, r"$C=1$"), (2.0, r"$C=2$")]:
        ax.axvline(c, color=colour(1), lw=1.0, ls="--")
    ax.set_yscale("log")
    ax.set_xlabel(r"$d_n / \log p_n$")
    ax.set_ylabel("count")
    ax.set_title("Normalised gap distribution")

    ax = axes[2]
    for i, c in enumerate([0.5, 1.0, 2.0]):
        st = density_stability(ps, gaps, c, splits=cfg.gap_stability_splits)
        ax.plot(np.arange(1, st.size + 1) / st.size, st, "o-", color=colour(i), lw=1.2, ms=3.5, label=rf"$C={c}$")
    ax.set_xlabel("fraction of the range")
    ax.set_ylabel("density within the block")
    ax.set_title("Stability across the range")
    ax.legend()

    fig.tight_layout()
    save(fig, "026-prime-gaps")

    _write(
        "026-prime-gaps",
        {
            "limit": cfg.gap_limit,
            "primes": int(ps.size),
            "thresholds": cs.tolist(),
            "exceedance": dens.tolist(),
            "prachar_density": prachar,
            "max_gap": int(gaps.max()),
            "max_gap_ratio": float((gaps / np.log(ps[:-1])).max()),
            "seconds": round(time.time() - t0, 1),
        },
    )


# --------------------------------------------------------------------------- 025
def probe_025(cfg) -> None:
    """Short Egyptian fractions."""
    import matplotlib.pyplot as plt

    from openmath.egyptian import (
        double_log_envelope,
        erdos_graham_envelope,
        greedy_length_table,
        max_min_length_table,
        vose_envelope,
    )

    t0 = time.time()
    bs, table = max_min_length_table(cfg.egyptian_b_max, max_terms=cfg.egyptian_max_terms)
    N = table[0]
    arg = table[1]
    print(f"  025: N(b) up to b={cfg.egyptian_b_max}, max N = {N.max()}")

    gbs, greedy = greedy_length_table(cfg.egyptian_b_max, step=4)
    mask = bs >= 10
    ratio = N[mask] / double_log_envelope(bs[mask])

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.3))

    ax = axes[0]
    ax.step(bs, N, where="post", color=colour(0), lw=1.2, label=r"$N(b)=\max_a N(a,b)$")
    ax.plot(gbs, greedy, color=colour(1), lw=1.1, ls="--", label="worst greedy length")
    ax.set_xlabel(r"$b$")
    ax.set_ylabel("number of terms")
    ax.set_title("Minimum length against the greedy algorithm")
    ax.legend()

    ax = axes[1]
    bb = np.linspace(10, cfg.egyptian_b_max, 300)
    ax.plot(bs, N, color=colour(0), lw=1.3, label=r"$N(b)$")
    ax.plot(bb, 2.2 * double_log_envelope(bb), color=colour(1), lw=1.2, ls="--", label=r"$c_2\log\log b$")
    ax.plot(bb, erdos_graham_envelope(bb), color=colour(2), lw=1.1, ls=":", label=r"$\log b/\log\log b$")
    ax.plot(bb, vose_envelope(bb), color=colour(3), lw=1.1, ls="-.", label=r"$\sqrt{\log b}$")
    ax.set_xlabel(r"$b$")
    ax.set_ylabel("terms")
    ax.set_title("Growth against the bounds")
    ax.legend(fontsize=7.5)

    ax = axes[2]
    ax.plot(bs[mask], ratio, color=colour(0), lw=1.3)
    ax.set_xlabel(r"$b$")
    ax.set_ylabel(r"$N(b)/\log\log b$")
    ax.set_title("Empirical constant in the double-logarithmic order")
    ax.set_ylim(0, max(6.0, float(ratio.max()) * 1.15))

    fig.tight_layout()
    save(fig, "025-egyptian-fractions")

    _write(
        "025-egyptian-fractions",
        {
            "b_max": cfg.egyptian_b_max,
            "b": bs.tolist(),
            "N": N.tolist(),
            "argmax": arg.tolist(),
            "max_N": int(N.max()),
            "implied_c2": ratio.tolist(),
            "seconds": round(time.time() - t0, 1),
        },
    )


# --------------------------------------------------------------------------- 017
def probe_017(cfg) -> None:
    """The irrationality exponent of pi."""
    import matplotlib.pyplot as plt
    import mpmath as mp

    from openmath.pi_exponent import approximation_table, empirical_exponent, semiconvergents

    t0 = time.time()
    with mp.workdps(cfg.pi_dps):
        table = approximation_table(mp.pi, terms=cfg.pi_terms, dps=cfg.pi_dps)
    from openmath.pi_exponent import to_float

    a = table["a"]
    # The first convergent of any x > 1 has q = 1; the definition of the
    # irrationality exponent requires q >= 2, so that index is dropped.
    keep = [i for i, qi in enumerate(table["q"]) if qi >= 2]
    q = to_float([table["q"][i] for i in keep])
    err = [table["err"][i] for i in keep]
    sq = to_float([table["squared_quality"][i] for i in keep])
    exponents = empirical_exponent(err, [table["q"][i] for i in keep])
    print(f"  017: {len(a)} convergents, max q = {q.max():.3e}")

    best = float(np.min(sq))
    # Catalan's constant for comparison: mu is not known to be 2.
    with mp.workdps(cfg.pi_dps):
        cat = mp.catalan
        cat_table = approximation_table(cat, terms=cfg.pi_terms, dps=cfg.pi_dps)
        sqrt2_table = approximation_table(mp.sqrt(2), terms=cfg.pi_terms, dps=cfg.pi_dps)

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.3))

    ax = axes[0]
    ax.loglog(q, sq, "o-", color=colour(0), lw=1.1, ms=3, label=r"$\pi$")
    for tbl, sty, col, lab in (
        (sqrt2_table, "s--", colour(2), r"$\sqrt{2}$"),
        (cat_table, "^:", colour(1), r"$G$ (Catalan)"),
    ):
        kk = [i for i, qi in enumerate(tbl["q"]) if qi >= 2]
        ax.loglog(to_float([tbl["q"][i] for i in kk]), to_float([tbl["squared_quality"][i] for i in kk]), sty, color=col, lw=1.0, ms=2.5, label=lab)
    ax.set_xlabel(r"$q$")
    ax.set_ylabel(r"$q^2\,|x-p/q|$")
    ax.set_title(r"Squared quality of convergents")
    ax.legend()

    ax = axes[1]
    ax.plot(q, exponents, "o-", color=colour(0), lw=1.1, ms=3, label=r"$\pi$")
    s2 = [i for i, qi in enumerate(sqrt2_table["q"]) if qi >= 2]
    ax.plot(to_float([sqrt2_table["q"][i] for i in s2]), empirical_exponent([sqrt2_table["err"][i] for i in s2], [sqrt2_table["q"][i] for i in s2]), "s--", color=colour(2), lw=1.0, ms=2.5, label=r"$\sqrt{2}$")
    ax.axhline(2.0, color=colour(1), lw=1.2, ls="--", label=r"$\nu=2$")
    ax.set_xscale("log")
    ax.set_xlabel(r"$q$")
    ax.set_ylabel(r"$-\log|x-p/q|/\log q$")
    ax.set_title("Empirical exponent along convergents")
    ax.legend()

    ax = axes[2]
    idx = int(np.argmax(a[1:]) + 1) if len(a) > 1 else 0
    ax.bar(np.arange(1, len(a) + 1), a, color=colour(0), width=0.8)
    ax.set_yscale("log")
    ax.set_xlabel("index")
    ax.set_ylabel("partial quotient")
    ax.set_title("Partial quotients of $\\pi$")

    fig.tight_layout()
    save(fig, "017-pi-exponent")

    _write(
        "017-pi-exponent",
        {
            "terms": int(len(a)),
            "partial_quotients": a.tolist(),
            "max_partial_quotient": int(a.max()),
            "max_q": float(q.max()),
            "min_squared_quality": best,
            "top_exponents": np.sort(exponents)[:10].tolist(),
            "seconds": round(time.time() - t0, 1),
        },
    )


# --------------------------------------------------------------------------- 003
def probe_003(cfg) -> None:
    """Quasi-Riemann hypothesis: the fixed zero-free half-plane."""
    import matplotlib.pyplot as plt

    from openmath.quasi_riemann import (
        CLASSICAL_C,
        ZERO_FREE_BOUNDARY,
        classical_zero_free_region,
        zero_count_check,
        zero_real_parts,
        zeta_zeros,
    )

    t0 = time.time()
    ts = zeta_zeros(cfg.zeta_zeros, dps=cfg.zeta_dps)
    re_parts = zero_real_parts(cfg.zeta_zeros, dps=cfg.zeta_dps)
    count = zero_count_check(float(ts.max()), cfg.zeta_zeros, dps=cfg.zeta_dps)
    print(f"  003: {ts.size} zeros, max Re = {re_parts.max():.6f}, t_max = {ts.max():.1f}")

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.3))

    ax = axes[0]
    tv = np.linspace(10.0, float(ts.max()), 400)
    ax.plot(classical_zero_free_region(tv, CLASSICAL_C), tv, color=colour(2), lw=1.3, label=r"classical $1-c/\log t$")
    ax.axvline(ZERO_FREE_BOUNDARY, color=colour(1), lw=1.6, label=r"$\mathrm{Re}\,s=7/8$ (theorem)")
    ax.axvline(1.0, color="#888888", lw=1.0, ls=":", label=r"$\mathrm{Re}\,s=1$")
    ax.plot(0.5, 0.0, "o", color=colour(0), ms=4, label=r"zeros $\mathrm{Re}\,\rho=1/2$")
    ax.set_xlim(0.45, 1.03)
    ax.set_xlabel(r"$\sigma = \mathrm{Re}\,s$")
    ax.set_ylabel(r"$t = \mathrm{Im}\,s$")
    ax.set_title("Zero-free region: classical shape against a fixed line")
    ax.legend(loc="upper left", fontsize=7)

    ax = axes[1]
    # The k-th zero has height t_k, and the smooth Riemann--von Mangoldt count
    # theta(T)/pi + 1 should be close to k.  The residual is the classical S(T).
    smooth = mp_siegel_count(ts)
    index = np.arange(1, ts.size + 1, dtype=float)
    ax.plot(ts, index - smooth, color=colour(0), lw=1.0)
    ax.axhline(0.0, color="#888888", lw=0.8, ls=":")
    ax.set_xlabel(r"$T$")
    ax.set_ylabel(r"$k - (\theta(t_k)/\pi + 1)$")
    ax.set_title(r"Riemann--von Mangoldt residual $S(T)$")
    ax.set_ylim(-3.0, 3.0)

    ax = axes[2]
    sig = np.linspace(0.4, 0.95, 60)
    dens = np.array([float(np.mean(re_parts > s)) for s in sig])
    ax.plot(sig, dens, color=colour(0), lw=1.4)
    ax.axvline(ZERO_FREE_BOUNDARY, color=colour(1), lw=1.4, ls="--", label=r"$7/8$")
    ax.axvline(0.5, color=colour(2), lw=1.0, ls=":", label=r"$1/2$")
    ax.set_xlabel(r"$\sigma$")
    ax.set_ylabel(r"fraction of zeros with $\mathrm{Re}\,\rho > \sigma$")
    ax.set_title("Empirical zero density")
    ax.set_ylim(-0.02, 1.02)
    ax.legend()

    fig.tight_layout()
    save(fig, "003-quasi-riemann")

    _write(
        "003-quasi-riemann",
        {
            "zeros": int(ts.size),
            "t_max": float(ts.max()),
            "max_real_part": float(re_parts.max()),
            "violations_of_7_8": int(np.sum(re_parts > ZERO_FREE_BOUNDARY)),
            "dps": cfg.zeta_dps,
            "zero_count": count,
            "seconds": round(time.time() - t0, 1),
        },
    )


def mp_siegel_count(ts: np.ndarray) -> np.ndarray:
    """The smooth Riemann--von Mangoldt count ``theta(T)/pi + 1`` at each height."""
    import mpmath as mp

    with mp.workdps(20):
        return np.array([float(mp.siegeltheta(float(t)) / mp.pi + 1) for t in ts])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", default="", help="comma-separated probe ids, e.g. 011,026")
    ap.add_argument("--config", default=None, help="path to a JSON config")
    args = ap.parse_args()

    cfg = cfgmod.load(Path(args.config) if args.config else None)
    apply_style()

    probes = {
        "011": probe_011,
        "012": probe_012,
        "021": probe_021,
        "026": probe_026,
        "025": probe_025,
        "017": probe_017,
        "003": probe_003,
    }
    wanted = [s.strip() for s in args.only.split(",") if s.strip()] or list(probes)
    for pid in wanted:
        if pid not in probes:
            print(f"unknown probe {pid}", file=sys.stderr)
            continue
        print(f"[{pid}]")
        probes[pid](cfg)


if __name__ == "__main__":
    main()
