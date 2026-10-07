"""Numerical limits for the probes.

All expensive parameters live here so that a run is reproducible from one file and
so that a reader can see exactly what range each figure covers.  The defaults are
chosen to finish in a few minutes on one laptop while keeping Monte Carlo and
truncation error well below the size of the effects being displayed.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

__all__ = ["Config", "load", "DEFAULT_PATH"]

DEFAULT_PATH = Path(__file__).resolve().parents[2] / "config.json"


@dataclass(frozen=True)
class Config:
    """Ranges and sample sizes for the seven probes."""

    # 011 Poisson--Dirichlet law for the prime predecessor.
    pd_prime_limit: int = 30_000_000
    pd_components: int = 8
    pd_simulations: int = 2_000_000
    pd_cdf_points: int = 120
    pd_windows: int = 6

    # 012 Joint Dickman law for consecutive integers.
    dickman_limit: int = 60_000_000
    dickman_grid: int = 40
    dickman_slices: int = 6

    # 021 Jacobsthal's function.
    jacobsthal_k_max: int = 10
    jacobsthal_extra: int = 2

    # 026 Positive lower density of large prime gaps.
    gap_limit: int = 200_000_000
    gap_thresholds: int = 40
    gap_stability_splits: int = 10

    # 025 Short Egyptian fractions.
    egyptian_b_max: int = 220
    egyptian_max_terms: int = 10

    # 017 Irrationality exponent of pi.
    pi_terms: int = 40
    pi_dps: int = 150

    # 003 Quasi-Riemann hypothesis.
    zeta_zeros: int = 600
    zeta_dps: int = 20
    character_q_max: int = 120
    character_samples: int = 30

    #: Seeds, one per stochastic probe, kept explicit for reproducibility.
    seeds: dict = field(default_factory=lambda: {"pd": 20260924, "char": 0})


def load(path: Path | None = None) -> Config:
    """Read a :class:`Config` from JSON, falling back to the defaults."""
    p = path or DEFAULT_PATH
    if not p.exists():
        return Config()
    data = json.loads(p.read_text())
    known = {f for f in Config.__dataclass_fields__}
    return Config(**{k: v for k, v in data.items() if k in known})


def dump(cfg: Config, path: Path | None = None) -> None:
    """Write ``cfg`` to JSON so a run can be reproduced exactly."""
    p = path or DEFAULT_PATH
    p.write_text(json.dumps(asdict(cfg), indent=2) + "\n")
