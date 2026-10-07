"""Configuration round-trip and the plotting helpers."""

from __future__ import annotations

import json

import numpy as np
import pytest

from openmath import config as cfgmod
from openmath.plotting import apply_style, colour, save


def test_config_defaults_are_sane():
    cfg = cfgmod.Config()
    assert cfg.pd_prime_limit > 1000
    assert cfg.pd_components >= 4
    assert cfg.gap_limit > 10**6
    assert cfg.egyptian_b_max >= 50
    assert cfg.pi_terms >= 10
    assert cfg.zeta_zeros >= 10


def test_config_round_trip(tmp_path):
    cfg = cfgmod.Config(pd_prime_limit=1234, gap_limit=5678)
    path = tmp_path / "config.json"
    cfgmod.dump(cfg, path)
    loaded = cfgmod.load(path)
    assert loaded.pd_prime_limit == 1234
    assert loaded.gap_limit == 5678
    assert loaded.egyptian_b_max == cfg.egyptian_b_max


def test_config_load_ignores_unknown_keys(tmp_path):
    path = tmp_path / "config.json"
    path.write_text(json.dumps({"pd_prime_limit": 99, "not_a_field": 1}))
    assert cfgmod.load(path).pd_prime_limit == 99


def test_config_load_missing_file_returns_defaults(tmp_path):
    assert cfgmod.load(tmp_path / "absent.json") == cfgmod.Config()


def test_colour_cycles():
    assert colour(0) == colour(0)
    assert colour(0) != colour(1)
    assert colour(0) == colour(7)


def test_apply_style_is_idempotent():
    apply_style()
    apply_style()
    import matplotlib.pyplot as plt

    assert plt.rcParams["axes.spines.top"] is False


def test_save_writes_png_and_pdf(tmp_path):
    import matplotlib.pyplot as plt

    apply_style()
    fig, ax = plt.subplots()
    ax.plot(np.arange(5), np.arange(5))
    path = save(fig, "unit-test", directory=tmp_path)
    assert path.name == "unit-test.png"
    assert path.exists()
    assert (tmp_path / "unit-test.pdf").exists()
