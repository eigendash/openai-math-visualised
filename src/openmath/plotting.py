"""Shared plotting conventions for the figures.

One place for the aesthetic choices, so that the seven figures read as one set: a
muted palette, serif type, thin spines, no chartjunk.  Matplotlib's cache directory
must be writable; set ``MPLCONFIGDIR`` if the default is not.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

__all__ = ["apply_style", "save", "colour", "FIGURE_DIR"]

FIGURE_DIR = Path(__file__).resolve().parents[2] / "results" / "figures"

#: A small qualitative palette, checked for print legibility.
_PALETTE = [
    "#2b4c7e",  # deep blue
    "#b4451f",  # rust
    "#3c7a5d",  # green
    "#8a6d1f",  # ochre
    "#6b3f74",  # plum
    "#2f6f77",  # teal
    "#8c8c8c",  # grey
]


def colour(i: int) -> str:
    """The ``i``-th palette colour, cycling."""
    return _PALETTE[i % len(_PALETTE)]


def apply_style() -> None:
    """Install the shared rcParams.  Idempotent."""
    plt.rcParams.update(
        {
            "figure.dpi": 130,
            "savefig.dpi": 200,
            "savefig.bbox": "tight",
            "font.family": "serif",
            "font.serif": ["DejaVu Serif", "Times New Roman", "Times"],
            "mathtext.fontset": "dejavuserif",
            "font.size": 9.5,
            "axes.titlesize": 10.5,
            "axes.titleweight": "normal",
            "axes.labelsize": 9.5,
            "axes.grid": True,
            "axes.axisbelow": True,
            "axes.edgecolor": "#4a4a4a",
            "axes.linewidth": 0.7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "grid.color": "#d8d8d8",
            "grid.linewidth": 0.6,
            "legend.frameon": False,
            "legend.fontsize": 8.5,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "xtick.color": "#333333",
            "ytick.color": "#333333",
        }
    )


def save(fig, name: str, directory: Path | None = None) -> Path:
    """Save ``fig`` as ``<name>.pdf`` and ``<name>.png`` under the figure directory."""
    out = directory or FIGURE_DIR
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{name}.png"
    fig.savefig(path)
    fig.savefig(out / f"{name}.pdf")
    plt.close(fig)
    return path
