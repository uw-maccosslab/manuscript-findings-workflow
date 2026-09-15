"""Render the group x cohort confound as an annotated sample-count heatmap.

A publication-ready figure (dual-exported SVG + 300 DPI PNG, plus a separate
colorbar legend image) that exposes the aliasing between diagnostic group and
collection cohort/site in the 40-sample experimental subset. Rows are diagnostic
group (ADD, PDD, PDCN, HCN); columns are collection cohort (Kerchner, Poston,
ADRC/PUC, PUC only). Each cell shows the integer sample count, shaded on a
perceptually-uniform sequential colormap (viridis) because this is a count
matrix, not a categorical encoding. Row labels carry their registry ``group``
color and column labels their registry ``cohort`` color (direct labels, so the
only separate legend is the count colorbar).

The block structure it must make visible: ADD is 10/10 Kerchner; PDCN is 10/10
Poston; the whole PD spectrum (PDD + PDCN) sits in Poston/PUC while the rest is
largely Kerchner -- so diagnostic group is almost perfectly aliased with site.

Everything is computed from data + the confound-stats JSON so the figure
regenerates from data alone (no hard-coded counts or statistics):

Inputs (all read-only, regenerable):
  data/metadata.csv            raw metadata (source of the crosstab; hashed for
                               the data version)
  results/confounding_stats.json  bias-corrected Cramer's V values (title/annotation)
  state/color_registry.json    group / cohort category -> color map

Run:  ./.venv/bin/python scripts/figures/plot_group_x_cohort_confound.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

# Import the project's shared figure helpers (seeded from the lib/ templates).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.colors import get_color
from common.figure_io import publication_style, save_figure

LOGGER = logging.getLogger("plot_group_x_cohort_confound")

# Fixed by the figure spec: ADD is the primary group of interest and leads; the
# PD spectrum (PDD, PDCN) is kept adjacent so its single-site block reads as one.
GROUP_ORDER: tuple[str, ...] = ("ADD", "PDD", "PDCN", "HCN")
COHORT_ORDER: tuple[str, ...] = ("Kerchner", "Poston", "ADRC/PUC", "PUC only")
EXPERIMENTAL_GROUPS = frozenset(GROUP_ORDER)

DEFAULT_OUTPUT_DIR = Path("figures") / "metadata" / "crosstabs"
DEFAULT_BASE_NAME = "group-x-cohort-confound"
DEFAULT_DATA = Path("data") / "metadata.csv"
DEFAULT_CONFOUNDS = Path("results") / "confounding_stats.json"
DEFAULT_REGISTRY = Path("state") / "color_registry.json"

CMAP_NAME = "viridis"


def data_version(path: Path) -> str:
    """Return a short content hash of the raw metadata, for provenance."""
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16]


def build_crosstab(data_path: Path) -> pd.DataFrame:
    """Compute the experimental group x cohort count matrix, fail-loud on gaps.

    Filters to the four experimental diagnostic groups (control pools excluded)
    and returns a DataFrame indexed by ``GROUP_ORDER`` with ``COHORT_ORDER``
    columns. Raises if a value falls outside the expected group/cohort sets, so
    an unexpected label can never be silently dropped from the picture.
    """
    df = pd.read_csv(data_path)
    for col in ("Condition", "Cohort"):
        if col not in df.columns:
            raise ValueError(f"{data_path} is missing the {col!r} column.")

    exp = df[df["Condition"].isin(EXPERIMENTAL_GROUPS)].copy()
    if len(exp) != 40:
        raise ValueError(
            f"Expected 40 experimental samples (ADD/PDD/PDCN/HCN); got {len(exp)}."
        )

    unexpected_cohorts = set(exp["Cohort"].unique()) - set(COHORT_ORDER)
    if unexpected_cohorts:
        raise ValueError(
            f"Unexpected cohort label(s) in experimental subset: {unexpected_cohorts}."
        )

    ct = pd.crosstab(exp["Condition"], exp["Cohort"])
    ct = (
        ct.reindex(index=GROUP_ORDER, columns=COHORT_ORDER).fillna(0).astype(int)
    )
    total = int(ct.to_numpy().sum())
    if total != 40:
        raise ValueError(f"Crosstab totals {total}, expected 40.")
    return ct


def load_cramers_v(confounds_path: Path) -> dict[str, float]:
    """Read the three bias-corrected Cramer's V values used in the annotation."""
    stats = json.loads(confounds_path.read_text())
    keys = ("group_x_cohort_V", "ADD_vs_other_cohort_V", "PD_vs_other_cohort_V")
    missing = [k for k in keys if k not in stats]
    if missing:
        raise ValueError(f"{confounds_path} is missing keys: {missing}.")
    return {k: float(stats[k]) for k in keys}


def make_figure(
    ct: pd.DataFrame, v_stats: dict[str, float], registry_path: Path
) -> tuple[plt.Figure, plt.Figure]:  # type: ignore[name-defined]
    """Draw the annotated heatmap and its companion colorbar legend figure."""
    counts = ct.to_numpy()
    n_rows, n_cols = counts.shape
    vmax = int(counts.max())
    norm = Normalize(vmin=0, vmax=vmax)
    cmap = matplotlib.colormaps[CMAP_NAME]

    row_colors = [
        get_color("group", g, registry_path=registry_path, persist=False)
        for g in GROUP_ORDER
    ]
    col_colors = [
        get_color("cohort", c, registry_path=registry_path, persist=False)
        for c in COHORT_ORDER
    ]

    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    ax.imshow(counts, cmap=cmap, norm=norm, aspect="equal")

    # Per-cell integer counts. Zero cells sit at the dark low end of viridis but
    # carry no samples, so their "0" is muted; non-zero cells get an adaptive
    # black/white glyph for contrast against the sequential fill.
    for i in range(n_rows):
        for j in range(n_cols):
            value = int(counts[i, j])
            if value == 0:
                ax.text(
                    j, i, "0", ha="center", va="center", fontsize=11,
                    color="#9aa0a6",
                )
            else:
                r, g, b, _ = cmap(norm(value))
                luminance = 0.299 * r + 0.587 * g + 0.114 * b
                ax.text(
                    j, i, str(value), ha="center", va="center",
                    fontsize=14, fontweight="bold",
                    color="black" if luminance > 0.6 else "white",
                )

    # Row marginal totals to the right of the grid -- the load-bearing n that
    # shows every diagnostic group is balanced at 10 samples. (Column totals are
    # left to the caption to avoid crowding the cohort tick labels below.)
    row_tot = counts.sum(axis=1)
    for i, t in enumerate(row_tot):
        ax.text(
            n_cols - 0.30, i, f"{t}", ha="left", va="center",
            fontsize=10.5, color="#333333", fontstyle="italic",
        )

    # Ticks: colored category labels act as the direct key (no categorical legend).
    ax.set_xticks(range(n_cols))
    ax.set_yticks(range(n_rows))
    ax.set_xticklabels(COHORT_ORDER)
    ax.set_yticklabels(GROUP_ORDER)
    # A thin dark halo (path effect) keeps every label legible on white while the
    # fill stays the exact registry color -- notably the pale ADRC/PUC yellow
    # (#F0E442), which would otherwise fall below the WCAG 3:1 contrast floor.
    halo = [path_effects.withStroke(linewidth=1.6, foreground="#222222")]
    for tick, color in zip(ax.get_xticklabels(), col_colors, strict=True):
        tick.set_color(color)
        tick.set_fontweight("bold")
        tick.set_path_effects(halo)
    for tick, color in zip(ax.get_yticklabels(), row_colors, strict=True):
        tick.set_color(color)
        tick.set_fontweight("bold")
        tick.set_path_effects(halo)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Thin white gridlines separate the cells so the blocks read cleanly.
    ax.set_xticks([x - 0.5 for x in range(n_cols + 1)], minor=True)
    ax.set_yticks([y - 0.5 for y in range(n_rows + 1)], minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", length=0)

    ax.set_xlabel("Collection cohort / site", labelpad=10)
    ax.set_ylabel("Diagnostic group", labelpad=10)

    ax.set_title(
        "Diagnostic group aliases with collection cohort (n=40)",
        pad=26,
    )
    # Load-bearing statistics kept terse, just under the title.
    ax.text(
        0.5,
        1.045,
        (
            f"bias-corrected Cramer's V (group x cohort) = "
            f"{v_stats['group_x_cohort_V']:.3f}   |   "
            f"ADD-vs-other = {v_stats['ADD_vs_other_cohort_V']:.3f}   |   "
            f"PD-vs-other = {v_stats['PD_vs_other_cohort_V']:.3f}"
        ),
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=9.5,
        color="#333333",
    )

    # --- companion colorbar legend, rendered as its own figure ---
    legend_fig, legend_ax = plt.subplots(figsize=(4.2, 0.9))
    mappable = ScalarMappable(norm=norm, cmap=cmap)
    cbar = legend_fig.colorbar(
        mappable, cax=legend_ax, orientation="horizontal"
    )
    cbar.set_ticks(range(0, vmax + 1))
    cbar.set_label("Sample count (viridis)", fontsize=10)
    cbar.ax.tick_params(labelsize=9)

    return fig, legend_fig


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(name)s: %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--confounds", type=Path, default=DEFAULT_CONFOUNDS)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--base-name", type=str, default=DEFAULT_BASE_NAME)
    args = parser.parse_args()

    dv = data_version(args.data)
    LOGGER.info("data_version=%s", dv)

    ct = build_crosstab(args.data)
    LOGGER.info("crosstab:\n%s", ct.to_string())
    v_stats = load_cramers_v(args.confounds)

    with publication_style():
        fig, legend_fig = make_figure(ct, v_stats, args.registry)

    artifacts = save_figure(
        fig, args.output_dir, args.base_name, legend_fig=legend_fig
    )
    LOGGER.info("wrote %s", artifacts.svg)
    LOGGER.info("wrote %s", artifacts.png)
    LOGGER.info("wrote %s", artifacts.legend_svg)
    LOGGER.info("wrote %s", artifacts.legend_png)


if __name__ == "__main__":
    main()
