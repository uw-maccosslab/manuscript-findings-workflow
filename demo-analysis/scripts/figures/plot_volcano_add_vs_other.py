"""Volcano plot for the PRIMARY differential-abundance result: ADD vs pooled Other.

Renders the CACHED moderated-model contrast table
(``results/differential_abundance/add-vs-other-moderated/contrast_table.csv``, 2340
proteins) produced by ``scripts/scratch/de_add_vs_other.py`` — a *load -> plot* figure
job; it recomputes NOTHING (conventions/results-cache.md). x = effect
(log2FC, ADD - Other); y = -log10(q) (BH-FDR).

Encoding (colors read from ``state/color_registry.json``; registry read, not mutated):
  * significant hits (q < 0.05) colored by DIRECTION via the ``primary_contrast``
    category -- higher in ADD (effect > 0) = ADD vermillion #D55E00; higher in Other
    (effect < 0) = Other gray #999999. Significance is ALSO encoded by position (above
    the q=0.05 line), so color adds only direction.
  * everything else = neutral background gray (no palette slot).
  * the 11 known plasma / EV AD markers are overlaid in a DISTINCT style so the reader
    sees they are NON-significant: color + shape by the three ``highlight_group``
    registry categories (same assignment as the dynamic-range AD-markers figure), black
    edge, drawn on top, directly labeled.
  * the strongest / most-significant hits are directly labeled (textalloc, collision
    free).

LOAD-BEARING HONESTY (title + caption, findings 0001/0004): EXPLORATORY. ADD is 100%
Kerchner cohort, so the effect is inseparable from a site effect and no cohort covariate
is possible (single-cohort). The delivered matrix is already batch-corrected on cohort,
so the q-values are ANTICONSERVATIVE (variance deflation; p-value null tail suppressed
to ~0.53x uniform) and the 363 q<0.05 hits are an INFLATED count -- read effects + CIs,
not q. Known plasma/EV AD markers (CD9/FLOT1/TSG101/ALIX/syntenin/GFAP/APP/L1CAM/NCAM1/
F13A1/F13B) are NOT significant; the top hits are intracellular proteins.

Run:  ./.venv/bin/python scripts/figures/plot_volcano_add_vs_other.py
"""

from __future__ import annotations

import argparse
import logging
import math
import sys
from pathlib import Path

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import textalloc as ta
from matplotlib.figure import Figure

# Reuse the project shared figure-IO + color-registry modules.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.colors import BACKGROUND_COLOR, get_color
from common.figure_io import publication_style, save_figure

LOGGER = logging.getLogger("plot_volcano_add_vs_other")

DEFAULT_TABLE = (
    Path("results") / "differential_abundance" / "add-vs-other-moderated"
    / "contrast_table.csv"
)
DEFAULT_OUTPUT_DIR = (
    Path("figures") / "analysis" / "differential-abundance" / "add-vs-other"
)
DEFAULT_BASE_NAME = "volcano-add-vs-other-moderated-q05"
DEFAULT_REGISTRY = Path("state") / "color_registry.json"

# Thresholds.
Q_PRIMARY = 0.05
Q_GUIDE = 0.10
LFC_GUIDE = 1.0

# --- Top hits to label (feature id -> on-canvas gene label). ------------------
TOP_HIT_LABELS: dict[str, str] = {
    "sp|P36507|MP2K2_HUMAN": "MEK2",
    "sp|P51812|KS6A3_HUMAN": "RSK2",
    "sp|Q16186|ADRM1_HUMAN": "ADRM1",
    "sp|O60763|USO1_HUMAN": "USO1",
    "sp|Q9P253|VPS18_HUMAN": "VPS18",
    "sp|P29350|PTN6_HUMAN": "SHP1",
    "sp|P19404|NDUV2_HUMAN": "NDUFV2",
    "sp|Q9Y262|EIF3L_HUMAN": "EIF3L",
    "sp|P26440|IVD_HUMAN": "IVD",
}

# --- Known EV / AD markers -> the three highlight_group registry categories. ---
# (Same buckets + colors as scripts/figures/plot_dynamic_range_ad_markers.py.)
G_EV = "EV-QC anchor (MISEV2018)"
G_PARENT = "AD-marker parent protein"
G_CARGO = "EV-cargo / candidate (some contested)"

# marker_class value (as written in the contrast table) -> highlight_group category.
MARKER_CLASS_TO_GROUP: dict[str, str] = {
    "EV-QC (MISEV2018 cytosolic)": G_EV,
    "EV-QC (MISEV2018 tetraspanin)": G_EV,
    "EV-QC (MISEV2018 flotillin)": G_EV,
    "AD glial marker (parent protein)": G_PARENT,
    "AD amyloid parent protein": G_PARENT,
    "contested NDEV capture marker": G_CARGO,
    "unbiased plasma-EV AD candidate (Nielsen 2021, unreplicated)": G_CARGO,
}
# Second channel: a distinct marker shape per highlight group (matches dynamic-range).
GROUP_MARKERS: dict[str, str] = {G_EV: "o", G_PARENT: "s", G_CARGO: "^"}
GROUP_ORDER = [G_EV, G_PARENT, G_CARGO]

# primary_contrast direction labels.
DIR_ADD = "higher in ADD"
DIR_OTHER = "higher in Other"


def _neglog10(q: pd.Series) -> pd.Series:
    """-log10(q) with a tiny floor so an exact zero can't become +inf."""
    floor = np.nextafter(0.0, 1.0)
    return -np.log10(np.clip(q.to_numpy(dtype=float), floor, 1.0))


def build_volcano(
    table: pd.DataFrame, registry_path: Path
) -> tuple[Figure, Figure, dict[str, str]]:
    """Build the volcano figure + its separate legend figure. Registry read-only."""
    df = table.copy()
    if df["q"].isna().any():
        # Constant/untestable features have NaN q (excluded from BH); drop for plotting.
        n_drop = int(df["q"].isna().sum())
        LOGGER.info("dropping %d features with NaN q (untestable)", n_drop)
        df = df[df["q"].notna()].copy()
    df["neglog10q"] = _neglog10(df["q"])

    # Colors from the registry (persist=False -> read, never mutate the shared file).
    c_add = get_color("primary_contrast", "ADD", registry_path=registry_path,
                      persist=False)
    c_other = get_color("primary_contrast", "Other", registry_path=registry_path,
                        persist=False)
    group_colors = {
        g: get_color("highlight_group", g, registry_path=registry_path, persist=False)
        for g in GROUP_ORDER
    }
    color_map: dict[str, str] = {
        DIR_ADD: c_add, DIR_OTHER: c_other, "not significant": BACKGROUND_COLOR,
        **group_colors,
    }

    is_marker = df["marker_class"].notna()
    markers = df[is_marker].copy()
    rest = df[~is_marker].copy()
    sig = rest["q"] < Q_PRIMARY
    up = sig & (rest["effect"] > 0)
    down = sig & (rest["effect"] < 0)
    nonsig = ~sig

    with publication_style():
        fig, ax = plt.subplots(figsize=(7.6, 7.2))

        # Background (non-significant) cloud.
        ax.scatter(rest.loc[nonsig, "effect"], rest.loc[nonsig, "neglog10q"],
                   s=10, c=BACKGROUND_COLOR, alpha=0.45, linewidths=0, zorder=1)
        # Significant hits, colored by direction (position already encodes sig.).
        ax.scatter(rest.loc[down, "effect"], rest.loc[down, "neglog10q"],
                   s=16, c=c_other, alpha=0.85, linewidths=0, zorder=2)
        ax.scatter(rest.loc[up, "effect"], rest.loc[up, "neglog10q"],
                   s=16, c=c_add, alpha=0.85, linewidths=0, zorder=2)

        # Threshold + guide lines.
        y05 = -math.log10(Q_PRIMARY)
        y10 = -math.log10(Q_GUIDE)
        ax.axhline(y05, color="#444444", lw=1.2, ls="-", zorder=3)
        ax.axhline(y10, color="#888888", lw=1.0, ls="--", zorder=3)
        ax.axvline(LFC_GUIDE, color="#bbbbbb", lw=0.9, ls=":", zorder=1)
        ax.axvline(-LFC_GUIDE, color="#bbbbbb", lw=0.9, ls=":", zorder=1)
        ax.axvline(0.0, color="#dddddd", lw=0.8, ls="-", zorder=0)

        xmax = float(np.nanmax(np.abs(df["effect"]))) * 1.08
        ax.set_xlim(-xmax, xmax)
        # Headroom above the highest point so the top-row hit labels have clear space.
        ymin_auto, _ = ax.get_ylim()
        ax.set_ylim(ymin_auto, float(df["neglog10q"].max()) + 0.45)
        ax.text(xmax, y05, "  q=0.05", va="center", ha="left", fontsize=8.5,
                color="#444444")
        ax.text(xmax, y10, "  q=0.10", va="center", ha="left", fontsize=8.5,
                color="#888888")
        ymin = ax.get_ylim()[0]
        ax.text(LFC_GUIDE, ymin, "|log2FC|=1", va="bottom", ha="center",
                fontsize=8, color="#999999")
        ax.text(-LFC_GUIDE, ymin, "|log2FC|=1", va="bottom", ha="center",
                fontsize=8, color="#999999")

        # Marker overlay: color + shape by highlight group, black edge, on top.
        markers = markers.assign(
            group=markers["marker_class"].map(MARKER_CLASS_TO_GROUP)
        )
        if markers["group"].isna().any():
            missing = sorted(
                markers.loc[markers["group"].isna(), "marker_class"].unique()
            )
            raise ValueError(f"Unmapped marker_class values: {missing}")
        for g in GROUP_ORDER:
            sub = markers[markers["group"] == g]
            ax.scatter(sub["effect"], sub["neglog10q"], s=95,
                       marker=GROUP_MARKERS[g], facecolor=group_colors[g],
                       edgecolor="black", linewidths=0.8, alpha=0.95, zorder=5)

        # Direct labels: top hits (bold black) + every marker (its gene, group color).
        hit_rows = df.set_index("feature").loc[list(TOP_HIT_LABELS)]
        lx = list(hit_rows["effect"]) + list(markers["effect"])
        ly = list(hit_rows["neglog10q"]) + list(markers["neglog10q"])
        marker_syms = list(markers["gene"].str.split("/").str[0])
        ltext = list(TOP_HIT_LABELS.values()) + marker_syms
        lcolor = ["black"] * len(TOP_HIT_LABELS) + [
            group_colors[g] for g in markers["group"]
        ]
        lsize = [9.5] * len(TOP_HIT_LABELS) + [8.0] * len(markers)
        lweight = ["bold"] * len(TOP_HIT_LABELS) + ["normal"] * len(markers)
        ta.allocate(
            ax, lx, ly, ltext,
            x_scatter=lx, y_scatter=ly,
            textsize=lsize, textcolor=lcolor, fontweight=lweight,
            linecolor="#999999", linewidth=0.6, draw_all=True,
            min_distance=0.015, max_distance=0.28, margin=0.01,
            avoid_label_lines_overlap=True,
            xlims=ax.get_xlim(), ylims=ax.get_ylim(),
        )

        ax.set_xlabel("log2 fold change (ADD − Other)")  # noqa: RUF001
        ax.set_ylabel(r"$-\log_{10}$(q)  (BH-FDR)")

        n = len(df)
        n_q05 = int((df["q"] < Q_PRIMARY).sum())
        n_q10 = int((df["q"] < Q_GUIDE).sum())
        title = (
            "Volcano — ADD vs Other (pooled PDD/PDCN/HCN), plasma-EV proteome"
            "  ·  moderated (limma eBayes)\n"
            "EXPLORATORY · single-cohort confound: ADD = 100% Kerchner, effect "
            "inseparable from site (findings 0001/0004)\n"
            "q anticonservative (matrix already batch-corrected) → read effect + "
            "CI, not q; hit counts inflated"
        )
        ax.set_title(title, fontsize=10.5, fontweight="bold", loc="left")

        # Terse load-bearing count box (N + hit counts + inflated caveat marker).
        ax.text(
            0.985, 0.985,
            f"n = {n} proteins  (ADD n=10 vs Other n=30)\n"
            f"q<0.05: {n_q05}   q<0.10: {n_q10}   (counts INFLATED — see caption)",
            transform=ax.transAxes, va="top", ha="right", fontsize=8.5,
            color="#333333",
            bbox={"boxstyle": "round,pad=0.35", "fc": "white", "ec": "#cccccc",
                  "alpha": 0.85},
        )

    # --- Separate legend figure (swatch key; never baked onto the plot). ---
    legend_fig = _build_legend(c_add, c_other, group_colors)
    return fig, legend_fig, color_map


def _build_legend(
    c_add: str, c_other: str, group_colors: dict[str, str]
) -> Figure:
    """Standalone swatch/line key for the volcano, as its own figure image."""
    with publication_style():
        lf, la = plt.subplots(figsize=(4.6, 3.0))
        la.axis("off")
        h = [
            mlines.Line2D([], [], marker="o", ls="", ms=8, mfc=c_add, mec="none",
                          label="q<0.05, higher in ADD"),
            mlines.Line2D([], [], marker="o", ls="", ms=8, mfc=c_other, mec="none",
                          label="q<0.05, higher in Other"),
            mlines.Line2D([], [], marker="o", ls="", ms=7, mfc=BACKGROUND_COLOR,
                          mec="none", label="not significant (q≥0.05)"),
            mlines.Line2D([], [], marker=GROUP_MARKERS[G_EV], ls="", ms=9,
                          mfc=group_colors[G_EV], mec="black", mew=0.8,
                          label="EV-QC anchor (MISEV2018)"),
            mlines.Line2D([], [], marker=GROUP_MARKERS[G_PARENT], ls="", ms=9,
                          mfc=group_colors[G_PARENT], mec="black", mew=0.8,
                          label="AD-marker parent protein"),
            mlines.Line2D([], [], marker=GROUP_MARKERS[G_CARGO], ls="", ms=9,
                          mfc=group_colors[G_CARGO], mec="black", mew=0.8,
                          label="EV-cargo / candidate (some contested)"),
            mlines.Line2D([], [], color="#444444", lw=1.2, ls="-", label="q = 0.05"),
            mlines.Line2D([], [], color="#888888", lw=1.0, ls="--", label="q = 0.10"),
            mlines.Line2D([], [], color="#bbbbbb", lw=0.9, ls=":",
                          label="|log2FC| = 1"),
        ]
        la.legend(handles=h, loc="center", frameon=False, fontsize=10,
                  handletextpad=0.6, labelspacing=0.8, title="volcano key",
                  title_fontsize=11)
    return lf


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(name)s: %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--base-name", type=str, default=DEFAULT_BASE_NAME)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    args = parser.parse_args()

    table = pd.read_csv(args.table)
    LOGGER.info("loaded %s: %d rows", args.table, len(table))

    fig, legend_fig, color_map = build_volcano(table, args.registry)
    artifacts = save_figure(
        fig, args.output_dir, args.base_name, legend_fig=legend_fig, dpi=300
    )
    LOGGER.info("color map: %s", color_map)
    for p in (artifacts.svg, artifacts.png, artifacts.legend_svg, artifacts.legend_png):
        LOGGER.info("wrote %s", p)


if __name__ == "__main__":
    main()
