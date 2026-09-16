"""Classic GSEA enrichment plot for the single strongest set: ADD vs pooled Other.

Renders the KEGG_2021_Human set **"Complement and coagulation cascades"**
(NES=+2.50, FDR q≈0, leading edge 27, up-in-ADD) from the CACHED preranked GSEA
(``results/gsea/fed7a349c550/``) produced by ``scripts/scratch/gsea_add_vs_other.py``.

This is a *load + DETERMINISTIC-recompute-of-the-walk -> plot* figure job. It does NOT
re-run the 1000-permutation GSEA (conventions/results-cache.md): the running
enrichment-score curve is a *deterministic* weighted Kolmogorov-Smirnov walk down the
ranked list, so it is recomputed directly from the exact preranked list
(``ranked_symbols.csv``) + the KEGG member set. NES / p / FDR are read verbatim from the
cached KEGG result row (``gsea_kegg.csv``) — never recomputed. The recomputed ES is
asserted against the cached ES so the walk is provably the same one the permutation test
scored.

The member set is resolved by LIBRARY NAME via ``gseapy.get_library("KEGG_2021_Human")``
(the exact library the analysis used; gseapy version pinned in meta.json), so the figure
is regenerable from the library name rather than a machine-specific gmt path.

Standard three-panel GSEA layout (shared x = rank position in the descending list):
  * top    = running enrichment score vs rank; peak marked = ES; zero reference line.
  * middle = gene-hit rug: a tick at each ranked position where a set member falls
             (leading-edge members — those at/before the peak — drawn emphasized).
  * bottom = the ranked metric (signed moderated t) as a filled profile vs rank, with
             the sign-change (metric=0) crossing marked.

Genes are ranked by signed moderated t (up-in-ADD at LEFT). The load-bearing point is
that the members concentrate at the LEFT (positive-metric / up-in-ADD) extreme.

Colors read from ``state/color_registry.json`` (read, never mutated): the enrichment
direction is up-in-ADD, so the ES curve, the hit ticks, and the positive-metric fill use
the ``primary_contrast`` ADD vermillion; the negative-metric fill uses the
``primary_contrast`` Other gray.

LOAD-BEARING HONESTY (title + caption, findings 0001/0004/0006): EXPLORATORY /
confound-bounded. Complement/coagulation elevation in plasma is a known
collection/pre-analytical artifact and ADD is 100% Kerchner cohort, so this is
inseparable from a collection-SITE effect; the FDR is anticonservative (finding 0006).

Run:  ./.venv/bin/python \
        scripts/figures/plot_gsea_enrichment_complement_coagulation_add_vs_other.py
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.patches import Patch

# Reuse the project shared figure-IO + color-registry modules.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.colors import get_color
from common.figure_io import publication_style, save_figure

LOGGER = logging.getLogger("plot_gsea_enrichment_complement_coagulation")

DEFAULT_RANKED = Path("results") / "gsea" / "fed7a349c550" / "ranked_symbols.csv"
DEFAULT_KEGG = Path("results") / "gsea" / "fed7a349c550" / "gsea_kegg.csv"
DEFAULT_OUTPUT_DIR = Path("figures") / "analysis" / "gsea" / "add-vs-other"
DEFAULT_BASE_NAME = "gsea-enrichment-complement-coagulation-add-vs-other"
DEFAULT_REGISTRY = Path("state") / "color_registry.json"

LIBRARY = "KEGG_2021_Human"
TERM = "Complement and coagulation cascades"

# GSEA classic weighting exponent (p=1); the analysis used gseapy prerank defaults.
WEIGHT_P = 1.0
# Tolerance for asserting the recomputed walk reproduces the cached ES.
ES_TOL = 1e-4


def _running_es(
    metric: np.ndarray, hits: np.ndarray
) -> tuple[np.ndarray, int, float]:
    """Weighted-KS running enrichment score down the (descending) ranked list.

    ``metric`` is the ranking statistic sorted descending; ``hits`` is the boolean
    membership mask in the same order. Returns (running ES array, peak index, ES).
    """
    n = metric.size
    n_hit = int(hits.sum())
    if n_hit == 0:
        raise ValueError("no set members present in the ranked list.")
    if n_hit == n:
        raise ValueError("every ranked gene is a member; ES undefined.")
    weights = np.abs(metric) ** WEIGHT_P
    n_r = float(weights[hits].sum())
    if n_r <= 0:
        raise ValueError("total hit weight is zero; cannot normalize the walk.")
    step = np.where(hits, weights / n_r, -1.0 / (n - n_hit))
    running = np.cumsum(step)
    peak = int(np.argmax(np.abs(running)))
    return running, peak, float(running[peak])


def build_figure(
    ranked: pd.DataFrame,
    kegg_row: pd.Series,
    registry_path: Path,
) -> tuple[Figure, Figure, dict[str, str]]:
    """Build the GSEA enrichment figure + its separate legend. Registry read-only."""
    # Exact preranked order: descending by the signed moderated-t ranking metric.
    ranked = ranked.sort_values("statistic", ascending=False).reset_index(drop=True)
    metric = ranked["statistic"].to_numpy(dtype=float)
    symbols = ranked["symbol"].tolist()

    # Member set by LIBRARY NAME (the exact library the analysis used).
    import gseapy as gp

    library = gp.get_library(name=LIBRARY)
    if TERM not in library:
        raise ValueError(f"term {TERM!r} not in library {LIBRARY!r}.")
    member_set = set(library[TERM])
    hits = np.array([s in member_set for s in symbols], dtype=bool)
    n_members_present = int(hits.sum())

    running, peak, es = _running_es(metric, hits)

    # Provably the same walk the permutation test scored.
    cached_es = float(kegg_row["ES"])
    if abs(es - cached_es) > ES_TOL:
        raise ValueError(
            f"recomputed ES {es:.6f} disagrees with cached ES {cached_es:.6f} "
            f"(tol {ES_TOL}); the walk does not reproduce the cached result."
        )

    nes = float(kegg_row["NES"])
    fdr = float(kegg_row["FDR q-val"])
    lead_genes = str(kegg_row["Lead_genes"]).split(";")
    n_lead = len(lead_genes)

    # Colors from the registry (persist=False -> read, never mutate the shared file).
    c_add = get_color("primary_contrast", "ADD", registry_path=registry_path,
                      persist=False)
    c_other = get_color("primary_contrast", "Other", registry_path=registry_path,
                        persist=False)
    color_map = {
        "running ES / hits / up-in-ADD (primary_contrast ADD)": c_add,
        "down-in-ADD metric fill (primary_contrast Other)": c_other,
    }

    x = np.arange(len(symbols))
    hit_x = x[hits]
    # Leading-edge members: hits at or before the peak (they build the max score).
    lead_mask = hits & (x <= peak)
    lead_x = x[lead_mask]
    zero_cross = int(np.searchsorted(-metric, 0.0))  # first index where metric < 0

    with publication_style():
        fig, axes = plt.subplots(
            3, 1, figsize=(8.2, 6.6), sharex=True,
            gridspec_kw={"height_ratios": [3.1, 0.5, 1.4], "hspace": 0.08},
        )
        ax_es, ax_rug, ax_metric = axes

        # --- Top: running enrichment score --------------------------------------
        ax_es.plot(x, running, color=c_add, lw=1.8, zorder=3)
        ax_es.axhline(0.0, color="#888888", lw=0.9, ls="-", zorder=1)
        ax_es.axvline(peak, color="#bbbbbb", lw=0.9, ls=":", zorder=1)
        ax_es.scatter([peak], [es], s=45, color=c_add, edgecolor="black",
                      linewidths=0.8, zorder=4)
        ax_es.annotate(
            f"peak ES = {es:.3f}\n(rank {peak})",
            xy=(peak, es), xytext=(peak + len(symbols) * 0.06, es - 0.10),
            fontsize=9, color="#222222", va="top",
            arrowprops={"arrowstyle": "-", "color": "#999999", "lw": 0.7},
        )
        ax_es.set_ylabel("running enrichment score")
        ax_es.set_ylim(min(-0.05, running.min() - 0.03), es + 0.12)
        ax_es.set_title(
            f"GSEA — KEGG '{TERM}'  ·  ADD vs pooled Other\n"
            f"NES = {nes:+.2f}   FDR {_fmt_q(fdr)}   leading edge {n_lead} "
            f"(of {n_members_present} in ranking)   ·   up-in-ADD at left\n"
            "EXPLORATORY · confounded (see 0001/0004/0006)",
            fontsize=9.5, fontweight="bold", loc="left",
        )

        # --- Middle: gene-hit rug ----------------------------------------------
        ax_rug.vlines(hit_x, 0.0, 1.0, color=c_add, lw=0.6, alpha=0.55, zorder=2)
        ax_rug.vlines(lead_x, 0.0, 1.0, color=c_add, lw=1.1, alpha=1.0, zorder=3)
        ax_rug.set_ylim(0, 1)
        ax_rug.set_yticks([])
        ax_rug.set_ylabel("hits", rotation=0, ha="right", va="center",
                          fontsize=9)
        for spine in ("left", "top", "right"):
            ax_rug.spines[spine].set_visible(False)

        # --- Bottom: ranked metric profile -------------------------------------
        pos = metric >= 0
        ax_metric.fill_between(x, 0.0, metric, where=pos, color=c_add,
                               linewidth=0, alpha=0.75, step="mid")
        ax_metric.fill_between(x, 0.0, metric, where=~pos, color=c_other,
                               linewidth=0, alpha=0.75, step="mid")
        ax_metric.axhline(0.0, color="#888888", lw=0.8, zorder=2)
        ax_metric.axvline(zero_cross, color="#bbbbbb", lw=0.9, ls=":", zorder=1)
        ax_metric.text(
            zero_cross, ax_metric.get_ylim()[1], "  metric = 0", fontsize=8,
            color="#777777", va="top", ha="left",
        )
        ax_metric.set_ylabel("signed moderated t")
        ax_metric.set_xlabel(
            f"gene rank in ADD-vs-Other list (0 = most up-in-ADD; "
            f"n = {len(symbols)} genes)"
        )
        ax_metric.set_xlim(0, len(symbols) - 1)

    legend_fig = _build_legend(c_add, c_other)
    return fig, legend_fig, color_map


def _fmt_q(q: float) -> str:
    """Render an FDR q-value, calling out an exact permutation-floor zero as ~0."""
    if q <= 0.0:
        return "q ≈ 0"
    if q < 1e-3:
        return f"q = {q:.1e}"
    return f"q = {q:.3f}"


def _build_legend(c_add: str, c_other: str) -> Figure:
    """Standalone key for the GSEA plot, as its own figure image."""
    with publication_style():
        lf, la = plt.subplots(figsize=(5.0, 2.6))
        la.axis("off")
        handles = [
            mlines.Line2D([], [], color=c_add, lw=1.8,
                          label="running enrichment score (up-in-ADD)"),
            mlines.Line2D([], [], marker="o", ls="", ms=8, mfc=c_add,
                          mec="black", mew=0.8, label="peak = enrichment score (ES)"),
            mlines.Line2D([], [], color=c_add, lw=1.1,
                          label="leading-edge member (rank ≤ peak)"),
            mlines.Line2D([], [], color=c_add, lw=0.6, alpha=0.55,
                          label="set member (gene-hit tick)"),
            Patch(facecolor=c_add, alpha=0.75,
                  label="metric > 0 (up-in-ADD)"),
            Patch(facecolor=c_other, alpha=0.75,
                  label="metric < 0 (up-in-Other)"),
        ]
        la.legend(handles=handles, loc="center", frameon=False, fontsize=10,
                  handletextpad=0.7, labelspacing=0.8,
                  title="GSEA enrichment key", title_fontsize=11)
    return lf


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(name)s: %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranked", type=Path, default=DEFAULT_RANKED)
    parser.add_argument("--kegg", type=Path, default=DEFAULT_KEGG)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--base-name", type=str, default=DEFAULT_BASE_NAME)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    args = parser.parse_args()

    ranked = pd.read_csv(args.ranked)
    LOGGER.info("loaded ranked list %s: %d genes", args.ranked, len(ranked))
    kegg = pd.read_csv(args.kegg)
    match = kegg[kegg["Term"] == TERM]
    if len(match) != 1:
        raise ValueError(
            f"expected exactly 1 KEGG row for {TERM!r}; found {len(match)}."
        )
    kegg_row = match.iloc[0]

    fig, legend_fig, color_map = build_figure(ranked, kegg_row, args.registry)
    artifacts = save_figure(
        fig, args.output_dir, args.base_name, legend_fig=legend_fig, dpi=300
    )
    LOGGER.info("color map: %s", color_map)
    for p in (artifacts.svg, artifacts.png, artifacts.legend_svg, artifacts.legend_png):
        LOGGER.info("wrote %s", p)


if __name__ == "__main__":
    main()
