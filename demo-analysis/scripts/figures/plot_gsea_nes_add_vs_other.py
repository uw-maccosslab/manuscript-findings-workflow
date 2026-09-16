"""NES bar chart of the significant preranked-GSEA gene sets: ADD vs pooled Other.

Renders the CACHED preranked-GSEA table
(``results/gsea/fed7a349c550/gsea_all_libraries.csv`` — gseapy prerank across
GO:BP / GO:MF / GO:CC + KEGG, ranked by moderated t) produced upstream. This is a
*load -> subset -> plot* figure job; it recomputes NOTHING (conventions/results-cache.md).

Headline the figure supports: abundance-weighted GSEA finds rich COORDINATED pathway
signal (39 gene sets at FDR<0.05) where a count-based ORA found nothing.

What is drawn:
  * one horizontal bar per gene set significant at FDR<0.05, bar length = NES, sorted by
    NES (positive on top). With 39 significant sets, the panel is capped at the TOP 24 by
    |NES| to stay legible (stated in the caption); the 15 weaker-|NES| sets are omitted.
  * DIRECTION encoded by color via the ``primary_contrast`` registry category (same
    assignment as the volcano): positive NES = up-in-ADD = ADD vermillion #D55E00;
    negative NES = up-in-Other = Other gray #999999.
  * each bar is labeled with a short LIBRARY tag ([GO:BP] / [GO:MF] / [GO:CC] / [KEGG])
    plus the (trimmed) gene-set Term, and annotated with its FDR q at the bar tip.
  * a vertical reference line at NES = 0.

LOAD-BEARING HONESTY (title + caption, findings 0001/0004/0006): EXPLORATORY,
confound-bounded. ADD = 100% Kerchner cohort on an extrapolated batch-corrected matrix,
so enrichment is inseparable from a collection-SITE effect; the up-in-ADD complement /
coagulation theme is a known plasma collection / pre-analytical artifact risk; FDR q's
are anticonservative.

Run:  ./.venv/bin/python scripts/figures/plot_gsea_nes_add_vs_other.py
"""

from __future__ import annotations

import argparse
import logging
import sys
import textwrap
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.figure import Figure

# Reuse the project shared figure-IO + color-registry modules.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.colors import get_color
from common.figure_io import publication_style, save_figure

LOGGER = logging.getLogger("plot_gsea_nes_add_vs_other")

DEFAULT_TABLE = Path("results") / "gsea" / "fed7a349c550" / "gsea_all_libraries.csv"
DEFAULT_OUTPUT_DIR = Path("figures") / "analysis" / "gsea" / "add-vs-other"
DEFAULT_BASE_NAME = "gsea-nes-fdr05-add-vs-other"
DEFAULT_REGISTRY = Path("state") / "color_registry.json"

FDR_THRESHOLD = 0.05
TOP_N = 24  # cap for legibility when > this many sets are significant

# gseapy library name -> short on-canvas tag.
LIBRARY_TAG: dict[str, str] = {
    "GO_Biological_Process_2021": "GO:BP",
    "GO_Molecular_Function_2021": "GO:MF",
    "GO_Cellular_Component_2021": "GO:CC",
    "KEGG_2021_Human": "KEGG",
}

DIR_ADD = "up in ADD (NES > 0)"
DIR_OTHER = "up in Other (NES < 0)"


def _clean_term(term: str) -> str:
    """Drop the trailing GO accession and squeeze whitespace for a readable label."""
    # e.g. "regulation of complement activation (GO:0030449)" -> "regulation of ..."
    if "(GO:" in term:
        term = term[: term.rindex("(GO:")].strip()
    return " ".join(term.split())


def _fmt_q(q: float) -> str:
    """Terse FDR annotation."""
    if q < 0.001:
        return "q<0.001"
    return f"q={q:.3f}"


def build_figure(
    table: pd.DataFrame, registry_path: Path
) -> tuple[Figure, Figure, dict[str, str], int, int]:
    """Build the NES bar chart + its separate legend figure. Registry read-only."""
    sig = table[table["FDR q-val"] < FDR_THRESHOLD].copy()
    n_sig = len(sig)
    if n_sig == 0:
        raise ValueError("No gene sets at FDR<0.05 — nothing to plot.")

    sig["abs_nes"] = sig["NES"].abs()
    capped = n_sig > TOP_N
    if capped:
        sig = sig.sort_values("abs_nes", ascending=False).head(TOP_N)
    n_shown = len(sig)

    # Sort for display: largest positive NES on top, largest negative at bottom.
    sig = sig.sort_values("NES", ascending=True).reset_index(drop=True)

    c_add = get_color("primary_contrast", "ADD", registry_path=registry_path,
                      persist=False)
    c_other = get_color("primary_contrast", "Other", registry_path=registry_path,
                        persist=False)
    color_map = {DIR_ADD: c_add, DIR_OTHER: c_other}

    unknown = set(sig["Library"]) - set(LIBRARY_TAG)
    if unknown:
        raise ValueError(f"Unmapped Library values: {sorted(unknown)}")

    with publication_style():
        fig, ax = plt.subplots(figsize=(9.2, 0.42 * n_shown + 1.9))

        y = range(n_shown)
        colors = [c_add if nes > 0 else c_other for nes in sig["NES"]]
        ax.barh(list(y), sig["NES"], color=colors, edgecolor="black",
                linewidth=0.4, height=0.72, zorder=2)

        ax.axvline(0.0, color="#333333", lw=1.1, zorder=3)

        nes_max = float(sig["NES"].max())
        nes_min = float(sig["NES"].min())
        span = nes_max - nes_min
        pad = span * 0.02
        # Term labels sit at the zero end, on the side opposite the bar direction.
        for i, row in sig.iterrows():
            tag = LIBRARY_TAG[row["Library"]]
            term = _clean_term(row["Term"])
            term = textwrap.shorten(term, width=52, placeholder="…")
            label = f"[{tag}] {term}"
            if row["NES"] > 0:
                ax.text(-pad, i, label, ha="right", va="center", fontsize=8.2)
                ax.text(row["NES"] + pad, i, _fmt_q(row["FDR q-val"]),
                        ha="left", va="center", fontsize=7.2, color="#555555")
            else:
                ax.text(pad, i, label, ha="left", va="center", fontsize=8.2)
                ax.text(row["NES"] - pad, i, _fmt_q(row["FDR q-val"]),
                        ha="right", va="center", fontsize=7.2, color="#555555")

        ax.set_yticks([])
        ax.set_ylim(-0.7, n_shown - 0.3)
        # Room on both sides for the term labels + q annotations.
        ax.set_xlim(nes_min - span * 0.55, nes_max + span * 0.30)
        ax.spines["left"].set_visible(False)
        ax.set_xlabel("Normalized enrichment score (NES)")

        title = (
            "Preranked GSEA — ADD vs Other (pooled PDD/PDCN/HCN), plasma-EV proteome\n"
            f"gene sets at FDR<0.05 ({n_sig} total"
            + (f"; top {n_shown} by |NES| shown" if capped else "")
            + ")  ·  bar = NES, colored by direction\n"
            "EXPLORATORY · confound-bounded: ADD = 100% Kerchner, enrichment "
            "inseparable from collection SITE; FDR anticonservative (findings "
            "0001/0004/0006)"
        )
        ax.set_title(title, fontsize=10.0, fontweight="bold", loc="left")

    legend_fig = _build_legend(c_add, c_other)
    return fig, legend_fig, color_map, n_sig, n_shown


def _build_legend(c_add: str, c_other: str) -> Figure:
    """Standalone key: direction colors + the library-tag glossary."""
    with publication_style():
        lf, la = plt.subplots(figsize=(4.4, 2.4))
        la.axis("off")
        handles = [
            mpatches.Patch(facecolor=c_add, edgecolor="black", linewidth=0.4,
                           label=DIR_ADD),
            mpatches.Patch(facecolor=c_other, edgecolor="black", linewidth=0.4,
                           label=DIR_OTHER),
        ]
        leg1 = la.legend(handles=handles, loc="upper left", frameon=False,
                         fontsize=10, title="direction (NES sign)", title_fontsize=11,
                         bbox_to_anchor=(0.0, 1.0))
        la.add_artist(leg1)
        tag_text = (
            "library tags\n"
            "[GO:BP]  GO Biological Process\n"
            "[GO:MF]  GO Molecular Function\n"
            "[GO:CC]  GO Cellular Component\n"
            "[KEGG]  KEGG pathways"
        )
        la.text(0.0, 0.52, tag_text, va="top", ha="left", fontsize=9.0)
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

    fig, legend_fig, color_map, n_sig, n_shown = build_figure(table, args.registry)
    artifacts = save_figure(
        fig, args.output_dir, args.base_name, legend_fig=legend_fig, dpi=300
    )
    LOGGER.info("FDR<0.05: %d sets; shown: %d", n_sig, n_shown)
    LOGGER.info("color map: %s", color_map)
    for p in (artifacts.svg, artifacts.png, artifacts.legend_svg, artifacts.legend_png):
        LOGGER.info("wrote %s", p)


if __name__ == "__main__":
    main()
