"""Four-arm profile of the up-in-ADD complement/coagulation program — cohort-colored.

Renders the CACHED four-arm complement result
(``results/four_arm_complement/dd884570c4fc/``, fingerprint ``dd884570c4fc``, data
version ``protein-data.csv sha256:90a0462eb03a02bd``) — a pure *load -> plot* figure job;
it recomputes NOTHING (conventions/results-cache.md). Everything is read straight from
the cache parquets:

  * ``sample_scores.parquet``  — per-sample summary z-score (``complement_score``, the
    KEGG-27 primary set), each sample's ``arm`` and collection ``cohort``  (n=10/arm).
  * ``per_arm_scores.parquet`` — per-arm ``mean_z`` + 95% CI (``ci_low``/``ci_high``),
    filtered to the primary set ``KEGG_complement_coag``.
  * ``arm_comparisons.parquet``— ADD-vs-each-arm Hedges g + 95% CI + BH-q (annotation),
    primary set ``KEGG_complement_coag``.

The read: the up-in-ADD KEGG "Complement and coagulation cascades" leading-edge program
(27 proteins, from finding 0011) is elevated SPECIFICALLY in the ADD arm (+0.53), while
PDD (-0.15), PDCN (-0.31) and HCN (-0.08) sit at/below the z=0 baseline. Each per-sample
point is COLORED BY ITS COLLECTION COHORT (Kerchner / Poston / ADRC/PUC / PUC only) —
this is LOAD-BEARING: the majority-Kerchner HCN arm sits at baseline, which weakly argues
against a pure Kerchner collection artifact. But the anti-artifact cue is soft:
ADD-vs-HCN is NOT BH-significant (g=0.80, q=0.081), so the annotation reads the q
honestly and marks it "n.s.". ADD-vs-PDD (g=1.33, q=0.009) and ADD-vs-PDCN (g=1.59,
q=0.005) are significant.

Confound bound (EXPLORATORY; state/METADATA.md, findings 0001/0004): the leading-edge is
defined from the ADD-vs-Other ranking, so ADD is extreme BY CONSTRUCTION — the
non-circular comparison is among PDD/PDCN/HCN; and ADD is single-cohort (100% Kerchner),
so ADD-disease vs ADD-Kerchner is not separable. n=10/arm (finding 0003).

Encoding (colors read from ``state/color_registry.json`` under the ``cohort`` category;
registry read, never mutated): point fill = collection cohort (Kerchner #56B4E9, Poston
#CC79A7, ADRC/PUC #F0E442, PUC only #000000). The per-arm mean is a BLACK diamond with a
white edge (a distinct shape, so it is never confused with a black "PUC only" sample
point) over a black 95% CI whisker; a dashed line marks z = 0. Four cohort colors, well
within the 8-color budget.

Run:  ./.venv/bin/python scripts/figures/plot_fourarm_complement_coag_per_arm_scores.py
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

# Reuse the project shared figure-IO + color-registry modules.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.colors import get_color
from common.figure_io import publication_style, save_figure

LOGGER = logging.getLogger("plot_fourarm_complement_coag_per_arm_scores")

DEFAULT_CACHE_DIR = Path("results") / "four_arm_complement" / "dd884570c4fc"
DEFAULT_OUTPUT_DIR = (
    Path("figures") / "analysis" / "four_arm_complement" / "complement-coagulation"
)
DEFAULT_BASE_NAME = "fourarm-complement-coag-per-arm-scores"
DEFAULT_REGISTRY = Path("state") / "color_registry.json"

# Fixed display order of the arms (ADD first) and the primary (KEGG-27) set.
ARM_ORDER = ["ADD", "PDD", "PDCN", "HCN"]
PRIMARY_SET = "KEGG_complement_coag"
SAMPLE_COL = "complement_score"

# Fixed legend order for the cohorts (only those actually present are shown).
COHORT_ORDER = ["Kerchner", "Poston", "ADRC/PUC", "PUC only"]

# BH significance threshold — only for reading the q honestly (n.s. marking).
Q_ALPHA = 0.05

_RNG_SEED = 20260915  # deterministic strip jitter (reproducible render)


def _cohort_colors(cohorts: list[str], registry_path: Path) -> dict[str, str]:
    """Hex color per cohort from the ``cohort`` registry category; read-only."""
    return {
        c: get_color("cohort", c, registry_path=registry_path, persist=False)
        for c in cohorts
    }


def build_figure(
    sample_scores: pd.DataFrame,
    per_arm: pd.DataFrame,
    comparisons: pd.DataFrame,
    registry_path: Path,
) -> tuple[Figure, Figure, dict[str, str]]:
    """Build the single-panel four-arm figure + its separate legend. Registry read-only."""
    for arm in ARM_ORDER:
        if arm not in set(sample_scores["arm"]):
            raise ValueError(f"arm {arm!r} missing from sample_scores.")

    pset = per_arm[per_arm["set"] == PRIMARY_SET].set_index("arm")
    if pset.empty:
        raise ValueError(f"primary set {PRIMARY_SET!r} missing from per_arm_scores.")
    cset = comparisons[comparisons["set"] == PRIMARY_SET].set_index("comparison")

    # Cohorts present, in fixed legend order (fail loud on an unknown cohort).
    present = set(sample_scores["cohort"])
    unknown = present - set(COHORT_ORDER)
    if unknown:
        raise ValueError(f"unexpected cohort(s) not in COHORT_ORDER: {sorted(unknown)}")
    cohorts = [c for c in COHORT_ORDER if c in present]
    colors = _cohort_colors(cohorts, registry_path)

    rng = np.random.default_rng(_RNG_SEED)

    ymin = float(sample_scores[SAMPLE_COL].min())
    ymax = float(sample_scores[SAMPLE_COL].max())
    pad = 0.14 * (ymax - ymin)
    ylo, yhi = ymin - pad, ymax + pad * 2.2  # head-room for the g/q annotations

    with publication_style():
        fig, ax = plt.subplots(figsize=(7.6, 5.6))

        ax.axhline(0.0, color="black", lw=0.9, ls="--", alpha=0.6, zorder=1)

        for x, arm in enumerate(ARM_ORDER):
            sub = sample_scores[sample_scores["arm"] == arm]
            pts = sub[SAMPLE_COL].to_numpy(float)
            pt_cols = [colors[c] for c in sub["cohort"]]
            jitter = rng.uniform(-0.16, 0.16, size=pts.size)
            ax.scatter(
                np.full(pts.size, x) + jitter, pts, s=42, c=pt_cols,
                alpha=0.9, edgecolor="#555555", linewidths=0.4, zorder=3,
            )

            m = float(pset.loc[arm, "mean_z"])
            lo = float(pset.loc[arm, "ci_low"])
            hi = float(pset.loc[arm, "ci_high"])
            ax.plot([x, x], [lo, hi], color="black", lw=1.8, zorder=5)
            ax.plot([x - 0.09, x + 0.09], [lo, lo], color="black", lw=1.8, zorder=5)
            ax.plot([x - 0.09, x + 0.09], [hi, hi], color="black", lw=1.8, zorder=5)
            ax.scatter(
                [x], [m], marker="D", s=120, color="black", edgecolor="white",
                linewidths=1.3, zorder=6,
            )
            ax.annotate(
                f"{m:+.2f}", xy=(x + 0.20, m), va="center", ha="left",
                fontsize=8.5, color="#222222", zorder=7,
            )

            # ADD-vs-arm effect (Hedges g + BH-q) atop the non-ADD arms; q read honestly.
            if arm != "ADD":
                key = f"ADD_vs_{arm}"
                if key in cset.index:
                    g = float(cset.loc[key, "hedges_g"])
                    q = float(cset.loc[key, "q_bh"])
                    sig = "" if q < Q_ALPHA else "  n.s."
                    ax.annotate(
                        f"g = {g:.2f}\nq = {q:.3f}{sig}",
                        xy=(x, yhi - 0.4 * pad), va="top", ha="center",
                        fontsize=7.5, color="#333333", zorder=7,
                    )

        ax.set_xticks(range(len(ARM_ORDER)))
        ax.set_xticklabels([f"{arm}\n(n=10)" for arm in ARM_ORDER], fontsize=10)
        ax.set_xlim(-0.6, len(ARM_ORDER) - 0.4)
        ax.set_ylim(ylo, yhi)
        ax.set_ylabel("summary z-score  (mean per-protein z, 27 proteins)")
        ax.set_title(
            "Complement/coagulation leading-edge (KEGG, 27 proteins) — four-arm profile",
            fontsize=12, fontweight="bold",
        )
        fig.tight_layout(rect=(0, 0.045, 1, 1))
        # Terse mandatory caveat marker on-canvas (bottom footnote, clear of the title
        # and the data); the full detail lives in the finding caption.
        fig.text(
            0.5, 0.012,
            "EXPLORATORY · confound-bounded: leading-edge defined from ADD-vs-Other "
            "(ADD extreme by construction); ADD is 100% Kerchner",
            ha="center", va="bottom", fontsize=7.5, color="#666666",
        )

    legend_fig = _build_legend(colors, cohorts)
    return fig, legend_fig, colors


def _build_legend(colors: dict[str, str], cohorts: list[str]) -> Figure:
    """Standalone key: cohort swatches + the glyph meanings, as its own figure image."""
    with publication_style():
        lf, la = plt.subplots(figsize=(5.2, 3.0))
        la.axis("off")
        handles = [
            mlines.Line2D([], [], marker="o", ls="", ms=10, mfc=colors[c],
                          mec="#555555", mew=0.6, label=f"cohort: {c}")
            for c in cohorts
        ]
        handles += [
            mlines.Line2D([], [], marker="D", ls="", ms=11, mfc="black", mec="white",
                          mew=1.3, label="per-arm mean z (± black 95% CI)"),
            mlines.Line2D([], [], marker="o", ls="", ms=8, mfc="#999999", mec="#555555",
                          mew=0.6, label="one sample (n=10/arm)"),
            mlines.Line2D([], [], color="black", lw=0.9, ls="--", alpha=0.6,
                          label="z = 0 reference"),
        ]
        la.legend(handles=handles, loc="center", frameon=False, fontsize=10,
                  handletextpad=0.7, labelspacing=0.9,
                  title="four-arm complement/coagulation key", title_fontsize=11)
    return lf


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(name)s: %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--base-name", type=str, default=DEFAULT_BASE_NAME)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    args = parser.parse_args()

    sample_scores = pd.read_parquet(args.cache_dir / "sample_scores.parquet")
    per_arm = pd.read_parquet(args.cache_dir / "per_arm_scores.parquet")
    comparisons = pd.read_parquet(args.cache_dir / "arm_comparisons.parquet")
    LOGGER.info(
        "loaded four_arm_complement cache %s: %d samples, %d per-arm rows, %d comparisons",
        args.cache_dir, len(sample_scores), len(per_arm), len(comparisons),
    )

    fig, legend_fig, colors = build_figure(
        sample_scores, per_arm, comparisons, args.registry
    )
    artifacts = save_figure(
        fig, args.output_dir, args.base_name, legend_fig=legend_fig, dpi=300
    )
    LOGGER.info("cohort colors: %s", colors)
    for p in (artifacts.svg, artifacts.png, artifacts.legend_svg, artifacts.legend_png):
        LOGGER.info("wrote %s", p)


if __name__ == "__main__":
    main()
