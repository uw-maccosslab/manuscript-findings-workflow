"""Four-arm signature-score plot — the ADD-specificity figure.

Renders the CACHED four-arm result (``results/four_arm/8d418bb8509f/``, fingerprint
``8d418bb8509f``, data version ``protein-data.csv sha256:90a0462eb03a02bd``) — a pure
*load -> plot* figure job; it recomputes NOTHING (conventions/results-cache.md). The
per-sample summary scores, per-arm means + 95% CI, and ADD-vs-arm effect sizes are read
straight from the cache parquets:

  * ``sample_scores.parquet``  — per-sample up-score / down-score + arm  (n=10/arm).
  * ``per_arm_scores.parquet`` — per-arm mean_z, ci_low, ci_high per set.
  * ``arm_comparisons.parquet``— ADD-vs-each-arm Hedges g + CI + BH-q (annotation).

The two-panel read: for the up-in-ADD (18 proteins) and down-in-ADD (22 proteins)
Boruta-confirmed signature subsets, the four arms ADD / PDD / PDCN / HCN are shown side
by side as a per-arm mean z-score with 95% CI (dot-and-whisker), overlaid on the 10
per-sample points (jittered strip) so the spread is visible. A horizontal reference line
sits at z=0. The signature is ADD-SPECIFIC: the down-set drops ONLY in ADD (-1.06),
while PDD (dementia control) sits at +0.47, essentially equal to non-demented PDCN
(+0.55) — so PDD != ADD (not dementia-generic). The up-set mirrors it.

CONFOUND made visible on the same axes (LOAD-BEARING; state/METADATA.md Table 1,
findings 0001/0004): each arm's Kerchner-cohort fraction is printed under its label
(ADD 100%, PDD 0%, PDCN 0%, HCN 70%) — the arms track the Kerchner-vs-Poston axis, so
"ADD-specific" is not separable from "Kerchner-site-specific". n=10/arm (finding 0003).
The up/down split is defined by ADD-vs-Other, so the non-circular comparison is among
PDD/PDCN/HCN.

Encoding (colors read from ``state/color_registry.json`` under the ``group`` category;
registry read, never mutated): arm identity -> Okabe-Ito color (ADD #D55E00, PDD #E69F00,
PDCN #009E73, HCN #0072B2). Color reinforces the x-position, which already encodes arm;
four categories, well within the 8-color budget.

Run:  ./.venv/bin/python scripts/figures/plot_four_arm_signature_scores.py
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

LOGGER = logging.getLogger("plot_four_arm_signature_scores")

DEFAULT_CACHE_DIR = Path("results") / "four_arm" / "8d418bb8509f"
DEFAULT_OUTPUT_DIR = Path("figures") / "analysis" / "four_arm" / "add-specificity"
DEFAULT_BASE_NAME = "four-arm-signature-scores-add-specificity"
DEFAULT_REGISTRY = Path("state") / "color_registry.json"

# Fixed display order of the arms (ADD first) and the two signature sets.
ARM_ORDER = ["ADD", "PDD", "PDCN", "HCN"]
SET_ORDER = ["up_in_ADD", "down_in_ADD"]

# n proteins per subset (from signature_direction.parquet: 18 up, 22 down) and each
# arm's Kerchner-cohort fraction (state/METADATA.md Table 1) — the load-bearing confound.
SET_TITLE = {
    "up_in_ADD": "up-in-ADD set  (18 proteins)",
    "down_in_ADD": "down-in-ADD set  (22 proteins)",
}
SET_SAMPLE_COL = {"up_in_ADD": "up_score", "down_in_ADD": "down_score"}
KERCHNER_FRACTION = {"ADD": "100%", "PDD": "0%", "PDCN": "0%", "HCN": "70%"}

_RNG_SEED = 20260915  # deterministic strip jitter (reproducible render)


def _arm_colors(registry_path: Path) -> dict[str, str]:
    """Hex color per arm from the ``group`` registry category; read-only."""
    return {
        arm: get_color("group", arm, registry_path=registry_path, persist=False)
        for arm in ARM_ORDER
    }


def build_four_arm_figure(
    sample_scores: pd.DataFrame,
    per_arm: pd.DataFrame,
    comparisons: pd.DataFrame,
    registry_path: Path,
) -> tuple[Figure, Figure, dict[str, str]]:
    """Build the two-panel four-arm figure + its separate legend. Registry read-only."""
    for arm in ARM_ORDER:
        if arm not in set(sample_scores["arm"]):
            raise ValueError(f"arm {arm!r} missing from sample_scores.")
    colors = _arm_colors(registry_path)
    rng = np.random.default_rng(_RNG_SEED)

    # Shared y-range across both panels so the arms are directly comparable.
    ymin = float(min(sample_scores["up_score"].min(), sample_scores["down_score"].min()))
    ymax = float(max(sample_scores["up_score"].max(), sample_scores["down_score"].max()))
    pad = 0.12 * (ymax - ymin)
    ylo, yhi = ymin - pad, ymax + pad * 1.6  # extra head-room for q-annotations

    with publication_style():
        fig, axes = plt.subplots(1, 2, figsize=(9.6, 5.4), sharey=True)

        for panel, (ax, set_name) in enumerate(zip(axes, SET_ORDER, strict=True)):
            col = SET_SAMPLE_COL[set_name]
            pset = per_arm[per_arm["set"] == set_name].set_index("arm")
            cset = comparisons[comparisons["set"] == set_name].set_index("comparison")

            ax.axhline(0.0, color="black", lw=0.9, ls="--", alpha=0.6, zorder=1)

            for x, arm in enumerate(ARM_ORDER):
                c = colors[arm]
                pts = sample_scores.loc[sample_scores["arm"] == arm, col].to_numpy(float)
                jitter = rng.uniform(-0.14, 0.14, size=pts.size)
                ax.scatter(
                    np.full(pts.size, x) + jitter, pts, s=28, color=c,
                    alpha=0.45, edgecolor="none", zorder=3,
                )
                m = float(pset.loc[arm, "mean_z"])
                lo = float(pset.loc[arm, "ci_low"])
                hi = float(pset.loc[arm, "ci_high"])
                ax.plot([x, x], [lo, hi], color="black", lw=1.6, zorder=5)
                ax.plot([x - 0.08, x + 0.08], [lo, lo], color="black", lw=1.6, zorder=5)
                ax.plot([x - 0.08, x + 0.08], [hi, hi], color="black", lw=1.6, zorder=5)
                ax.scatter(
                    [x], [m], marker="o", s=90, color=c, edgecolor="black",
                    linewidths=1.1, zorder=6,
                )
                # Terse mean label beside each dot-and-whisker.
                ax.annotate(
                    f"{m:+.2f}", xy=(x + 0.16, m), va="center", ha="left",
                    fontsize=8, color="#222222", zorder=7,
                )
                # ADD-vs-arm effect (Hedges g + BH-q) atop the non-ADD arms.
                if arm != "ADD":
                    key = f"ADD_vs_{arm}"
                    if key in cset.index:
                        g = float(cset.loc[key, "hedges_g"])
                        q = float(cset.loc[key, "q_bh"])
                        ax.annotate(
                            f"g={g:+.1f}\nq={q:.0e}", xy=(x, yhi - 0.5 * pad),
                            va="top", ha="center", fontsize=6.5, color="#444444",
                            zorder=7,
                        )

            ax.set_xticks(range(len(ARM_ORDER)))
            ax.set_xticklabels(
                [f"{arm}\nKerchner {KERCHNER_FRACTION[arm]}" for arm in ARM_ORDER],
                fontsize=9,
            )
            ax.set_xlim(-0.6, len(ARM_ORDER) - 0.4)
            ax.set_ylim(ylo, yhi)
            ax.set_title(SET_TITLE[set_name], fontsize=12, fontweight="bold")
            if panel == 0:
                ax.set_ylabel("summary score (mean per-protein z)")

        fig.suptitle(
            "Four-arm signature scores — ADD-specific  ·  EXPLORATORY / "
            "confound-bounded (arms track Kerchner cohort)",
            fontsize=11, fontweight="bold", y=1.005,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.97))

    legend_fig = _build_legend(colors)
    return fig, legend_fig, colors


def _build_legend(colors: dict[str, str]) -> Figure:
    """Standalone key: arm swatches + the glyph meanings, as its own figure image."""
    with publication_style():
        lf, la = plt.subplots(figsize=(5.4, 2.6))
        la.axis("off")
        handles = [
            mlines.Line2D([], [], marker="o", ls="", ms=11, mfc=colors[arm],
                          mec="black", mew=1.0, label=f"{arm}  (Kerchner "
                          f"{KERCHNER_FRACTION[arm]})")
            for arm in ARM_ORDER
        ]
        handles += [
            mlines.Line2D([], [], marker="o", ls="", ms=11, mfc="white", mec="black",
                          mew=1.1, label="filled dot = arm mean z (± black 95% CI)"),
            mlines.Line2D([], [], marker="o", ls="", ms=7, mfc="#888888", mec="none",
                          alpha=0.5, label="small point = one sample (n=10/arm)"),
            mlines.Line2D([], [], color="black", lw=0.9, ls="--", alpha=0.6,
                          label="dashed line = z = 0"),
        ]
        la.legend(handles=handles, loc="center", frameon=False, fontsize=10,
                  handletextpad=0.7, labelspacing=0.9,
                  title="four-arm signature-score key", title_fontsize=11)
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
        "loaded four_arm cache %s: %d samples, %d per-arm rows, %d comparisons",
        args.cache_dir, len(sample_scores), len(per_arm), len(comparisons),
    )

    fig, legend_fig, colors = build_four_arm_figure(
        sample_scores, per_arm, comparisons, args.registry
    )
    artifacts = save_figure(
        fig, args.output_dir, args.base_name, legend_fig=legend_fig, dpi=300
    )
    LOGGER.info("arm colors: %s", colors)
    for p in (artifacts.svg, artifacts.png, artifacts.legend_svg, artifacts.legend_png):
        LOGGER.info("wrote %s", p)


if __name__ == "__main__":
    main()
