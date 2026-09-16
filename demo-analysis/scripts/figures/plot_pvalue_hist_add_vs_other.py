"""Render the raw-p-value calibration histogram for the PRIMARY ADD-vs-Other DE run.

A publication-ready figure (dual-exported SVG + 300 DPI PNG, plus a separate swatch
legend image) of the distribution of the 2340 **raw** per-feature p-values from the
primary differential-abundance contrast — ADD (n=10) vs pooled Other {PDD, PDCN, HCN}
(n=30), ``method="moderated"`` (limma-style empirical Bayes) — with the uniform null
drawn at density 1.0 for reference.

This is the honesty / calibration check on the DE run. The observed shape is NOT a clean
"uniform body + spike at 0": there is a strong monotone spike at p≈0 (real small-p
enrichment) BUT the null region is DEPLETED (the p>0.5 tail sits below the uniform line,
not on it). A spike plus a suppressed null tail is the signature of ANTICONSERVATIVE
p-values — here from variance deflation (the DE was run on an already-batch-corrected
matrix, finding 0004) compounded by the ADD/Kerchner single-cohort confound (finding
0001). Consequence: the BH-q hit count is inflated; trust effect sizes + CIs over the
q-values. The on-canvas annotation carries only the two load-bearing numbers (the bin[0]
enrichment and the p>0.5 tail level); the full reading lives in the finding prose /
caption, not on the plot.

Input (read-only, produced by scripts/scratch/de_add_vs_other.py):
  results/differential_abundance/add-vs-other-moderated/contrast_table.csv
      the moderated contrast table; column ``p`` = raw p-values (2340 features)
  state/color_registry.json   Okabe-Ito palette (colors read; NOT persisted here —
                              a p-value-histogram label is figure-local)

Run:  ./.venv/bin/python scripts/figures/plot_pvalue_hist_add_vs_other.py
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd

# `common.*` and the sibling `pvalue_hist` module live under scripts/ -- make scripts/
# importable regardless of how this file is invoked (matches plot_cv_distribution.py).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.figure_io import save_figure

from figures.pvalue_hist import plot_pvalue_histogram

LOGGER = logging.getLogger("plot_pvalue_hist_add_vs_other")

DEFAULT_TABLE = (
    Path("results")
    / "differential_abundance"
    / "add-vs-other-moderated"
    / "contrast_table.csv"
)
DEFAULT_OUTPUT_DIR = (
    Path("figures") / "analysis" / "differential-abundance" / "add-vs-other"
)
DEFAULT_BASE_NAME = "pvalue-hist-add-vs-other"
DEFAULT_REGISTRY = Path("state") / "color_registry.json"

# Figure-local distribution label (NOT persisted; see module docstring). Names the
# contrast + method so the legend entry ("... (n=2340, π0≈0.53)") is self-describing.
DIST_LABEL = "ADD vs Other (moderated)"

TITLE = "Raw p-value calibration — ADD vs Other (moderated DE)"
# Mandatory short caveat marker: the load-bearing shape read, kept to one line so it
# reads off the canvas without becoming a paragraph (the interpretation lives in the
# finding prose / caption, per the annotation budget).
CAVEAT = "spike at p≈0 with a DEPLETED null tail → anticonservative p-values"


def data_version(table: Path) -> str:
    """Short content hash of the input contrast table, for provenance."""
    return hashlib.sha256(table.read_bytes()).hexdigest()[:16]


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(name)s: %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--base-name", type=str, default=DEFAULT_BASE_NAME)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    args = parser.parse_args()

    dv = data_version(args.table)
    LOGGER.info("data_version(contrast_table.csv)=%s", dv)

    frame = pd.read_csv(args.table)
    if "p" not in frame.columns:
        raise ValueError(
            f"{args.table} has no 'p' column; got {list(frame.columns)}. "
            f"This figure needs the RAW p-values."
        )
    pvals = frame["p"].to_numpy(dtype=float)
    LOGGER.info(
        "loaded %d rows; finite p=%d, NaN p=%d, range=[%.3g, %.3g]",
        len(frame),
        int(np.isfinite(pvals).sum()),
        int(np.isnan(pvals).sum()),
        float(np.nanmin(pvals)),
        float(np.nanmax(pvals)),
    )

    # Single distribution -> the template draws the shaded histogram, the gray uniform
    # reference at 1.0, and (show_pi0=True) the π0 null-floor line. Colors are not
    # persisted: the label is figure-local, so nothing is written to the color registry.
    plot = plot_pvalue_histogram(
        {DIST_LABEL: pvals},
        n_bins=20,
        show_pi0=True,
        title=None,  # title + caveat added below so the caveat is a smaller note
        legend_title="Contrast (raw p)",
        registry_path=args.registry,
        persist_colors=False,
    )

    # Report the load-bearing calibration numbers (also placed on the canvas below).
    m = plot.result.counts[DIST_LABEL]
    finite = pvals[np.isfinite(pvals)]
    counts, _edges = np.histogram(finite, bins=np.linspace(0, 1, 21))
    bin0_density = counts[0] / m / 0.05  # density of p in [0, 0.05]
    tail_density = float(np.mean(finite > 0.5)) / 0.5  # avg density over p > 0.5
    pi0 = plot.result.pi0[DIST_LABEL]
    LOGGER.info(
        "calibration: n=%d, bin[0] density=%.3f (%.2fx uniform), "
        "p>0.5 tail density=%.3f (%.2fx uniform), π0≈%.3f",
        m,
        bin0_density,
        bin0_density,
        tail_density,
        tail_density,
        pi0,
    )
    LOGGER.info("color_map=%s", plot.color_map)

    fig = plot.figure
    ax = fig.axes[0]
    ax.set_ylim(0.0, max(7.2, bin0_density * 1.08))

    # (1) bin[0] enrichment: point at the tall first bar (real small-p signal).
    ax.annotate(
        f"p∈[0,0.05]: {bin0_density:.1f}× uniform",  # noqa: RUF001
        xy=(0.05, bin0_density),
        xytext=(0.20, bin0_density * 0.90),
        fontsize=10,
        color="#333333",
        arrowprops={"arrowstyle": "->", "color": "#333333", "linewidth": 1.0},
    )
    # (2) label the uniform-null reference line (drawn by the template at density 1.0).
    ax.text(
        0.28,
        1.06,
        "uniform null (density = 1)",
        fontsize=9,
        color="gray",
        va="bottom",
        ha="left",
    )
    # (3) the depleted null tail: annotate the region below the uniform line. The
    #     template already marks the π0 floor line (== the p>0.5 tail level at λ=0.5)
    #     with "π0≈" at the right; this note names what that suppression means.
    ax.annotate(
        f"p>0.5 tail ≈ {tail_density:.2f}× uniform\n(null region suppressed)",  # noqa: RUF001
        xy=(0.75, tail_density),
        xytext=(0.42, tail_density + 1.9),
        fontsize=10,
        color="#333333",
        ha="left",
        va="bottom",
        arrowprops={"arrowstyle": "->", "color": "#333333", "linewidth": 1.0},
    )

    # Short bold title + a smaller italic caveat marker (the load-bearing shape read).
    fig.subplots_adjust(top=0.83)
    fig.suptitle(TITLE, y=0.99, fontsize=14, weight="bold")
    fig.text(
        0.5,
        0.90,
        CAVEAT,
        ha="center",
        va="top",
        fontsize=9,
        color="#444444",
        style="italic",
    )

    artifacts = save_figure(
        fig, args.output_dir, args.base_name, legend_fig=plot.legend_figure
    )
    LOGGER.info("wrote %s", artifacts.svg)
    LOGGER.info("wrote %s", artifacts.png)
    LOGGER.info("wrote %s", artifacts.legend_svg)
    LOGGER.info("wrote %s", artifacts.legend_png)


if __name__ == "__main__":
    main()
