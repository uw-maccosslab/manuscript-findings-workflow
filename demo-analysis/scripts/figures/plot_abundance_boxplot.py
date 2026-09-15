"""Render the Stage-3 QC per-sample abundance box plot for the delivered state.

A publication-ready figure (dual-exported SVG + 300 DPI PNG, plus a separate
acquisition-order colorbar legend image) showing, for every experimental sample,
the distribution of its 2340 protein log2 abundances as one box. Boxes are drawn
in the Dataset's row order — already sorted by ``Run Order`` (acquisition order)
by the loader — and colored by run-order position on a sequential colormap, so
any residual left-to-right abundance drift over the run is visible and the
per-sample median line (black) shows whether the cohort is level.

SINGLE PANEL by necessity. This dataset arrives already log2, median-normalized,
and batch-corrected on Cohort; the pre-processing states (raw / normalized) are
UNAVAILABLE (state/DATA_DESCRIPTION.md, finding 0004), so there is nothing
earlier to stack. Controls (CP/KP) are absent (finding 0003): the 40 samples are
all experimental, so there is no separate control figure. The expectation for
already-normalized data is that the per-sample medians sit flat/aligned across
the run; this figure confirms that visually.

Inputs (all read-only, regenerable):
  results/qc_states/delivered_log/   the prepared delivered-state Dataset
                                     (load via common.dataset_io.load_dataset;
                                     materialized once by scripts/scratch/qc_prep.py)
  state/color_registry.json          palette source (no categorical dimension is
                                     drawn here, so the registry is not mutated)

Run:  ./.venv/bin/python scripts/figures/plot_abundance_boxplot.py
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np

# Import the project's shared figure helpers + the seeded abundance-boxplot module.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.dataset_io import load_dataset

from figures.abundance_boxplot import (
    AbundanceBoxplotScaleWarning,
    compute_sample_medians,
    save_abundance_boxplots,
)

LOGGER = logging.getLogger("plot_abundance_boxplot")

# The single available processing state (finding 0004). With one panel the per-panel
# left title is the figure's only title, so it carries the full processing-state
# statement (no separate suptitle, which would overlap the panel title on a short
# single-panel figure).
STATE_LABEL = (
    "delivered: log2, median-normalized, batch-corrected on cohort "
    "— single available state (finding 0004)"
)

DEFAULT_STATE_DIR = Path("results") / "qc_states" / "delivered_log"
DEFAULT_OUTPUT_DIR = Path("figures") / "qc" / "abundance-boxplot"
DEFAULT_BASE_NAME = "abundance-boxplot-delivered"
DEFAULT_REGISTRY = Path("state") / "color_registry.json"


def data_version(state_dir: Path) -> str:
    """Short content hash of the prepared state (abundances + sidecar); provenance."""
    hasher = hashlib.sha256()
    for name in ("abundances.npy", "dataset.json", "metadata.parquet"):
        hasher.update((state_dir / name).read_bytes())
    return hasher.hexdigest()[:16]


def log_median_flatness(medians: np.ndarray) -> None:
    """Report per-sample median spread + run-order drift (the QC read of the figure)."""
    finite = medians[np.isfinite(medians)]
    positions = np.arange(finite.size)
    slope = float(np.polyfit(positions, finite, 1)[0])
    drift = slope * (finite.size - 1)
    corr = float(np.corrcoef(positions, finite)[0, 1])
    LOGGER.info(
        "per-sample median (log2): min=%.4f max=%.4f spread=%.4f std=%.4f",
        finite.min(),
        finite.max(),
        finite.max() - finite.min(),
        finite.std(),
    )
    LOGGER.info(
        "run-order drift: slope=%.5f/sample, total=%.4f log2 over the run, "
        "corr(position, median)=%.4f",
        slope,
        drift,
        corr,
    )


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(name)s: %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state-dir", type=Path, default=DEFAULT_STATE_DIR)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--base-name", type=str, default=DEFAULT_BASE_NAME)
    args = parser.parse_args()

    LOGGER.info("data_version=%s", data_version(args.state_dir))

    dataset = load_dataset(args.state_dir)
    n_samples, n_features = dataset.abundances.shape
    LOGGER.info(
        "loaded %s: %d samples x %d features, scale=%s",
        args.state_dir,
        n_samples,
        n_features,
        dataset.scale,
    )
    log_median_flatness(compute_sample_medians(dataset))

    # Capture the template's scale warning so we can report explicitly that a log2
    # state renders silently (the guard fires only on a non-log scale).
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", AbundanceBoxplotScaleWarning)
        artifacts = save_abundance_boxplots(
            {STATE_LABEL: dataset},
            args.output_dir,
            args.base_name,
            feature_type="protein",
            title=None,
            registry_path=args.registry,
            persist_colors=False,
        )
    scale_warnings = [
        w for w in caught if issubclass(w.category, AbundanceBoxplotScaleWarning)
    ]
    if scale_warnings:
        for w in scale_warnings:
            LOGGER.warning("scale warning: %s", w.message)
    else:
        LOGGER.info("no scale warning (state is log2, in LOG_SCALES) — silent")

    LOGGER.info("wrote %s", artifacts.svg)
    LOGGER.info("wrote %s", artifacts.png)
    LOGGER.info("wrote %s", artifacts.legend_svg)
    LOGGER.info("wrote %s", artifacts.legend_png)


if __name__ == "__main__":
    main()
