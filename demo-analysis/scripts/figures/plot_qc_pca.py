"""Render the Stage-3 QC PCA family for the delivered proteomics state.

Two figures from ONE script / ONE run, both on the ``delivered_log`` QC-state Dataset
(40 experimental samples x 2340 proteins; scale ``log2``; median-normalized AND
batch-corrected on Cohort via limma ``removeBatchEffect`` protecting Condition):

  (A) ``pca-delivered-by-cohort``  -- points colored by Cohort (collection site):
      the BATCH AXIS. Because the delivered matrix was batch-corrected on Cohort, the
      cohorts should be largely mixed/overlapping; a residual cohort split would signal
      incomplete correction.
  (B) ``pca-delivered-by-group``   -- points colored by Condition (diagnostic group
      ADD/PDD/PDCN/HCN): the BIOLOGY. Do the diagnostic groups separate? (Caveat, for
      the caption: group is confounded with cohort per finding 0001, and ADD/PDCN are
      single-cohort so their cohort correction is extrapolated per finding 0004 -- so
      any group separation cannot be attributed to biology vs residual site effect.)

SINGLE-STATE PCA, not a processing-state series: the pre-normalization /
pre-batch-correction matrices are unavailable (finding 0004), so the usual
raw -> normalized -> batch-corrected series collapses to this one delivered state.
Controls (CP/KP) are absent from the matrix (finding 0003): all 40 samples are
experimental, so no reference greying is needed.

Both figures render PC1/PC2 and PC3/PC4 panels with % variance explained on the axes and
per-PC Kruskal-Wallis group-separation p-values on the marginals. Colors route through
``state/color_registry.json`` (cohort / group dimensions, both <= 4 categories, within
the >8-category guard). Each figure dual-exports SVG + 300 DPI PNG plus a separate
swatch-legend image. Seeds from the ``pca-plot`` lib template via ``figures.pca``.

Run:  ./.venv/bin/python scripts/figures/plot_qc_pca.py
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

# The project root is two levels up from this file (…/scripts/figures/plot_qc_pca.py);
# put …/scripts on sys.path so ``common.*`` and ``figures.*`` resolve.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.dataset_io import load_dataset

from figures.pca import save_pca

LOGGER = logging.getLogger("plot_qc_pca")

# The single delivered QC processing state (there is no state series for this dataset).
DEFAULT_STATE_DIR = Path("results") / "qc_states" / "delivered_log"
DEFAULT_OUTPUT_DIR = Path("figures") / "qc" / "pca"
DEFAULT_REGISTRY = Path("state") / "color_registry.json"

# The processing state, stated verbatim in every title/caption per the QC convention.
STATE_LABEL = "delivered: log2, median-normalized, batch-corrected on cohort"

FEATURE_TYPE = "protein"

# One entry per figure in the family: (base_name, color_by column, registry category,
# legend title, a short "what this coloring is" phrase for the suptitle).
FIGURES = (
    (
        "pca-delivered-by-cohort",
        "Cohort",
        "cohort",
        "Cohort (collection site)",
        "colored by cohort (batch axis)",
    ),
    (
        "pca-delivered-by-group",
        "Condition",
        "group",
        "Diagnostic group (Condition)",
        "colored by diagnostic group (biology)",
    ),
)


def data_version(state_dir: Path) -> str:
    """Return a short content hash of the loaded QC-state, for figure provenance.

    Hashes the state's abundance matrix, metadata table, and sidecar so a change to any
    of them re-versions the figure (the staleness machinery keys on this).
    """
    h = hashlib.sha256()
    for name in ("abundances.npy", "metadata.parquet", "dataset.json"):
        path = state_dir / name
        if not path.exists():
            raise FileNotFoundError(f"QC-state file missing: {path}")
        h.update(path.read_bytes())
    return h.hexdigest()[:16]


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(name)s: %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state-dir", type=Path, default=DEFAULT_STATE_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()

    dv = data_version(args.state_dir)
    LOGGER.info("data_version=%s (state=%s)", dv, args.state_dir)

    dataset = load_dataset(args.state_dir)
    n = int(dataset.abundances.shape[0])
    LOGGER.info(
        "loaded Dataset: %d samples x %d %ss, scale=%s",
        n,
        int(dataset.abundances.shape[1]),
        FEATURE_TYPE,
        dataset.scale,
    )

    for base_name, color_by, category, legend_title, coloring in FIGURES:
        title = f"PCA of experimental samples (n={n}), {coloring}\n{STATE_LABEL}"
        artifacts = save_pca(
            dataset,
            color_by,
            args.output_dir,
            base_name,
            category=category,
            feature_type=FEATURE_TYPE,
            title=title,
            legend_title=legend_title,
            registry_path=args.registry,
            dpi=args.dpi,
        )
        LOGGER.info("wrote %s", artifacts.svg)
        LOGGER.info("wrote %s", artifacts.png)
        LOGGER.info("wrote %s", artifacts.legend_svg)
        LOGGER.info("wrote %s", artifacts.legend_png)


if __name__ == "__main__":
    main()
