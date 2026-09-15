"""Stage-3 QC "prep-once": materialize the canonical QC state matrices ONCE.

Loads the Stage-4 analysis subset (40 samples x 2340 proteins, via
``scripts/scratch/project_loader.py``) and writes the processing-state ``Dataset``s
that every downstream QC figure job consumes, so PCA / sample-correlation /
abundance-boxplot / CV / dynamic-range each just ``common.dataset_io.load_dataset(...)``
the state they need instead of re-loading and re-deriving it (see the "prep-once"
pattern described in this project's ``CLAUDE.md``).

**THIS DATASET HAS ONLY ONE DELIVERED PROCESSING STATE.** Per
``state/DATA_DESCRIPTION.md`` (and finding 0004), ``data/protein-data.csv`` arrives
already log2, already median-normalized, and already batch-corrected on ``Cohort``
(limma ``removeBatchEffect``, protecting ``Condition``). The pre-normalization /
pre-batch-correction matrices were **not delivered** — they do not exist in this
project. So the usual raw -> normalized -> batch-corrected QC state series
**collapses to a single delivered state**: there is nothing upstream of it to
materialize. This script therefore:

  * does **NOT** call ``common.normalize.normalize()`` / ``log2_transform()`` — the
    data is already log2 (linear-scale input required; would raise / this is not a
    "re-normalize to be safe" situation, it is a correctness violation to try);
  * does **NOT** call ``common.batch_correct.combat_correct()`` — the data is already
    batch-corrected; running it again would double-correct;
  * DOES call ``common.batch_correct.assess_batch_confounding()`` on the delivered
    matrix (``Cohort`` vs ``Condition``) to RE-CONFIRM the confound already recorded
    as findings 0001 / 0004. This is a corroborating check, not a new finding.

Writes, under ``--output-dir`` (default ``results/qc_states/``):

  * ``delivered_log/``      — the analysis-subset ``Dataset`` as-is, ``scale="log2"``.
                               Consumed by: PCA, sample-correlation, abundance-boxplot.
  * ``delivered_linear/``   — ``common.normalize.to_linear(delivered_log)``,
                               ``scale="linear"``. This is the de-logged IMAGE of the
                               already-normalized, already-batch-corrected matrix, NOT
                               raw/unnormalized intensity. Consumed by: CV,
                               dynamic-range (both require a linear scale).
  * ``batch_confounding_report.json`` — the ``ConfoundingReport``(s) from
                               ``assess_batch_confounding``.
  * ``README.md``           — the single-state explanation above, for anyone landing
                               in ``results/qc_states/`` without this docstring.

No stochastic step runs in this script (loading, de-logging, and a contingency-table
association statistic are all deterministic), so no random seed is required or
recorded.

Run:
    ./.venv/bin/python scripts/scratch/qc_prep.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import sys
import warnings
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# `common.*` lives at scripts/common/ -- make `scripts/` importable regardless of how
# this file is invoked (direct script, `-m`, or imported by pytest), matching the
# pattern used by scripts/scratch/project_loader.py and
# scripts/figures/plot_group_x_cohort_confound.py.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.batch_correct import (
    ConfoundingReport,
    assess_batch_confounding,
)
from common.data_loading import Dataset
from common.dataset_io import save_dataset
from common.normalize import to_linear

# project_loader.py sits next to this file in scripts/scratch/; running this script
# directly (or collecting it via pytest, which has no __init__.py in scratch/) puts
# scripts/scratch/ on sys.path[0] automatically, so a bare import resolves -- the same
# pattern scripts/scratch/test_project_loader.py already relies on.
from project_loader import (
    ANALYSIS_DROP_FEATURE_IDS,
    load_analysis_subset,
    load_paired_dataset,
)

__script_meta__: dict[str, object] = {
    "task": "qc-prep",
    "kind": "analysis",
    "provides": [],
    "uses": [
        "project_loader",
        "common.data_loading",
        "common.dataset_io",
        "common.normalize",
        "common.batch_correct",
    ],
    "seeded_from": None,
    "description": (
        "Stage-3 QC 'prep-once': materializes the canonical QC state Dataset(s) under "
        "results/qc_states/ so figure jobs load->subset->plot instead of re-deriving. "
        "This project has only ONE delivered processing state (log2, pre-normalized, "
        "pre-batch-corrected on Cohort -- state/DATA_DESCRIPTION.md, finding 0004), "
        "so no normalize()/combat_correct() runs here; it writes delivered_log "
        "(as-is) and delivered_linear (to_linear de-logged image, for the "
        "linear-scale CV / dynamic-range figures), plus a "
        "batch_confounding_report.json that RE-CONFIRMS the Cohort<->Condition "
        "confound already recorded as findings 0001/0004 (no new finding is "
        "created)."
    ),
}

logger = logging.getLogger(__name__)

DEFAULT_DATA = Path("data/protein-data.csv")
DEFAULT_METADATA = Path("data/metadata.csv")
DEFAULT_PREPARED_DIR = Path("results/prepared")
DEFAULT_OUTPUT_DIR = Path("results/qc_states")
DEFAULT_BATCH_COLUMN = "Cohort"
DEFAULT_COVARIATE_COLUMNS: tuple[str, ...] = ("Condition",)
# Matches common.batch_correct.assess_batch_confounding's own default warn threshold
# (kept as an explicit, overridable constant here rather than imported, since it is not
# part of that module's public `provides`).
DEFAULT_WARN_THRESHOLD = 0.5
EXPECTED_N_SAMPLES = 40
# data/metadata.csv classifies 48 rows into 40 experimental + 8 control
# (CP01-04, KP01-04 pooled references); see state/METADATA.md.
EXPECTED_N_CONTROL = 8

DELIVERED_LOG_DIRNAME = "delivered_log"
DELIVERED_LINEAR_DIRNAME = "delivered_linear"
CONFOUNDING_REPORT_FILENAME = "batch_confounding_report.json"
README_FILENAME = "README.md"

# Findings this script's confounding assessment RE-CONFIRMS -- it must never create a
# new finding; this constant is only for the corroboration note in the written report
# and the log line, not for writing to findings/.
CORROBORATES_FINDINGS: tuple[str, ...] = ("0001", "0004")

_README_TEXT = """\
# QC processing states

This project's proteomics data has **only one delivered processing state**. There is
no raw -> normalized -> batch-corrected series to materialize here, because the
pre-normalization and pre-batch-correction matrices were never delivered to this
project (see `state/DATA_DESCRIPTION.md` and finding
[0004](../../findings/0004-batch-correction-extrapolated-single-cohort.md)).

## What "delivered" means

`data/protein-data.csv` arrives:

- **log2**-scale (not linear);
- **already median-normalized** upstream (per-sample medians span 0.15 log2 units);
- **already batch-corrected on `Cohort`** upstream, via
  `limma::removeBatchEffect(..., design = ~Condition, batch = Cohort)` (protecting
  `Condition`) -- see finding 0004. The pre-correction matrix (`nobatchProt` in the
  scientist's R pipeline) is **not available** in this project.

Consequently `scripts/scratch/qc_prep.py` does **not** call
`common.normalize.normalize()` / `log2_transform()` (the data is already log2;
those functions require a **linear**-scale input and would refuse to run -- this is a
correctness guard, not an oversight) and does **not** run ComBat /
`removeBatchEffect` again (the matrix is already corrected; re-running it would
double-correct).

## The two states materialized here

- **`delivered_log/`** (`scale="log2"`) -- the Stage-4 analysis-subset `Dataset` (40
  samples x 2340 proteins) exactly as loaded, no transform applied. Consumed by: PCA,
  sample-correlation, abundance-boxplot.
- **`delivered_linear/`** (`scale="linear"`) --
  `common.normalize.to_linear(delivered_log)`, the mathematical de-log (`2**x`) of
  the delivered matrix. Consumed by: CV, dynamic-range.

**`delivered_linear` is NOT raw / unnormalized intensity.** It is the de-logged image
of the already-normalized, already-batch-corrected delivered matrix -- a linear-scale
view materialized only because coefficient-of-variation and dynamic-range are
undefined on a log scale (a per-feature mean near 0 makes `std/mean` explode, and
"detected" is a linear-scale concept). Do not read or label it as raw intensity.

## Controls are absent

The 8 pooled-reference control samples (`CP01`-`CP04`, `KP01`-`KP04`) exist only in
`data/metadata.csv`; they have no columns in `data/protein-data.csv` (finding
[0003](../../findings/0003-contrasts-class-imbalanced-small-n.md)). Both states above
therefore contain **only the 40 experimental samples** -- every QC figure built from
them is experimental-only; there is no pooled-reference / run-drift comparison
available from this matrix.

## Batch confounding re-confirmation

`batch_confounding_report.json` (written by the same script) re-runs
`common.batch_correct.assess_batch_confounding(delivered_log, batch_column="Cohort",
covariate_columns=["Condition"])` on the delivered matrix. This **RE-CONFIRMS** the
Cohort <-> Condition confound already recorded in findings
[0001](../../findings/0001-group-aliased-with-cohort.md) and
[0004](../../findings/0004-batch-correction-extrapolated-single-cohort.md) -- it is a
corroborating check computed for provenance, not a new finding.

## Regenerating

```
./.venv/bin/python scripts/scratch/qc_prep.py
```

See `--help` for parameters. Figure-generator jobs should
`common.dataset_io.load_dataset(...)` the state they need from this directory rather
than re-deriving it.
"""


@dataclass(frozen=True)
class DeliveredStates:
    """The two canonical QC states materialized from the delivered analysis subset.

    Attributes
    ----------
    log:
        The delivered analysis-subset Dataset, unchanged, ``scale="log2"``.
    linear:
        ``to_linear(log)`` -- the de-logged image, ``scale="linear"``. NOT raw
        intensity (see the module docstring).
    """

    log: Dataset
    linear: Dataset


def sha256_16(path: Path) -> str:
    """First 16 hex characters of ``path``'s sha256 -- a short data-version stamp."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def build_delivered_states(analysis_subset: Dataset) -> DeliveredStates:
    """Return the two canonical QC states from the delivered analysis-subset Dataset.

    ``analysis_subset`` MUST already be on scale ``"log2"`` -- this project's data
    arrives pre-logged, pre-normalized, pre-batch-corrected (state/DATA_DESCRIPTION.md,
    finding 0004), so there is nothing to normalize or correct here. Fails loud if the
    scale is wrong, if the input contains a non-finite value, or if the de-logged
    result is not strictly positive / shape-preserving (each would mean the delivered
    matrix or ``to_linear`` broke an invariant this whole prep-once step depends on).

    Returns
    -------
    DeliveredStates
        ``log`` is ``analysis_subset`` itself (not copied further -- ``to_linear``
        already returns an independent Dataset for ``linear``); ``linear`` is its
        de-logged image.
    """
    if analysis_subset.scale != "log2":
        raise ValueError(
            f"build_delivered_states expects the delivered analysis subset on scale "
            f"'log2' (this project's data arrives pre-logged -- see "
            f"state/DATA_DESCRIPTION.md); got scale={analysis_subset.scale!r}. This "
            f"project must never re-normalize or re-log the delivered matrix."
        )
    if not np.isfinite(analysis_subset.abundances).all():
        n_bad = int((~np.isfinite(analysis_subset.abundances)).sum())
        raise ValueError(
            f"delivered_log abundances contain {n_bad} non-finite value(s); the "
            f"Stage-3 loader (project_loader.load_paired_dataset) should already "
            f"guarantee a complete matrix (missing='error')."
        )

    delivered_linear = to_linear(analysis_subset)

    if delivered_linear.scale != "linear":
        raise ValueError(  # pragma: no cover - guaranteed by to_linear's own contract
            f"to_linear returned scale={delivered_linear.scale!r}, expected 'linear'."
        )
    if delivered_linear.abundances.shape != analysis_subset.abundances.shape:
        raise ValueError(
            f"to_linear changed the matrix shape "
            f"{analysis_subset.abundances.shape} -> "
            f"{delivered_linear.abundances.shape}; expected it preserved."
        )
    if not bool(np.all(delivered_linear.abundances > 0)):
        raise ValueError(
            "delivered_linear contains non-positive value(s); the de-logged image of "
            "a finite log2 matrix must be strictly positive."
        )
    return DeliveredStates(log=analysis_subset, linear=delivered_linear)


def assess_confounding(
    dataset: Dataset,
    *,
    batch_column: str,
    covariate_columns: Sequence[str],
    warn_threshold: float,
) -> tuple[list[ConfoundingReport], list[dict[str, str]]]:
    """Run ``assess_batch_confounding``, capturing any warning(s) it raises.

    Returns ``(reports, captured_warnings)``: ``reports`` is one
    :class:`~common.batch_correct.ConfoundingReport` per covariate (input order);
    ``captured_warnings`` is a JSON-serializable summary
    (``[{"category": ..., "message": ...}, ...]``) of every warning emitted (a
    :class:`~common.batch_correct.BatchConfoundingWarning` on perfect/strong
    confounding). Any exception the underlying call raises (e.g. a missing or
    NaN-containing column) is NOT caught here -- it propagates (fail loud).
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        reports = assess_batch_confounding(
            dataset,
            batch_column,
            covariate_columns,
            warn_threshold=warn_threshold,
        )
    captured = [
        {"category": w.category.__name__, "message": str(w.message)} for w in caught
    ]
    return reports, captured


def report_to_jsonable(report: ConfoundingReport) -> dict[str, object]:
    """Convert one :class:`ConfoundingReport` into a JSON-serializable ``dict``."""
    crosstab_int = report.crosstab.astype(int)
    crosstab: dict[str, dict[str, int]] = {
        str(batch_level): {
            str(cov_level): int(count) for cov_level, count in row.items()
        }
        for batch_level, row in crosstab_int.to_dict(orient="index").items()
    }
    return {
        "batch_column": report.batch_column,
        "covariate_column": report.covariate_column,
        "cramers_v": report.cramers_v,
        "perfectly_confounded": report.perfectly_confounded,
        "message": report.message,
        "crosstab": crosstab,
    }


def write_confounding_report(
    reports: Sequence[ConfoundingReport],
    captured_warnings: Sequence[dict[str, str]],
    *,
    output_path: Path,
    data_version: dict[str, str],
) -> Path:
    """Write the confounding assessment as JSON to ``output_path``; return it.

    Includes an explicit ``corroborates_findings`` / ``note`` field so the file is
    self-explanatory about NOT being a new finding -- it re-confirms findings 0001 and
    0004 on the delivered (already-corrected) matrix.
    """
    payload: dict[str, object] = {
        "data_version": dict(data_version),
        "reports": [report_to_jsonable(r) for r in reports],
        "warnings_captured": list(captured_warnings),
        "corroborates_findings": list(CORROBORATES_FINDINGS),
        "note": (
            "This is a RE-CONFIRMATION of the Cohort<->Condition confound already "
            "recorded as findings 0001 (group-aliased-with-cohort) and 0004 "
            "(batch-correction-extrapolated-single-cohort). It produces no new "
            "correction (the delivered matrix is already batch-corrected -- see "
            "state/DATA_DESCRIPTION.md) and is not itself a new finding."
        ),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return output_path


def write_readme(output_dir: Path) -> Path:
    """Write the single-processing-state README into ``output_dir``; return its path."""
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / README_FILENAME
    path.write_text(_README_TEXT, encoding="utf-8")
    return path


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--prepared-dir", type=Path, default=DEFAULT_PREPARED_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--batch-column", type=str, default=DEFAULT_BATCH_COLUMN)
    parser.add_argument(
        "--covariate-columns",
        type=str,
        nargs="+",
        default=list(DEFAULT_COVARIATE_COLUMNS),
        help="Covariate(s) of interest to check against --batch-column.",
    )
    parser.add_argument("--warn-threshold", type=float, default=DEFAULT_WARN_THRESHOLD)
    parser.add_argument("--expected-n-samples", type=int, default=EXPECTED_N_SAMPLES)
    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
    )
    return parser.parse_args(argv)


def run(args: argparse.Namespace) -> None:
    """Execute the QC prep-once pipeline for the given parsed arguments."""
    logger.info(
        "qc_prep: data=%s metadata=%s prepared_dir=%s output_dir=%s batch_column=%s "
        "covariate_columns=%s warn_threshold=%s expected_n_samples=%s",
        args.data,
        args.metadata,
        args.prepared_dir,
        args.output_dir,
        args.batch_column,
        args.covariate_columns,
        args.warn_threshold,
        args.expected_n_samples,
    )
    logger.info(
        "no stochastic step runs in this script (deterministic load/de-log/crosstab); "
        "no random seed is required or recorded"
    )

    data_version = {
        "data": sha256_16(args.data),
        "metadata": sha256_16(args.metadata),
    }
    logger.info("data_version=%s", data_version)

    full_dataset = load_paired_dataset(
        args.data,
        args.metadata,
        prepared_dir=args.prepared_dir,
        expected_n_control=EXPECTED_N_CONTROL,
    )
    logger.info(
        "loaded full paired dataset: %d sample(s) x %d feature(s), scale=%s",
        *full_dataset.abundances.shape,
        full_dataset.scale,
    )

    analysis_subset = load_analysis_subset(full_dataset)
    expected_n_features = full_dataset.abundances.shape[1] - len(
        ANALYSIS_DROP_FEATURE_IDS
    )
    expected_shape = (args.expected_n_samples, expected_n_features)
    if analysis_subset.abundances.shape != expected_shape:
        raise ValueError(
            f"analysis subset shape {analysis_subset.abundances.shape} does not "
            f"match the expected {expected_shape} "
            f"({args.expected_n_samples} experimental samples x "
            f"{full_dataset.abundances.shape[1]} raw features minus "
            f"{len(ANALYSIS_DROP_FEATURE_IDS)} analysis-excluded features)."
        )
    logger.info(
        "loaded analysis subset: %d sample(s) x %d feature(s), scale=%s",
        *analysis_subset.abundances.shape,
        analysis_subset.scale,
    )

    states = build_delivered_states(analysis_subset)
    logger.info(
        "delivered_log: %d sample(s) x %d feature(s), scale=%s",
        *states.log.abundances.shape,
        states.log.scale,
    )
    logger.info(
        "delivered_linear: %d sample(s) x %d feature(s), scale=%s",
        *states.linear.abundances.shape,
        states.linear.scale,
    )

    log_path = save_dataset(states.log, args.output_dir / DELIVERED_LOG_DIRNAME)
    linear_path = save_dataset(
        states.linear, args.output_dir / DELIVERED_LINEAR_DIRNAME
    )
    logger.info("wrote %s", log_path)
    logger.info("wrote %s", linear_path)

    reports, captured_warnings = assess_confounding(
        states.log,
        batch_column=args.batch_column,
        covariate_columns=args.covariate_columns,
        warn_threshold=args.warn_threshold,
    )
    for report in reports:
        logger.info(
            "confounding %s x %s: Cramer's V=%.3f perfectly_confounded=%s",
            report.batch_column,
            report.covariate_column,
            report.cramers_v,
            report.perfectly_confounded,
        )
        logger.info("%s", report.message)
    for captured in captured_warnings:
        logger.warning("%s: %s", captured["category"], captured["message"])
    logger.info(
        "this RE-CONFIRMS the %s<->%s confound already recorded as findings %s -- "
        "no new finding is created by this script",
        args.batch_column,
        list(args.covariate_columns),
        ", ".join(CORROBORATES_FINDINGS),
    )

    report_path = write_confounding_report(
        reports,
        captured_warnings,
        output_path=args.output_dir / CONFOUNDING_REPORT_FILENAME,
        data_version=data_version,
    )
    logger.info("wrote %s", report_path)

    readme_path = write_readme(args.output_dir)
    logger.info("wrote %s", readme_path)


def main(argv: Sequence[str] | None = None) -> None:
    args = _parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, str(args.log_level)),
        format="%(levelname)s %(name)s: %(message)s",
    )
    run(args)


if __name__ == "__main__":
    main()
