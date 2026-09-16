"""BORUTA all-relevant feature selection: ADD vs Other, plasma-EV study — COMPUTE.

Call site for the project copy of the vetted ``lib/analysis/boruta`` template
(``scripts/scratch/boruta.py``, ``boruta_select``). This is the *compute* half of the
compute-once / render-many split (conventions/results-cache.md): it runs Boruta ONCE on
the whole experimental matrix and caches the frozen ``BorutaResult`` via
``result_io.save_cached_result`` under ``results/boruta/<fingerprint>/`` so the figure
generator and the summary read the cache instead of re-running the (slow) forest.

What it does:

  1. loads the 40-sample x 2340-protein log2 analysis subset via ``project_loader``
     (controls already excluded; already median-normalized + ComBat/limma
     batch-corrected on Cohort; complete, no NaN — state/DATA_DESCRIPTION.md);
  2. builds the binary target **ADD vs Other** on ``Condition`` (ADD=1 -> "ADD"
     positive, {PDD,PDCN,HCN}=0 -> "Other") as a *categorical* metadata column, so the
     template infers a **classification** task from the target dtype (a numeric 0/1
     column would mis-infer regression);
  3. runs ``boruta_select`` on the WHOLE matrix — no CV, no external null: Boruta's
     shadow (permuted-copy) features ARE its built-in null (conventions/statistics.md
     §B; feeding it the ComBat-corrected matrix is legitimate — this is selection, not a
     per-feature significance test);
  4. CACHES the ``BorutaResult`` and prints the ``ResultMeta`` (fingerprint / path /
     params) for registration in ``results/manifest.md``.

TUNING (recorded — every non-default enters ``provenance.params`` + the fingerprint):
the default run is ``max_iter=100, alpha=0.05, perc=100`` (perc=100 = the MAX shadow
importance, the strictest noise floor). On this p>>n matrix (2340 features, 40 samples,
10 positives) the max-shadow floor is stringent, so if too few features Confirm, raise
``--max-iter`` (200) and/or lower ``--perc`` toward 95/90; if the Confirmed set is too
large for a tractable downstream enrichment, lower ``--alpha`` to 0.01. Whatever is run
is what is cached and reported.

CAVEATS (state/DATA_DESCRIPTION.md; findings 0001, 0004):

  * **EXPLORATORY.** Small n (10 positives). Boruta's shadow null addresses
    relevance-vs-noise, not sampling stability at this n.
  * **COHORT CONFOUND (hard).** ADD is 100% Kerchner cohort (finding 0001): among these
    40 samples "ADD" and "Kerchner" are the SAME partition. Boruta confirms features
    relevant to the ADD-vs-Other partition, which is ALIASED with collection site — so
    "all-relevant to ADD" is not "all-relevant to disease". The shadow null does NOT
    resolve this confound.

Non-interactive, parameterized (argparse), seeded. Fail-loud at every boundary.

Run::

    ./.venv/bin/python scripts/scratch/boruta_add_vs_other.py
    ./.venv/bin/python scripts/scratch/boruta_add_vs_other.py --max-iter 200 --perc 90
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

# Make `scripts/` and `scripts/scratch/` importable regardless of invocation, mirroring
# classify_add_vs_other.py / de_add_vs_other.py.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from boruta_selection import boruta_select
from common.data_loading import Dataset
from common.result_io import ResultMeta, save_cached_result
from project_loader import load_analysis_subset, load_paired_dataset

logger = logging.getLogger(__name__)

__script_meta__: dict[str, object] = {
    "task": "boruta-add-vs-other",
    "kind": "analysis",
    "provides": [
        "OUTCOME_COLUMN",
        "TARGET_COLUMN",
        "ADD_LEVEL",
        "OTHER_LEVELS",
        "DATA_VERSION",
        "build_target",
        "run_boruta",
        "boruta_params",
        "main",
    ],
    "uses": [
        "project_loader",
        "boruta_selection",
        "common.data_loading",
        "common.result_io",
    ],
    "seeded_from": {
        "boruta": {"template": "boruta", "version": "0.1"},
        "cache": {"template": "result-io", "version": "0.2"},
    },
    "description": (
        "BORUTA all-relevant selection (COMPUTE + cache): ADD (10, positive) vs pooled "
        "Other {PDD,PDCN,HCN} (30) on the delivered log2 40x2340 analysis subset. "
        "Builds a categorical ADD/Other target on Condition so classification is "
        "inferred, runs boruta_select ONCE on the whole matrix (shadow built-in null; "
        "no CV, no "
        "external null), and caches the BorutaResult via result_io.save_cached_result "
        "(analysis 'boruta') returning the ResultMeta for registration. Tuning knobs "
        "(max_iter/alpha/perc) enter provenance.params + the fingerprint. EXPLORATORY; "
        "small n; ADD==Kerchner cohort confound (finding 0001) is a hard bound the "
        "shadow null cannot resolve."
    ),
}

OUTCOME_COLUMN = "Condition"
TARGET_COLUMN = "add_vs_other"
ADD_LEVEL = "ADD"
OTHER_LEVELS: tuple[str, ...] = ("PDD", "PDCN", "HCN")
ADD_TARGET = "ADD"
OTHER_TARGET = "Other"

# Data-version stamp for the cache fingerprint (state/DATA_DESCRIPTION.md).
DATA_VERSION = "protein-data.csv sha256:90a0462eb03a02bd"

DEFAULT_SEED = 20260915
DEFAULT_CACHE_ROOT = Path("results")

# Boruta defaults (mirror the template); overridable on the command line.
DEFAULT_MAX_ITER = 100
DEFAULT_ALPHA = 0.05
DEFAULT_PERC = 100.0


def _assert_expected_split(dataset: Dataset) -> None:
    """Fail loud unless Condition is exactly 10 ADD + 30 Other (a data-change guard)."""
    conditions = dataset.metadata[OUTCOME_COLUMN].astype(str)
    seen = set(conditions.unique())
    expected = {ADD_LEVEL, *OTHER_LEVELS}
    if seen != expected:
        raise ValueError(
            f"Condition levels {sorted(seen)} != expected {sorted(expected)} "
            f"(controls should already be excluded upstream)."
        )
    n_add = int((conditions == ADD_LEVEL).sum())
    n_other = int(conditions.isin(OTHER_LEVELS).sum())
    if (n_add, n_other) != (10, 30):
        raise ValueError(
            f"ADD-vs-Other split is {n_add} ADD / {n_other} Other; expected 10 / 30."
        )
    logger.info("target: %d ADD (positive) vs %d Other (negative)", n_add, n_other)


def build_target(dataset: Dataset) -> Dataset:
    """Attach a categorical ``add_vs_other`` target (ADD vs Other) to the metadata.

    ADD -> "ADD" (positive, the 1 partition), {PDD,PDCN,HCN} -> "Other" (the 0
    partition). Stored as a pandas ``category`` dtype (NOT numeric 0/1) so
    ``boruta_select`` infers a **classification** task from the target dtype rather than
    mis-inferring regression from a low-cardinality numeric column.
    """
    conditions = dataset.metadata[OUTCOME_COLUMN].astype(str)
    target = conditions.where(conditions == ADD_LEVEL, OTHER_TARGET)
    target = target.where(conditions != ADD_LEVEL, ADD_TARGET)
    metadata = dataset.metadata.copy()
    metadata[TARGET_COLUMN] = pd.Categorical(
        target, categories=[OTHER_TARGET, ADD_TARGET]
    )
    return Dataset(
        abundances=dataset.abundances,
        feature_names=dataset.feature_names,
        feature_metadata=dataset.feature_metadata,
        metadata=metadata,
        scale=dataset.scale,
    )


def boruta_params(
    *, max_iter: int, alpha: float, perc: float, seed: int
) -> dict[str, object]:
    """The identity params for the cache fingerprint (every numbers-changing knob)."""
    return {
        "target": TARGET_COLUMN,
        "binarize": f"{ADD_LEVEL}=1 vs {{{','.join(OTHER_LEVELS)}}}=0",
        "positive_label": ADD_LEVEL,
        "negative_levels": sorted(OTHER_LEVELS),
        "task": "classification",
        "method": "boruta_all_relevant",
        "estimator": "random_forest_classifier(class_weight=balanced)",
        "n_estimators": "auto",
        "max_iter": max_iter,
        "alpha": alpha,
        "perc": perc,
        "two_step": True,
        "run_null": False,
        "shadow_null": "built-in (permuted shadow features)",
        "cv": "none (selection on whole matrix)",
        "seed": seed,
    }


def run_boruta(
    dataset: Dataset, *, max_iter: int, alpha: float, perc: float, seed: int
) -> object:
    """Run Boruta all-relevant selection once on the whole experimental matrix."""
    result = boruta_select(
        dataset,
        target=TARGET_COLUMN,
        task="classification",
        max_iter=max_iter,
        alpha=alpha,
        perc=perc,
        random_state=seed,
    )
    logger.info(
        "boruta: %d Confirmed / %d Tentative / %d Rejected over %d features "
        "(%d iterations, shadow_threshold=%.4f)",
        result.n_confirmed,
        result.n_tentative,
        result.n_rejected,
        result.n_features,
        result.n_iter,
        result.shadow_threshold,
    )
    return result


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("data/protein-data.csv"))
    parser.add_argument("--metadata", type=Path, default=Path("data/metadata.csv"))
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--max-iter", type=int, default=DEFAULT_MAX_ITER)
    parser.add_argument("--alpha", type=float, default=DEFAULT_ALPHA)
    parser.add_argument("--perc", type=float, default=DEFAULT_PERC)
    parser.add_argument("--cache-root", type=Path, default=DEFAULT_CACHE_ROOT)
    args = parser.parse_args(argv)

    dataset = load_analysis_subset(load_paired_dataset(args.data, args.metadata))
    logger.info(
        "loaded analysis subset: %d samples x %d features, scale=%s",
        dataset.abundances.shape[0],
        dataset.abundances.shape[1],
        dataset.scale,
    )
    _assert_expected_split(dataset)
    dataset = build_target(dataset)

    result = run_boruta(
        dataset,
        max_iter=args.max_iter,
        alpha=args.alpha,
        perc=args.perc,
        seed=args.seed,
    )

    params = boruta_params(
        max_iter=args.max_iter, alpha=args.alpha, perc=args.perc, seed=args.seed
    )
    tuned = (args.max_iter, args.alpha, args.perc) != (
        DEFAULT_MAX_ITER,
        DEFAULT_ALPHA,
        DEFAULT_PERC,
    )
    label = (
        f"add-vs-other Boruta all-relevant "
        f"(max_iter={args.max_iter}, alpha={args.alpha}, perc={args.perc}"
        f"{'; TUNED' if tuned else '; defaults'})"
    )
    meta: ResultMeta = save_cached_result(
        result,
        cache_root=args.cache_root,
        analysis="boruta",
        data_version=DATA_VERSION,
        params=params,
        seed=args.seed,
        label=label,
    )
    logger.info(
        "cached boruta result -> id=%s path=%s",
        meta.fingerprint,
        Path(args.cache_root) / "boruta" / meta.fingerprint,
    )
    logger.info("RESULT_META fingerprint=%s label=%s", meta.fingerprint, meta.label)


if __name__ == "__main__":
    main()
