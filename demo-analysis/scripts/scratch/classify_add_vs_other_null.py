"""PRIMARY classification: ADD vs Other (pooled) — LABEL-SHUFFLE NULL pass (Task 1).

The required second pass for the two first-pass ADD-vs-Other linear classifiers
(``scripts/scratch/classify_add_vs_other.py``): re-run BOTH models with
``run_null=True`` on the SAME outcome / binarize / CV design / seed as the first pass,
applying the two tuning fixes the first pass flagged, and cache each as a NEW result
(new fingerprint) that SUPERSEDES the corresponding first-pass cache.

The label-shuffle null is:

  * the **leakage detector** — if nested performance does NOT collapse toward the null
    when the target is permuted, there is either leakage or no real signal;
  * the gate that lifts the finding from the first-pass ``exploratory`` cap toward
    ``validated`` eligibility (the coefficients/weights become licensed *once above
    null*).

It is NOT a fix for the cohort confound: ADD is 100% Kerchner among these 40 samples
(finding 0001), so a unit-level label shuffle permutes the very partition that is
simultaneously "ADD" and "Kerchner". A clean above-null result rules out leakage / pure
overfitting; it does NOT license "disease, not site."

Tuning fixes applied (first-pass flags — reported, not silently worked around):

  * **Linear SVM** — the default broad C grid let 13/25 outer folds tune ``C`` *below*
    the ~1e-3 plateau (``TuningNoiseWarning``), deflating the nested AUC by tuning noise
    rather than data. Narrow ``c_grid`` to ``(1e-4 .. 1e-1)`` (half-decade steps) so the
    grid brackets the knee and the low-C tail the inner folds cannot resolve is gone.
  * **Elastic net** — the saga solver did not fully converge at ``max_iter=5000`` on
    the p>>n matrix (``ConvergenceWarning``), so the reported coefficient magnitudes
    were not yet trustworthy. Raise ``max_iter`` to 20000 to clear it.

Both fixes change the cache fingerprint (``c_grid`` / ``max_iter`` are in ``params``),
as does ``run_null`` flipping True — so these are genuinely new results, kept alongside
the first-pass ones (keep-all; the findings-manager marks first-pass rows superseded).

Run::

    ./.venv/bin/python scripts/scratch/classify_add_vs_other_null.py
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import warnings
from pathlib import Path
from typing import Any, Literal

import numpy as np

# Make `scripts/` and `scripts/scratch/` importable regardless of invocation.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from classification import ClassificationResult, classify
from classification_svm import SVMClassificationResult, classify_svm
from classify_add_vs_other import (
    DATA_VERSION,
    N_REPEATS,
    N_SPLITS,
    OUTCOME_COLUMN,
    STABILITY_REPEATS,
    SVM_TOP_K,
    _en_params,
    _svm_params,
    build_binarize,
    build_binarize_svm,
)
from common.data_loading import Dataset
from common.result_io import ResultMeta, result_fingerprint, save_cached_result
from project_loader import load_analysis_subset, load_paired_dataset

logger = logging.getLogger(__name__)

__script_meta__: dict[str, object] = {
    "task": "classify-add-vs-other-null",
    "kind": "analysis",
    "provides": [
        "EN_NULL_MAX_ITER",
        "SVM_NULL_C_GRID",
        "N_PERMUTATIONS",
        "NULL_REPEATS",
        "en_null_params",
        "svm_null_params",
        "run_elastic_net_null",
        "run_linear_svm_null",
        "null_summary",
        "main",
    ],
    "uses": [
        "project_loader",
        "classification",
        "classification_svm",
        "classify_add_vs_other",
        "common.data_loading",
        "common.result_io",
    ],
    "seeded_from": {
        "elastic_net": {"template": "classification", "version": "0.3"},
        "linear_svm": {"template": "classification-svm", "version": "0.3"},
        "cache": {"template": "result-io", "version": "0.2"},
        "first_pass": "scripts/scratch/classify_add_vs_other.py",
    },
    "description": (
        "Label-shuffle NULL pass (run_null=True) for the ADD (10, positive) vs pooled "
        "Other {PDD,PDCN,HCN} (30) elastic-net + linear-SVM classifiers, SAME "
        "outcome/binarize/CV/seed (20260915) as the first pass, with the two flagged "
        "tuning fixes: SVM c_grid narrowed to (1e-4..1e-1) to bracket the ~1e-3 knee "
        "(kills the 13/25 sub-plateau folds), EN max_iter raised to 20000 to clear the "
        "saga ConvergenceWarnings. Caches each as a new fingerprint (SUPERSEDES "
        "first-pass cbd809ac6168 EN / cba84421b85e SVM); reports recomputed observed "
        "nested AUC, the null AUC distribution (mean, 95%), the empirical permutation "
        "p, and whether observed collapses to the null. Above-null rules out "
        "leakage/overfitting only; it does NOT break the ADD==Kerchner single-cohort "
        "confound (the null permutes that shared partition). Exploratory; small n."
    ),
}

# --- Tuning fixes (the two the first pass flagged) --------------------------------- #

#: Elastic net: raise the saga iteration cap to clear ConvergenceWarnings on the p>>n
#: matrix so the reported coefficient magnitudes are trustworthy (first-pass used 5000).
EN_NULL_MAX_ITER = 20000

#: Linear SVM: narrow the C grid to bracket the ~1e-3 plateau knee (half-decade steps).
#: The default broad grid reached down to 1e-5, and 13/25 outer folds tuned below the
#: plateau (TuningNoiseWarning), deflating the nested AUC. This grid keeps two points
#: below the knee (1e-4, 3e-4) and up to 1e-1 above it — strictly increasing (required).
SVM_NULL_C_GRID: tuple[float, ...] = (1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1)

# --- Null knobs (template defaults; recorded for provenance) ----------------------- #
N_PERMUTATIONS = 1000
NULL_REPEATS = 3
# groups=None -> resolves to a plain row ("samples") shuffle
NULL_PERMUTATION: Literal["units", "within_units"] = "units"


def en_null_params(seed: int) -> dict[str, object]:
    """EN cache params = first-pass EN params + the null + the max_iter tuning fix.

    ``max_iter`` is added to the fingerprint because it changes the coefficient
    magnitudes (it was not a knob in the first-pass params dict; the first pass ran the
    template default 5000).
    """
    params = _en_params(seed)
    params["run_null"] = True
    params["n_permutations"] = N_PERMUTATIONS
    params["null_repeats"] = NULL_REPEATS
    params["max_iter"] = EN_NULL_MAX_ITER
    return params


def svm_null_params(seed: int) -> dict[str, object]:
    """SVM cache params = first-pass SVM params + the null + the narrowed c_grid."""
    params = _svm_params(seed)
    params["run_null"] = True
    params["n_permutations"] = N_PERMUTATIONS
    params["null_repeats"] = NULL_REPEATS
    params["c_grid"] = list(SVM_NULL_C_GRID)
    return params


def run_elastic_net_null(
    dataset: Dataset, *, seed: int
) -> tuple[ClassificationResult, list[warnings.WarningMessage]]:
    """Elastic-net logistic, NULL pass (run_null=True, max_iter=20000)."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = classify(
            dataset,
            outcome=OUTCOME_COLUMN,
            binarize=build_binarize(),
            generalization_target="individuals",
            n_splits=N_SPLITS,
            n_repeats=N_REPEATS,
            stability_repeats=STABILITY_REPEATS,
            run_null=True,
            n_permutations=N_PERMUTATIONS,
            null_repeats=NULL_REPEATS,
            null_permutation=NULL_PERMUTATION,
            max_iter=EN_NULL_MAX_ITER,
            random_state=seed,
        )
    return result, list(caught)


def run_linear_svm_null(
    dataset: Dataset, *, seed: int
) -> tuple[SVMClassificationResult, list[warnings.WarningMessage]]:
    """Linear SVM, NULL pass (run_null=True, narrowed c_grid); same seed as EN."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = classify_svm(
            dataset,
            outcome=OUTCOME_COLUMN,
            binarize=build_binarize_svm(),
            generalization_target="individuals",
            c_grid=SVM_NULL_C_GRID,
            top_k=SVM_TOP_K,
            n_splits=N_SPLITS,
            n_repeats=N_REPEATS,
            stability_repeats=STABILITY_REPEATS,
            run_null=True,
            n_permutations=N_PERMUTATIONS,
            null_repeats=NULL_REPEATS,
            null_permutation=NULL_PERMUTATION,
            random_state=seed,
        )
    return result, list(caught)


def _null_stats(
    observed: float, nulls: np.ndarray, null_p: float, nested_auc: float
) -> dict[str, Any]:
    """The null read: observed (fixed-hp) vs the permutation distribution + verdict.

    ``observed`` is the fixed-hyperparameter procedure's AUC (what ``null_p`` belongs
    to); ``nested_auc`` is the tune-in-fold nested estimate reported as headline
    performance. Both are given so the reader does not conflate them.
    """
    null_mean = float(np.mean(nulls))
    lo, hi = (float(x) for x in np.percentile(nulls, [2.5, 97.5]))
    null_max = float(np.max(nulls))
    margin = observed - null_mean
    collapses = bool(observed <= hi)  # observed inside the null's upper 95% => collapse
    return {
        "nested_cv_auc_tune_in_fold": round(nested_auc, 4),
        "observed_auc_fixed_hyperparam": round(observed, 4),
        "null_auc_mean": round(null_mean, 4),
        "null_auc_95pct_interval": [round(lo, 4), round(hi, 4)],
        "null_auc_max": round(null_max, 4),
        "observed_minus_null_mean": round(margin, 4),
        "empirical_permutation_p": round(null_p, 5),
        "n_permutations": int(nulls.size),
        "collapses_to_null": collapses,
        "verdict": (
            "COLLAPSES to the null (observed inside the null 95% band) — leakage or no "
            "real signal"
            if collapses
            else "STAYS CLEARLY ABOVE the null — rules out leakage / pure overfitting "
            "(does NOT resolve the ADD==Kerchner cohort confound)"
        ),
    }


def _warn_lines(caught: list[warnings.WarningMessage]) -> list[str]:
    """Compact, de-duplicated 'Category: message' lines for the captured warnings."""
    seen: dict[str, int] = {}
    for w in caught:
        key = f"{w.category.__name__}: {w.message}"
        seen[key] = seen.get(key, 0) + 1
    return [f"{k} (x{n})" if n > 1 else k for k, n in seen.items()]


def null_summary(
    en_result: ClassificationResult,
    svm_result: SVMClassificationResult,
    en_meta: ResultMeta,
    svm_meta: ResultMeta,
    en_warnings: list[warnings.WarningMessage],
    svm_warnings: list[warnings.WarningMessage],
    *,
    seed: int,
    cache_root: Path,
) -> dict[str, Any]:
    assert en_result.null_aucs is not None and en_result.observed_auc is not None
    assert en_result.null_p is not None
    assert svm_result.null_aucs is not None and svm_result.observed_auc is not None
    assert svm_result.null_p is not None
    return {
        "analysis": "PRIMARY classification ADD vs Other — NULL pass (run_null=True)",
        "data_version": DATA_VERSION,
        "seed": seed,
        "outcome": "ADD (10, positive) vs Other {PDD,PDCN,HCN} (30, negative)",
        "cv": (
            f"stratified {N_SPLITS}-fold x {N_REPEATS} repeats at the individual (row) "
            "level; no groups (subjects are singletons)"
        ),
        "null_permutation": "units (groups=None -> plain row/'samples' shuffle)",
        "supersedes": {
            "elastic_net_first_pass": "cbd809ac6168",
            "linear_svm_first_pass": "cba84421b85e",
        },
        "elastic_net": {
            "cache_fingerprint": en_meta.fingerprint,
            "cache_path": str(cache_root / "classification" / en_meta.fingerprint),
            "max_iter": EN_NULL_MAX_ITER,
            "best_c": en_result.best_c,
            "best_l1_ratio": en_result.best_l1_ratio,
            "n_features_selected_nonzero": len(en_result.coefficients),
            "cv_balanced_accuracy": round(en_result.cv_balanced_accuracy, 4),
            "cv_pr_auc": round(en_result.cv_average_precision, 4),
            "null_read": _null_stats(
                en_result.observed_auc,
                en_result.null_aucs,
                en_result.null_p,
                en_result.cv_auc,
            ),
            "warnings": _warn_lines(en_warnings),
        },
        "linear_svm": {
            "cache_fingerprint": svm_meta.fingerprint,
            "cache_path": str(
                cache_root / "classification_svm" / svm_meta.fingerprint
            ),
            "c_grid": list(SVM_NULL_C_GRID),
            "best_c": svm_result.best_c,
            "plateau_start_c": svm_result.plateau_start_c,
            "n_folds_sub_plateau": svm_result.n_folds_sub_plateau,
            "n_support_vectors": svm_result.n_support_vectors,
            "cv_balanced_accuracy": round(svm_result.cv_balanced_accuracy, 4),
            "cv_pr_auc": round(svm_result.cv_average_precision, 4),
            "null_read": _null_stats(
                svm_result.observed_auc,
                svm_result.null_aucs,
                svm_result.null_p,
                svm_result.cv_auc,
            ),
            "warnings": _warn_lines(svm_warnings),
        },
        "caveats": [
            "The null rules out leakage / pure overfitting ONLY. It does NOT resolve "
            "the cohort confound: ADD is 100% Kerchner among these 40 samples (finding "
            "0001), so the unit-level label shuffle permutes the SAME partition that "
            "is both ADD and Kerchner. Above-null != disease-not-site.",
            "Matrix is already limma removeBatchEffect-corrected on Cohort (finding "
            "0004); ADD is single-cohort so its correction was EXTRAPOLATED, not "
            "estimated.",
            "Small n (10 positives): wide CV CIs; EXPLORATORY.",
        ],
    }


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("data/protein-data.csv"))
    parser.add_argument("--metadata", type=Path, default=Path("data/metadata.csv"))
    parser.add_argument("--seed", type=int, default=20260915)
    parser.add_argument("--cache-root", type=Path, default=Path("results"))
    parser.add_argument(
        "--summary-dir", type=Path, default=Path("results/classification_add_vs_other")
    )
    args = parser.parse_args(argv)

    dataset = load_analysis_subset(load_paired_dataset(args.data, args.metadata))
    logger.info(
        "loaded analysis subset: %d samples x %d features, scale=%s",
        dataset.abundances.shape[0],
        dataset.abundances.shape[1],
        dataset.scale,
    )

    en_params = en_null_params(args.seed)
    svm_params = svm_null_params(args.seed)
    en_fp = result_fingerprint(
        analysis="classification",
        data_version=DATA_VERSION,
        params=en_params,
        seed=args.seed,
    )
    svm_fp = result_fingerprint(
        analysis="classification_svm",
        data_version=DATA_VERSION,
        params=svm_params,
        seed=args.seed,
    )
    logger.info("target fingerprints: EN=%s SVM=%s", en_fp, svm_fp)
    if (args.cache_root / "classification" / en_fp).exists():
        logger.warning(
            "EN null result %s already cached — reuse it (not recomputing).", en_fp
        )
    if (args.cache_root / "classification_svm" / svm_fp).exists():
        logger.warning(
            "SVM null result %s already cached — reuse it (not recomputing).", svm_fp
        )

    logger.info("running EN null pass (max_iter=%d)...", EN_NULL_MAX_ITER)
    en_result, en_warnings = run_elastic_net_null(dataset, seed=args.seed)
    logger.info("running SVM null pass (narrowed c_grid)...")
    svm_result, svm_warnings = run_linear_svm_null(dataset, seed=args.seed)

    en_meta = save_cached_result(
        en_result,
        cache_root=args.cache_root,
        analysis="classification",
        data_version=DATA_VERSION,
        params=en_params,
        seed=args.seed,
        label="add-vs-other elastic-net NULL pass (run_null=True, max_iter=20000)",
    )
    svm_meta = save_cached_result(
        svm_result,
        cache_root=args.cache_root,
        analysis="classification_svm",
        data_version=DATA_VERSION,
        params=svm_params,
        seed=args.seed,
        label="add-vs-other linear-SVM NULL pass (run_null=True, c_grid 1e-4..1e-1)",
    )

    summary = null_summary(
        en_result,
        svm_result,
        en_meta,
        svm_meta,
        en_warnings,
        svm_warnings,
        seed=args.seed,
        cache_root=args.cache_root,
    )
    args.summary_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.summary_dir / "null_pass.json"
    with summary_path.open("w") as fh:
        json.dump(summary, fh, indent=2)
    logger.info("wrote summary -> %s", summary_path)
    logger.info("EN  null read: %s", summary["elastic_net"]["null_read"]["verdict"])
    logger.info("SVM null read: %s", summary["linear_svm"]["null_read"]["verdict"])
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
