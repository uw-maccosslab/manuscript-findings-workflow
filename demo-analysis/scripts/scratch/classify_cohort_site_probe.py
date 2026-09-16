"""COHORT CONFOUND PROBE (Task 2): can the proteome predict collection SITE among the
NON-ADD samples?

The ADD classifier (scripts/scratch/classify_add_vs_other.py) separates ADD from Other
at ~AUC 0.94-0.96, but ADD is 100% Kerchner cohort (finding 0001), so that separation is
inseparable from collection site. This probe asks the site question DIRECTLY, where it
is *not* aliased with disease: among the 30 NON-ADD samples, can the SAME leakage-safe
elastic-net classifier predict Kerchner vs Poston collection site?

Design:

  * Subset to the 30 non-ADD samples (Condition in {PDD, PDCN, HCN}).
  * Binary SITE label from ``Cohort``: **Kerchner (n=7, all HCN, positive) vs Poston
    (n=20, negative)**. The 3 samples from other sites (ADRC/PUC n=2, PUC only n=1) are
    **DROPPED** (LevelMap lists only Kerchner/Poston, so unlisted levels are excluded) —
    NOT folded into a "non-Kerchner" negative. Rationale: this keeps the negative class
    a single coherent collection site (Poston), so the contrast is a clean two-site
    Kerchner-vs-Poston probe rather than Kerchner-vs-heterogeneous-mix; the 3 dropped
    samples are two further near-singleton sites that would muddy "site" as a construct.
    Final: 7 positive vs 20 negative = 27 samples.
  * ``classify`` (elastic net), same leakage-safe nested CV as the ADD run: stratified
    5-fold x 5-repeat, in-fold StandardScaler, class_weight="balanced",
    generalization_target="individuals", run_null=True (this tiny-n probe NEEDS its own
    null), same seed.

INTERPRETATION (state BOTH directions; do not overclaim either way):

  * The matrix is ALREADY limma removeBatchEffect-corrected on Cohort, and the
    within-HCN Kerchner-vs-Poston mean delta was verified ~0 post-correction. So a
    **NULL-level** site AUC here means "the correction removed detectable within-non-ADD
    site structure" — it does NOT clear the ADD result, because ADD is single-cohort and
    its correction was EXTRAPOLATED, not estimated (finding 0004).
  * An **ABOVE-null** site AUC means residual site structure survives correction even
    where it was estimable — strong evidence the ADD separation is substantially site,
    not disease.

Also (cheap, informative): the OVERLAP between the ADD-classifier's top selected
features and this site probe's top selected features. High overlap = the ADD signal and
the site signal ride on the same proteins (damning for a disease interpretation).

Everything EXPLORATORY; tiny n (7 positives); imbalanced; a bounded probe, not a powered
analysis.

Run::

    ./.venv/bin/python scripts/scratch/classify_cohort_site_probe.py
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from classification import ClassificationResult, LevelMap, classify
from classify_add_vs_other import DATA_VERSION, _parse_gene
from common.data_loading import Dataset
from common.result_io import (
    ResultMeta,
    load_cached_result,
    result_fingerprint,
    save_cached_result,
)
from project_loader import load_analysis_subset, load_paired_dataset
from scipy.stats import hypergeom

logger = logging.getLogger(__name__)

__script_meta__: dict[str, object] = {
    "task": "classify-cohort-site-probe",
    "kind": "analysis",
    "provides": [
        "SITE_COLUMN",
        "POSITIVE_SITE",
        "NEGATIVE_SITE",
        "N_SPLITS",
        "N_REPEATS",
        "subset_non_add",
        "run_site_probe",
        "site_params",
        "feature_overlap",
        "null_read",
        "main",
    ],
    "uses": [
        "project_loader",
        "classification",
        "classify_add_vs_other",
        "common.data_loading",
        "common.result_io",
    ],
    "seeded_from": {
        "elastic_net": {"template": "classification", "version": "0.3"},
        "cache": {"template": "result-io", "version": "0.2"},
    },
    "description": (
        "Cohort/site confound probe: among the 30 NON-ADD samples, elastic-net "
        "classification of collection site Kerchner (7, positive) vs Poston (20, "
        "negative), dropping the 3 other-site samples (ADRC/PUC + PUC only). Same "
        "leakage-safe nested CV as the ADD run (stratified 5-fold x 5-repeat, in-fold "
        "scaling, class_weight balanced, target=individuals, run_null=True, seed "
        "20260915). Caches the result (new fingerprint) and reports the site-probe AUC "
        "vs its own null plus the OVERLAP of its top selected features with the ADD "
        "classifier's top selected features. Interpretation is confound-bounded BOTH "
        "ways (null-level != ADD cleared, because ADD's on-Cohort correction was "
        "extrapolated; above-null = residual site structure survives correction). "
        "EXPLORATORY, tiny n (7 positives), imbalanced."
    ),
}

SITE_COLUMN = "Cohort"
CONDITION_COLUMN = "Condition"
ADD_LEVEL = "ADD"
POSITIVE_SITE = "Kerchner"
NEGATIVE_SITE = "Poston"
DROPPED_SITES: tuple[str, ...] = ("ADRC/PUC", "PUC only")

N_SPLITS = 5
N_REPEATS = 5
STABILITY_REPEATS = 10
N_PERMUTATIONS = 1000
NULL_REPEATS = 3

#: Match the ADD-side EN null pass (max_iter=20000) so both elastic-net models are on
#: equal, converged footing — the selected-feature overlap compares the two, so an
#: asymmetric convergence setting would confound the comparison. The default 5000 leaves
#: saga unconverged (ConvergenceWarning) on this p>>n matrix.
EN_MAX_ITER = 20000

#: The current ADD-classifier EN caches to compare selected-feature overlap against.
#: The NULL-pass EN (Task 1) is preferred (superseding, tuning-fixed); the first-pass EN
#: is the fallback if the null pass has not run yet.
ADD_EN_FIRST_PASS_FP = "cbd809ac6168"


def subset_non_add(dataset: Dataset) -> Dataset:
    """Return a Dataset of the non-ADD samples only (features/scale untouched)."""
    conditions = dataset.metadata[CONDITION_COLUMN].astype(str).to_numpy()
    mask = conditions != ADD_LEVEL
    n_non_add = int(mask.sum())
    if n_non_add != 30:
        raise ValueError(f"Expected 30 non-ADD samples; got {n_non_add}.")
    subset = Dataset(
        abundances=dataset.abundances[mask, :],
        feature_names=dataset.feature_names,
        feature_metadata=dataset.feature_metadata,
        metadata=dataset.metadata.iloc[mask].copy(),
        scale=dataset.scale,
    )
    counts = subset.metadata[SITE_COLUMN].astype(str).value_counts().to_dict()
    logger.info("non-ADD subset cohort counts: %s", counts)
    if counts.get(POSITIVE_SITE) != 7 or counts.get(NEGATIVE_SITE) != 20:
        raise ValueError(
            f"Expected Kerchner=7 / Poston=20 among non-ADD; got {counts}."
        )
    return subset


def site_params(seed: int) -> dict[str, object]:
    """Cache fingerprint params — every knob that changes the numbers."""
    return {
        "sample_subset": "non_ADD (Condition in {PDD,PDCN,HCN})",
        "outcome": SITE_COLUMN,
        "binarize": "LevelMap",
        "positive_label": POSITIVE_SITE,
        "negative_label": NEGATIVE_SITE,
        "dropped_levels": list(DROPPED_SITES),
        "method": "elastic_net_logistic",
        "generalization_target": "individuals",
        "cv_scheme": "stratified_kfold",
        "n_splits": N_SPLITS,
        "n_repeats": N_REPEATS,
        "grouped": False,
        "stability_repeats": STABILITY_REPEATS,
        "tuning_metric": "roc_auc",
        "select": "best",
        "class_weight": "balanced",
        "max_iter": EN_MAX_ITER,
        "feature_list": None,
        "run_null": True,
        "n_permutations": N_PERMUTATIONS,
        "null_repeats": NULL_REPEATS,
        "null_permutation": "units",
        "seed": seed,
    }


def run_site_probe(
    subset: Dataset, *, seed: int
) -> tuple[ClassificationResult, list[warnings.WarningMessage]]:
    """Elastic-net Kerchner-vs-Poston site classifier with its own null."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = classify(
            subset,
            outcome=SITE_COLUMN,
            binarize=LevelMap(positive=(POSITIVE_SITE,), negative=(NEGATIVE_SITE,)),
            generalization_target="individuals",
            n_splits=N_SPLITS,
            n_repeats=N_REPEATS,
            stability_repeats=STABILITY_REPEATS,
            run_null=True,
            n_permutations=N_PERMUTATIONS,
            null_repeats=NULL_REPEATS,
            null_permutation="units",
            max_iter=EN_MAX_ITER,
            random_state=seed,
        )
    return result, list(caught)


def null_read(result: ClassificationResult) -> dict[str, Any]:
    assert result.null_aucs is not None and result.observed_auc is not None
    assert result.null_p is not None
    nulls = result.null_aucs
    null_mean = float(np.mean(nulls))
    lo, hi = (float(x) for x in np.percentile(nulls, [2.5, 97.5]))
    observed = float(result.observed_auc)
    collapses = bool(observed <= hi)
    return {
        "nested_cv_auc_tune_in_fold": round(result.cv_auc, 4),
        "nested_cv_auc_sd_over_folds": round(result.cv_auc_sd, 4),
        "cv_balanced_accuracy": round(result.cv_balanced_accuracy, 4),
        "cv_pr_auc": round(result.cv_average_precision, 4),
        "observed_auc_fixed_hyperparam": round(observed, 4),
        "null_auc_mean": round(null_mean, 4),
        "null_auc_95pct_interval": [round(lo, 4), round(hi, 4)],
        "null_auc_max": round(float(np.max(nulls)), 4),
        "observed_minus_null_mean": round(observed - null_mean, 4),
        "empirical_permutation_p": round(float(result.null_p), 5),
        "n_permutations": int(nulls.size),
        "collapses_to_null": collapses,
        "verdict": (
            "SITE AUC COLLAPSES to its null — the on-Cohort correction removed "
            "detectable within-non-ADD site structure. This does NOT clear the ADD "
            "result: ADD is single-cohort, so its correction was EXTRAPOLATED, not "
            "estimated (finding 0004)."
            if collapses
            else "SITE AUC STAYS ABOVE its null — residual site structure survives the "
            "on-Cohort correction even where it was estimable. Strong evidence the ADD "
            "separation is substantially SITE, not disease."
        ),
    }


def _load_add_en_coeffs(cache_root: Path) -> tuple[pd.DataFrame, str]:
    """Load the ADD-classifier EN selected-feature table (null pass preferred).

    Scans results/classification/ for the NULL-pass EN cache (run_null=True in its
    meta.json); falls back to the known first-pass fingerprint. Returns
    (coefficients, fingerprint).
    """
    en_dir = cache_root / "classification"
    chosen: str | None = None
    for sub in sorted(p for p in en_dir.iterdir() if p.is_dir()):
        meta_path = sub / "meta.json"
        if not meta_path.is_file():
            continue
        meta = json.loads(meta_path.read_text())
        params = meta.get("params", {})
        if params.get("run_null") is True:
            chosen = sub.name
            break
    if chosen is None:
        chosen = ADD_EN_FIRST_PASS_FP
        logger.warning(
            "no run_null EN cache found; falling back to first-pass EN %s for overlap.",
            chosen,
        )
    result = load_cached_result(
        ClassificationResult,
        cache_root=cache_root,
        analysis="classification",
        fingerprint=chosen,
    )
    return result.coefficients, chosen


def feature_overlap(
    site_coeffs: pd.DataFrame,
    add_coeffs: pd.DataFrame,
    add_fp: str,
    *,
    n_features_total: int,
    top_ns: tuple[int, ...] = (25, 50),
) -> dict[str, Any]:
    """Overlap of the two classifiers' top selected features (by |coef|) + full-set
    Jaccard + a descriptive hypergeometric enrichment p.

    Both tables list only features NON-ZERO in the respective final EN model (the
    selected set). Overlap on the whole selected sets and on top-N by |coef| is
    reported; the hypergeometric p treats each top-N as a draw of N from the
    ``n_features_total`` background (descriptive — the two selections are not
    independent draws, so read it as a rough enrichment gauge, not an inferential test).
    """
    site_selected = [str(f) for f in site_coeffs["feature"].tolist()]
    add_selected = [str(f) for f in add_coeffs["feature"].tolist()]
    site_set = set(site_selected)
    add_set = set(add_selected)

    full_inter = site_set & add_set
    union = site_set | add_set
    jaccard = len(full_inter) / len(union) if union else 0.0

    per_top: list[dict[str, Any]] = []
    for n in top_ns:
        s_top = [str(f) for f in site_coeffs["feature"].head(n).tolist()]
        a_top = [str(f) for f in add_coeffs["feature"].head(n).tolist()]
        s_top_set, a_top_set = set(s_top), set(a_top)
        inter = s_top_set & a_top_set
        k = len(inter)
        # P(X >= k) drawing len(s_top_set) from n_features_total with len(a_top_set)
        # flagged "success" (the ADD top-N). Descriptive.
        p_enrich = float(
            hypergeom.sf(k - 1, n_features_total, len(a_top_set), len(s_top_set))
        )
        per_top.append(
            {
                "top_n": n,
                "site_top_n": len(s_top_set),
                "add_top_n": len(a_top_set),
                "n_overlap": k,
                "overlap_genes": sorted(_parse_gene(f) for f in inter),
                "hypergeometric_p_enrichment_descriptive": round(p_enrich, 6),
            }
        )

    return {
        "add_en_cache_fingerprint": add_fp,
        "n_features_total_background": n_features_total,
        "n_site_selected_nonzero": len(site_set),
        "n_add_selected_nonzero": len(add_set),
        "full_selected_set_overlap": len(full_inter),
        "full_selected_set_jaccard": round(jaccard, 4),
        "per_top_n": per_top,
        "note": (
            "High overlap => the ADD signal and the site signal ride on the same "
            "proteins (damning for a disease interpretation). Low overlap is NOT "
            "exoneration: EN keeps one of a correlated cluster, so the two minimal-"
            "optimal sets can pick different members of the same underlying structure."
        ),
    }


def _site_stability_top(result: ClassificationResult, n: int) -> list[dict[str, Any]]:
    cols = [
        "feature",
        "coef",
        "selection_frequency",
        "sign_consistency",
        "coef_median",
        "coef_q25",
        "coef_q75",
    ]
    return [
        {
            "feature": str(r["feature"]),
            "gene": _parse_gene(str(r["feature"])),
            "coef": round(float(r["coef"]), 4),
            "coef_sign_toward": (
                POSITIVE_SITE if float(r["coef"]) > 0 else NEGATIVE_SITE
            ),
            "selection_frequency": round(float(r["selection_frequency"]), 3),
            "sign_consistency": round(float(r["sign_consistency"]), 3),
        }
        for _, r in result.coefficients[cols].head(n).iterrows()
    ]


def _warn_lines(caught: list[warnings.WarningMessage]) -> list[str]:
    seen: dict[str, int] = {}
    for w in caught:
        key = f"{w.category.__name__}: {w.message}"
        seen[key] = seen.get(key, 0) + 1
    return [f"{k} (x{n})" if n > 1 else k for k, n in seen.items()]


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("data/protein-data.csv"))
    parser.add_argument("--metadata", type=Path, default=Path("data/metadata.csv"))
    parser.add_argument("--seed", type=int, default=20260915)
    parser.add_argument("--cache-root", type=Path, default=Path("results"))
    parser.add_argument(
        "--summary-dir", type=Path, default=Path("results/classification_cohort_site")
    )
    parser.add_argument("--top-n", type=int, default=25)
    args = parser.parse_args(argv)

    dataset = load_analysis_subset(load_paired_dataset(args.data, args.metadata))
    subset = subset_non_add(dataset)

    params = site_params(args.seed)
    fp = result_fingerprint(
        analysis="classification_cohort_site",
        data_version=DATA_VERSION,
        params=params,
        seed=args.seed,
    )
    logger.info("target site-probe fingerprint: %s", fp)
    if (args.cache_root / "classification_cohort_site" / fp).exists():
        logger.warning("site-probe result %s already cached — reuse it.", fp)

    logger.info("running site probe (Kerchner vs Poston, non-ADD, run_null=True)...")
    result, caught = run_site_probe(subset, seed=args.seed)

    meta: ResultMeta = save_cached_result(
        result,
        cache_root=args.cache_root,
        analysis="classification_cohort_site",
        data_version=DATA_VERSION,
        params=params,
        seed=args.seed,
        label="cohort-site probe Kerchner-vs-Poston among non-ADD (EN, run_null=True)",
    )

    add_coeffs, add_fp = _load_add_en_coeffs(args.cache_root)
    overlap = feature_overlap(
        result.coefficients,
        add_coeffs,
        add_fp,
        n_features_total=int(dataset.abundances.shape[1]),
    )

    summary: dict[str, Any] = {
        "analysis": "COHORT/SITE confound probe (non-ADD): Kerchner vs Poston",
        "data_version": DATA_VERSION,
        "seed": args.seed,
        "cache_fingerprint": meta.fingerprint,
        "cache_path": str(
            args.cache_root / "classification_cohort_site" / meta.fingerprint
        ),
        "design": (
            "non-ADD subset (30) -> Cohort LevelMap Kerchner (7, positive) vs Poston "
            "(20, negative); 3 other-site samples (ADRC/PUC=2, PUC only=1) DROPPED; "
            "n=27. Elastic net, stratified 5-fold x 5-repeat, in-fold scaling, "
            "class_weight balanced, target=individuals."
        ),
        "n_samples": result.n_samples,
        "n_positive_kerchner": result.n_positive,
        "n_negative_poston": result.n_negative,
        "n_dropped_other_sites": result.n_dropped_unassigned,
        "best_c": result.best_c,
        "best_l1_ratio": result.best_l1_ratio,
        "n_features_selected_nonzero": len(result.coefficients),
        "null_read": null_read(result),
        "feature_overlap_with_add_classifier": overlap,
        "site_probe_stability_top": _site_stability_top(result, args.top_n),
        "warnings": _warn_lines(caught),
        "caveats": [
            "EXPLORATORY: tiny n (7 Kerchner positives), imbalanced; a bounded probe, "
            "not a powered analysis. Wide CIs.",
            "The matrix is ALREADY batch-corrected on Cohort (finding 0004). A "
            "null-level site AUC means the correction removed detectable within-nonADD "
            "site structure — it does NOT clear the ADD result (ADD single-cohort => "
            "correction extrapolated, not estimated).",
            "An above-null site AUC means residual site structure survives correction "
            "where it was estimable — evidence the ADD separation is substantially "
            "site.",
            "Low feature overlap is NOT exoneration (EN keeps one of a correlated "
            "cluster; the two minimal-optimal sets can pick different members of the "
            "same structure).",
        ],
    }

    args.summary_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.summary_dir / "site_probe.json"
    with summary_path.open("w") as fh:
        json.dump(summary, fh, indent=2)
    logger.info("wrote summary -> %s", summary_path)
    logger.info("site-probe null read: %s", summary["null_read"]["verdict"])
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
