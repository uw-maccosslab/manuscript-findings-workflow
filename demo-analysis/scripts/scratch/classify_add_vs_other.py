"""PRIMARY classification: ADD vs Other (pooled), plasma-EV study — FIRST PASS.

Call site for the project copies of the two leakage-safe *linear* classifier lib/
templates:

  * ``scripts/scratch/classification.py``     — elastic-net logistic (``classify``),
    seeded from ``lib/analysis/classification`` @ 0.3;
  * ``scripts/scratch/classification_svm.py`` — soft-margin linear SVM
    (``classify_svm``), seeded from ``lib/analysis/classification-svm`` @ 0.3.

Both are run on the SAME seed and the SAME analyzed samples so their outer folds are
identical (the templates build the outer CV byte-identically) and the two nested-CV
AUCs can be compared **paired per fold** (conventions/statistics.md, "Choosing between
the two linear classifiers"). Neither runs the label-shuffle null (``run_null=False``):
this is the FIRST PASS — results first, the null is the immediate next step.

What it does:

  1. loads the 40-sample x 2340-protein log2 analysis subset via ``project_loader``
     (controls already excluded; yeast spike-in + two floored-sentinel features dropped
     upstream — state/DATA_DESCRIPTION.md);
  2. builds the binary outcome **ADD (n=10, positive) vs Other {PDD,PDCN,HCN} (n=30,
     negative)** via a ``LevelMap`` on ``Condition``;
  3. runs ``classify`` (elastic net) and ``classify_svm`` (linear SVM), both
     ``run_null=False``, stratified k-fold at the INDIVIDUAL (row) level (no groups: the
     40 subjects are singletons, so a held-out row IS a held-out individual —
     Stage-4 checkpoint), same ``random_state``;
  4. CACHES each result via ``result_io.save_cached_result`` under ``results/``
     (analysis ``classification`` / ``classification_svm``), returning the two
     ``ResultMeta`` for registration in ``results/manifest.md``;
  5. computes the PAIRED per-fold AUC comparison (EN vs SVM), the imbalance-appropriate
     metrics (ROC-AUC + per-fold AUCs, balanced accuracy, PR-AUC, MCC), the coefficient/
     weight stability reads, and the DE cross-checks; writes a JSON summary to
     ``results/classification_add_vs_other/first_pass_comparison.json``.

CRITICAL CAVEATS baked into the reporting (state/DATA_DESCRIPTION.md; findings 0001 +
0004):

  * **EXPLORATORY, capped.** Performance is NOT yet tested against a label-shuffle
    null, so the finding is capped ``exploratory`` and every coefficient/weight is
    flagged "not tested against a null" until the null pass runs (``run_null=True``).
    The null is also the leakage detector.
  * **COHORT CONFOUND (hard).** ADD is 100% Kerchner cohort (0001). Among these 40
    samples, "ADD" and "Kerchner" are the SAME partition, so any above-null performance
    CANNOT be attributed to disease vs collection-site — the model may be detecting
    Kerchner as much as ADD. The label-shuffle null does NOT resolve this (it permutes
    the very partition that is both ADD and Kerchner).
  * **Small n (10 positives).** Wide CV CIs; treat the point estimates cautiously.

Non-interactive, parameterized (argparse), seeded. Fail-loud at every boundary.

Run::

    ./.venv/bin/python scripts/scratch/classify_add_vs_other.py
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

# Make `scripts/` and `scripts/scratch/` importable regardless of invocation (direct
# script / -m / pytest), mirroring project_loader.py / de_add_vs_other.py.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from classification import (
    DEFAULT_C_GRID as EN_C_GRID,
)
from classification import (
    DEFAULT_L1_RATIOS as EN_L1_RATIOS,
)
from classification import (
    ClassificationResult,
    LevelMap,
    classify,
)
from classification_svm import (
    DEFAULT_C_GRID as SVM_C_GRID,
)
from classification_svm import (
    LevelMap as SvmLevelMap,
)
from classification_svm import (
    SVMClassificationResult,
    classify_svm,
)
from common.data_loading import Dataset
from common.result_io import ResultMeta, save_cached_result
from project_loader import load_analysis_subset, load_paired_dataset
from sklearn.metrics import matthews_corrcoef, roc_auc_score

logger = logging.getLogger(__name__)

__script_meta__: dict[str, object] = {
    "task": "classify-add-vs-other",
    "kind": "analysis",
    "provides": [
        "OUTCOME_COLUMN",
        "ADD_LEVEL",
        "OTHER_LEVELS",
        "DATA_VERSION",
        "build_binarize",
        "run_elastic_net",
        "run_linear_svm",
        "fold_auc_table",
        "paired_fold_comparison",
        "pooled_mcc",
        "de_cross_check",
        "main",
    ],
    "uses": [
        "project_loader",
        "classification",
        "classification_svm",
        "common.data_loading",
        "common.result_io",
    ],
    "seeded_from": {
        "elastic_net": {"template": "classification", "version": "0.3"},
        "linear_svm": {"template": "classification-svm", "version": "0.3"},
        "cache": {"template": "result-io", "version": "0.2"},
    },
    "description": (
        "PRIMARY plasma-EV classification (FIRST PASS, no null): ADD (10, positive) vs "
        "pooled Other {PDD,PDCN,HCN} (30) on the delivered log2 40x2340 analysis "
        "subset. Runs the elastic-net logistic classifier and the linear SVM on the "
        "same seed / same stratified individual-level folds so their nested-CV AUCs "
        "compare paired per fold; caches each result via result_io.save_cached_result "
        "(analysis classification / classification_svm) and reports ROC-AUC (+ "
        "per-fold AUCs), balanced accuracy, PR-AUC, pooled MCC, the coefficient/weight "
        "stability reads, and DE cross-checks. run_null=False -> exploratory, "
        "coefficients not yet licensed. Cohort confound (ADD == Kerchner among these "
        "40) is a hard bound the null cannot resolve."
    ),
}

OUTCOME_COLUMN = "Condition"
ADD_LEVEL = "ADD"
OTHER_LEVELS: tuple[str, ...] = ("PDD", "PDCN", "HCN")

# Data-version stamp for the cache fingerprint (state/DATA_DESCRIPTION.md).
DATA_VERSION = "protein-data.csv sha256:90a0462eb03a02bd"

DEFAULT_SEED = 20260915
DEFAULT_CACHE_ROOT = Path("results")
DEFAULT_SUMMARY_DIR = Path("results/classification_add_vs_other")
DEFAULT_DE_TABLE = Path(
    "results/differential_abundance/add-vs-other-moderated/contrast_table.csv"
)

# Seed-noise band on a ~40-sample nested AUC (conventions/statistics.md): a paired
# per-fold mean difference inside this band is a TIE, decided by the reading needed,
# not the higher mean.
SEED_NOISE_AUC = 0.03

# Common CV / analysis knobs (shared by both classifiers so their folds coincide).
N_SPLITS = 5
N_REPEATS = 5
STABILITY_REPEATS = 10
SVM_TOP_K = 20


def build_binarize() -> LevelMap:
    """The ADD-vs-Other binarization (elastic-net module): ADD positive, Other neg."""
    return LevelMap(positive=(ADD_LEVEL,), negative=OTHER_LEVELS)


def build_binarize_svm() -> SvmLevelMap:
    """The ADD-vs-Other binarization (SVM module's own LevelMap type)."""
    return SvmLevelMap(positive=(ADD_LEVEL,), negative=OTHER_LEVELS)


def _assert_expected_split(dataset: Dataset) -> None:
    """Fail loud unless Condition is exactly 10 ADD + 30 Other (a data change guard)."""
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
    logger.info("outcome: %d ADD (positive) vs %d Other (negative)", n_add, n_other)


def run_elastic_net(dataset: Dataset, *, seed: int) -> ClassificationResult:
    """Elastic-net logistic classifier, FIRST PASS (run_null=False)."""
    result = classify(
        dataset,
        outcome=OUTCOME_COLUMN,
        binarize=build_binarize(),
        generalization_target="individuals",
        n_splits=N_SPLITS,
        n_repeats=N_REPEATS,
        stability_repeats=STABILITY_REPEATS,
        run_null=False,
        random_state=seed,
    )
    logger.info(
        "elastic net: nested AUC=%.3f (sd %.3f), balacc=%.3f, PR-AUC=%.3f, "
        "best_C=%.4g best_l1=%.2f, selected=%d features",
        result.cv_auc,
        result.cv_auc_sd,
        result.cv_balanced_accuracy,
        result.cv_average_precision,
        result.best_c,
        result.best_l1_ratio,
        len(result.coefficients),
    )
    return result


def run_linear_svm(dataset: Dataset, *, seed: int) -> SVMClassificationResult:
    """Linear-SVM classifier, FIRST PASS (run_null=False), same seed as the EN run."""
    result = classify_svm(
        dataset,
        outcome=OUTCOME_COLUMN,
        binarize=build_binarize_svm(),
        generalization_target="individuals",
        top_k=SVM_TOP_K,
        n_splits=N_SPLITS,
        n_repeats=N_REPEATS,
        stability_repeats=STABILITY_REPEATS,
        run_null=False,
        random_state=seed,
    )
    logger.info(
        "linear SVM: nested AUC=%.3f (sd %.3f), balacc=%.3f, PR-AUC=%.3f, best_C=%.4g, "
        "support vectors=%d (neg %d / pos %d), plateau_start=%.4g, folds<plateau=%s",
        result.cv_auc,
        result.cv_auc_sd,
        result.cv_balanced_accuracy,
        result.cv_average_precision,
        result.best_c,
        result.n_support_vectors,
        result.n_support_negative,
        result.n_support_positive,
        result.plateau_start_c if result.plateau_start_c is not None else float("nan"),
        result.n_folds_sub_plateau,
    )
    return result


@dataclass(frozen=True)
class _FoldAuc:
    repeat: int
    fold: int
    test_indices: tuple[int, ...]
    auc: float


def fold_auc_table(
    y_true_list: list[np.ndarray],
    y_score_list: list[np.ndarray],
    repeats: list[int],
    folds: list[int],
    test_indices_list: list[np.ndarray],
) -> list[_FoldAuc]:
    """Per-outer-fold AUC, keyed by (repeat, fold) and its held-out test indices."""
    out: list[_FoldAuc] = []
    for y_true, y_score, repeat, fold, test_idx in zip(
        y_true_list, y_score_list, repeats, folds, test_indices_list, strict=True
    ):
        out.append(
            _FoldAuc(
                repeat=int(repeat),
                fold=int(fold),
                test_indices=tuple(int(i) for i in test_idx),
                auc=float(roc_auc_score(y_true, y_score)),
            )
        )
    return out


def _en_fold_aucs(result: ClassificationResult) -> list[_FoldAuc]:
    return fold_auc_table(
        [fp.y_true for fp in result.fold_predictions],
        [fp.y_prob for fp in result.fold_predictions],
        [fp.repeat for fp in result.fold_predictions],
        [fp.fold for fp in result.fold_predictions],
        [fp.test_indices for fp in result.fold_predictions],
    )


def _svm_fold_aucs(result: SVMClassificationResult) -> list[_FoldAuc]:
    return fold_auc_table(
        [fp.y_true for fp in result.fold_predictions],
        [fp.y_score for fp in result.fold_predictions],
        [fp.repeat for fp in result.fold_predictions],
        [fp.fold for fp in result.fold_predictions],
        [fp.test_indices for fp in result.fold_predictions],
    )


def paired_fold_comparison(
    en_folds: list[_FoldAuc], svm_folds: list[_FoldAuc]
) -> dict[str, Any]:
    """Pair EN vs SVM outer folds by (repeat, fold), verify identical test indices,
    and summarize the per-fold AUC difference against the seed-noise band."""
    en_by_key = {(f.repeat, f.fold): f for f in en_folds}
    svm_by_key = {(f.repeat, f.fold): f for f in svm_folds}
    if set(en_by_key) != set(svm_by_key):
        raise ValueError(
            "Elastic-net and SVM outer folds do not share the same (repeat, fold) "
            "keys; the paired comparison requires the same seed and CV design."
        )
    diffs: list[float] = []
    per_fold: list[dict[str, Any]] = []
    for key in sorted(en_by_key):
        en_f = en_by_key[key]
        svm_f = svm_by_key[key]
        if en_f.test_indices != svm_f.test_indices:
            raise ValueError(
                f"Fold {key} has different held-out test indices between the two "
                f"models; the outer CV is not shared (cannot pair)."
            )
        diff = en_f.auc - svm_f.auc
        diffs.append(diff)
        per_fold.append(
            {
                "repeat": key[0],
                "fold": key[1],
                "auc_elastic_net": round(en_f.auc, 4),
                "auc_linear_svm": round(svm_f.auc, 4),
                "diff_en_minus_svm": round(diff, 4),
            }
        )
    diff_arr = np.asarray(diffs, dtype=float)
    mean_diff = float(diff_arr.mean())
    n_en_win = int((diff_arr > SEED_NOISE_AUC).sum())
    n_svm_win = int((diff_arr < -SEED_NOISE_AUC).sum())
    n_tie = int((np.abs(diff_arr) <= SEED_NOISE_AUC).sum())

    # Wilcoxon signed-rank on the paired per-fold differences (a p if one is wanted).
    from scipy.stats import wilcoxon

    if np.allclose(diff_arr, 0.0):
        wilcoxon_p: float | None = 1.0
    else:
        try:
            wilcoxon_p = float(wilcoxon(diff_arr).pvalue)
        except ValueError:
            wilcoxon_p = None

    if abs(mean_diff) <= SEED_NOISE_AUC:
        verdict = (
            f"TIE — the mean paired AUC difference ({mean_diff:+.4f}) is inside the "
            f"seed-noise band (+/-{SEED_NOISE_AUC:.2f} at n~40). Neither model wins on "
            f"performance; choose by the reading needed (sparse minimal set -> elastic "
            f"net; dense module weights -> SVM)."
        )
    elif mean_diff > 0:
        verdict = (
            f"ELASTIC NET higher by {mean_diff:+.4f} (beyond the "
            f"+/-{SEED_NOISE_AUC:.2f} seed-noise band) on the shared folds — check "
            f"sign-consistency across a few seeds before calling it a real preference."
        )
    else:
        verdict = (
            f"LINEAR SVM higher by {-mean_diff:+.4f} (beyond the "
            f"+/-{SEED_NOISE_AUC:.2f} seed-noise band) on the shared folds — check "
            f"sign-consistency across a few seeds before calling it a real preference."
        )
    return {
        "n_folds": int(diff_arr.size),
        "mean_paired_diff_en_minus_svm": round(mean_diff, 4),
        "median_paired_diff_en_minus_svm": round(float(np.median(diff_arr)), 4),
        "sd_paired_diff": round(float(diff_arr.std()), 4),
        "seed_noise_band": SEED_NOISE_AUC,
        "n_folds_en_wins": n_en_win,
        "n_folds_svm_wins": n_svm_win,
        "n_folds_tie_within_band": n_tie,
        "wilcoxon_signed_rank_p": (
            None if wilcoxon_p is None else round(wilcoxon_p, 4)
        ),
        "verdict": verdict,
        "per_fold": per_fold,
    }


def pooled_mcc(
    y_true_list: list[np.ndarray],
    y_pred_list: list[np.ndarray],
) -> float:
    """Matthews correlation over the pooled held-out predictions of every outer fold.

    Each sample appears once per outer repeat, so the pool has n_samples * n_repeats
    entries; MCC on the pool is a single imbalance-robust aggregate of the held-out
    hard-label agreement.
    """
    y_true = np.concatenate(y_true_list)
    y_pred = np.concatenate(y_pred_list)
    return float(matthews_corrcoef(y_true, y_pred))


def _en_pooled_mcc(result: ClassificationResult) -> float:
    return pooled_mcc(
        [fp.y_true for fp in result.fold_predictions],
        [(fp.y_prob >= 0.5).astype(int) for fp in result.fold_predictions],
    )


def _svm_pooled_mcc(result: SVMClassificationResult) -> float:
    return pooled_mcc(
        [fp.y_true for fp in result.fold_predictions],
        [(fp.y_score >= 0.0).astype(int) for fp in result.fold_predictions],
    )


def _parse_gene(feature_id: str) -> str:
    """Return the protein short-name (``NAME`` of ``db|ACC|NAME_ORGANISM``), else id."""
    parts = feature_id.split("|")
    if len(parts) >= 3:
        name = parts[2]
        return name.rsplit("_", 1)[0] if "_" in name else name
    return feature_id


def de_cross_check(
    coefficients: pd.DataFrame,
    de_table: pd.DataFrame,
    *,
    coef_column: str,
    top_n: int,
    q_threshold: float = 0.05,
) -> dict[str, Any]:
    """Overlap of the top classifier features with the DE hits + direction agreement.

    ``coef_column`` is the signed all-data coefficient/weight column. The classifier's
    sign is *toward ADD* (the positive class); the DE effect is log2FC(ADD - Other), so
    a shared feature's signs agree iff both are the same sign.
    """
    de_by_feature = de_table.set_index("feature")
    de_hits = set(de_by_feature.index[de_by_feature["q"] < q_threshold].tolist())

    top = coefficients.head(top_n)
    rows: list[dict[str, Any]] = []
    n_in_de_hits = 0
    n_sign_checked = 0
    n_sign_agree = 0
    for _, r in top.iterrows():
        feature = str(r["feature"])
        coef = float(r[coef_column])
        in_de = feature in de_by_feature.index
        is_hit = feature in de_hits
        de_effect = (
            float(de_by_feature.loc[feature, "effect"]) if in_de else float("nan")
        )
        de_q = float(de_by_feature.loc[feature, "q"]) if in_de else float("nan")
        agree: bool | None = None
        if in_de and np.isfinite(de_effect) and coef != 0.0:
            agree = bool(np.sign(coef) == np.sign(de_effect))
            n_sign_checked += 1
            n_sign_agree += int(agree)
        if is_hit:
            n_in_de_hits += 1
        rows.append(
            {
                "feature": feature,
                "gene": _parse_gene(feature),
                "classifier_coef": round(coef, 4),
                "coef_sign_toward": "ADD" if coef > 0 else "Other",
                "in_de_table": in_de,
                "de_bh_q_lt_0.05": is_hit,
                "de_log2fc_add_minus_other": (
                    None if not np.isfinite(de_effect) else round(de_effect, 4)
                ),
                "de_bh_q": None if not np.isfinite(de_q) else de_q,
                "sign_agrees_with_de": agree,
            }
        )
    return {
        "top_n": top_n,
        "n_top_in_de_hits_q_lt_0.05": n_in_de_hits,
        "n_de_hits_q_lt_0.05_total": len(de_hits),
        "sign_agreement_shared_features": (
            None
            if n_sign_checked == 0
            else {
                "n_checked": n_sign_checked,
                "n_agree": n_sign_agree,
                "fraction_agree": round(n_sign_agree / n_sign_checked, 3),
            }
        ),
        "top_features": rows,
    }


def _en_params(seed: int) -> dict[str, object]:
    return {
        "outcome": OUTCOME_COLUMN,
        "binarize": "LevelMap",
        "positive_label": ADD_LEVEL,
        "negative_levels": sorted(OTHER_LEVELS),
        "method": "elastic_net_logistic",
        "generalization_target": "individuals",
        "cv_scheme": "stratified_kfold",
        "n_splits": N_SPLITS,
        "n_repeats": N_REPEATS,
        "grouped": False,
        "stability_repeats": STABILITY_REPEATS,
        "tuning_metric": "roc_auc",
        "select": "best",
        "c_grid": list(EN_C_GRID),
        "l1_ratios": list(EN_L1_RATIOS),
        "feature_list": None,
        "run_null": False,
        "null_permutation": "units",
        "seed": seed,
    }


def _svm_params(seed: int) -> dict[str, object]:
    return {
        "outcome": OUTCOME_COLUMN,
        "binarize": "LevelMap",
        "positive_label": ADD_LEVEL,
        "negative_levels": sorted(OTHER_LEVELS),
        "method": "linear_svm",
        "generalization_target": "individuals",
        "cv_scheme": "stratified_kfold",
        "n_splits": N_SPLITS,
        "n_repeats": N_REPEATS,
        "grouped": False,
        "stability_repeats": STABILITY_REPEATS,
        "tuning_metric": "roc_auc",
        "select": "best",
        "c_grid": list(SVM_C_GRID),
        "top_k": SVM_TOP_K,
        "feature_list": None,
        "run_null": False,
        "null_permutation": "units",
        "seed": seed,
    }


def _en_stability_top(result: ClassificationResult, n: int) -> list[dict[str, Any]]:
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
            "coef_sign_toward": "ADD" if float(r["coef"]) > 0 else "Other",
            "selection_frequency": round(float(r["selection_frequency"]), 3),
            "sign_consistency": round(float(r["sign_consistency"]), 3),
            "coef_iqr": [
                round(float(r["coef_q25"]), 4),
                round(float(r["coef_q75"]), 4),
            ],
        }
        for _, r in result.coefficients[cols].head(n).iterrows()
    ]


def _svm_stability_top(result: SVMClassificationResult, n: int) -> list[dict[str, Any]]:
    cols = [
        "feature",
        "coef",
        "top_k_frequency",
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
            "coef_sign_toward": "ADD" if float(r["coef"]) > 0 else "Other",
            "top_k_frequency": round(float(r["top_k_frequency"]), 3),
            "sign_consistency": round(float(r["sign_consistency"]), 3),
            "coef_iqr": [
                round(float(r["coef_q25"]), 4),
                round(float(r["coef_q75"]), 4),
            ],
        }
        for _, r in result.coefficients[cols].head(n).iterrows()
    ]


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("data/protein-data.csv"))
    parser.add_argument("--metadata", type=Path, default=Path("data/metadata.csv"))
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--cache-root", type=Path, default=DEFAULT_CACHE_ROOT)
    parser.add_argument("--summary-dir", type=Path, default=DEFAULT_SUMMARY_DIR)
    parser.add_argument("--de-table", type=Path, default=DEFAULT_DE_TABLE)
    parser.add_argument("--top-n", type=int, default=25)
    args = parser.parse_args(argv)

    dataset = load_analysis_subset(load_paired_dataset(args.data, args.metadata))
    logger.info(
        "loaded analysis subset: %d samples x %d features, scale=%s",
        dataset.abundances.shape[0],
        dataset.abundances.shape[1],
        dataset.scale,
    )
    _assert_expected_split(dataset)

    en_result = run_elastic_net(dataset, seed=args.seed)
    svm_result = run_linear_svm(dataset, seed=args.seed)

    en_meta: ResultMeta = save_cached_result(
        en_result,
        cache_root=args.cache_root,
        analysis="classification",
        data_version=DATA_VERSION,
        params=_en_params(args.seed),
        seed=args.seed,
        label="add-vs-other elastic-net first-pass (no null)",
    )
    svm_meta: ResultMeta = save_cached_result(
        svm_result,
        cache_root=args.cache_root,
        analysis="classification_svm",
        data_version=DATA_VERSION,
        params=_svm_params(args.seed),
        seed=args.seed,
        label="add-vs-other linear-SVM first-pass (no null)",
    )
    logger.info(
        "cached: elastic net -> %s (%s); linear SVM -> %s (%s)",
        en_meta.fingerprint,
        Path(args.cache_root) / "classification" / en_meta.fingerprint,
        svm_meta.fingerprint,
        Path(args.cache_root) / "classification_svm" / svm_meta.fingerprint,
    )

    en_folds = _en_fold_aucs(en_result)
    svm_folds = _svm_fold_aucs(svm_result)
    comparison = paired_fold_comparison(en_folds, svm_folds)

    de_table = pd.read_csv(args.de_table)
    en_cross = de_cross_check(
        en_result.coefficients, de_table, coef_column="coef", top_n=args.top_n
    )
    svm_cross = de_cross_check(
        svm_result.coefficients, de_table, coef_column="coef", top_n=args.top_n
    )

    summary: dict[str, Any] = {
        "analysis": "PRIMARY classification ADD vs Other (FIRST PASS, no null)",
        "data_version": DATA_VERSION,
        "seed": args.seed,
        "outcome": "ADD (10, positive) vs Other {PDD,PDCN,HCN} (30, negative)",
        "cv": (
            f"stratified {N_SPLITS}-fold x {N_REPEATS} repeats at the individual (row) "
            "level; no groups (subjects are singletons)"
        ),
        "run_null": False,
        "phase_cap": "exploratory (performance not yet tested against a label-shuffle "
        "null; coefficients/weights not yet licensed)",
        "elastic_net": {
            "cache_fingerprint": en_meta.fingerprint,
            "cache_path": str(
                Path(args.cache_root) / "classification" / en_meta.fingerprint
            ),
            "cv_roc_auc": round(en_result.cv_auc, 4),
            "cv_roc_auc_sd_over_folds": round(en_result.cv_auc_sd, 4),
            "cv_balanced_accuracy": round(en_result.cv_balanced_accuracy, 4),
            "cv_pr_auc": round(en_result.cv_average_precision, 4),
            "pooled_oof_mcc": round(_en_pooled_mcc(en_result), 4),
            "best_c": en_result.best_c,
            "best_l1_ratio": en_result.best_l1_ratio,
            "n_features_selected_nonzero": len(en_result.coefficients),
            "per_fold_auc": [round(f.auc, 4) for f in en_folds],
            "stability_top": _en_stability_top(en_result, args.top_n),
        },
        "linear_svm": {
            "cache_fingerprint": svm_meta.fingerprint,
            "cache_path": str(
                Path(args.cache_root) / "classification_svm" / svm_meta.fingerprint
            ),
            "cv_roc_auc": round(svm_result.cv_auc, 4),
            "cv_roc_auc_sd_over_folds": round(svm_result.cv_auc_sd, 4),
            "cv_balanced_accuracy": round(svm_result.cv_balanced_accuracy, 4),
            "cv_pr_auc": round(svm_result.cv_average_precision, 4),
            "pooled_oof_mcc": round(_svm_pooled_mcc(svm_result), 4),
            "best_c": svm_result.best_c,
            "plateau_start_c": svm_result.plateau_start_c,
            "n_folds_sub_plateau": svm_result.n_folds_sub_plateau,
            "n_support_vectors": svm_result.n_support_vectors,
            "n_support_negative": svm_result.n_support_negative,
            "n_support_positive": svm_result.n_support_positive,
            "per_fold_auc": [round(f.auc, 4) for f in svm_folds],
            "stability_top": _svm_stability_top(svm_result, args.top_n),
        },
        "paired_per_fold_comparison": comparison,
        "de_cross_check_elastic_net": en_cross,
        "de_cross_check_linear_svm": svm_cross,
        "caveats": [
            "EXPLORATORY, capped: performance NOT yet tested against a label-shuffle "
            "null; coefficients/weights flagged 'not tested against a null' until the "
            "null pass (run_null=True), which is also the leakage detector.",
            "COHORT CONFOUND (hard): ADD is 100% Kerchner cohort (finding 0001). Among "
            "these 40 samples ADD and Kerchner are the SAME partition, so above-null "
            "performance CANNOT be attributed to disease vs collection-site. The "
            "label-shuffle null does NOT resolve this (it permutes the shared "
            "ADD/Kerchner partition).",
            "Small n (10 positives): wide CV CIs; treat point estimates cautiously.",
            "Matrix is already limma removeBatchEffect-corrected on Cohort (finding "
            "0004); the pre-correction matrix is unavailable.",
        ],
    }

    args.summary_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.summary_dir / "first_pass_comparison.json"
    with summary_path.open("w") as fh:
        json.dump(summary, fh, indent=2)
    logger.info("wrote summary -> %s", summary_path)
    logger.info(
        "PAIRED: mean EN-SVM AUC diff=%+.4f over %d folds -> %s",
        comparison["mean_paired_diff_en_minus_svm"],
        comparison["n_folds"],
        comparison["verdict"],
    )


if __name__ == "__main__":
    main()
