"""QC: verify sample SEX (primary) and AGE (secondary) metadata against the proteome.

Independent integrity check that the ``Sex``/``isMale`` and ``Age`` annotations in
``data/metadata.csv`` are consistent with the actual protein abundances, and flags any
sample that looks misannotated or swapped.

Design (see the task brief and ``state/DATA_DESCRIPTION.md``):

SEX (primary)
    1. Search the 2340 analysis features for Y-chromosome / sex-specific proteins by
       UniProt entry-name token (DDX3Y, RS4Y1/2, UTY, KDM5D, USP9Y, NLGNY, IF1AY, ZFY,
       TYB4Y/TMSB4Y, PRKY, ...). Report which are actually PRESENT.
    2. If >=1 Y-linked protein is present -> a transparent threshold/abundance call.
       If NONE are present (the case for this plasma dataset) -> fall back to a
       leakage-safe multivariate sex classifier: StandardScaler + regularized logistic
       regression, evaluated by leave-one-out CV (all learning inside each fold), with
       a label-shuffle permutation null (also the leakage detector) and a bootstrap CI
       on AUC. The primary model is ridge-end elastic net (``l1_ratio=0``); an
       ``l1_ratio=0.5`` elastic net is reported beside it.
    3. Per-sample out-of-fold P(male); high-confidence mismatches only are flagged.

AGE (secondary, low-powered with n=40 -- deliberately conservative)
    Leakage-safe RidgeCV LOO predicted age vs annotated age: Pearson/Spearman with CI,
    a permutation null, and GROSS standardized-residual outliers only. Catches gross
    swaps, not +/-5-year errors.

No leakage: every scaler / alpha-selection / coefficient is fit inside the training
fold only, never on the full matrix. CV target = per-sample sex/age (samples are
independent; no repeated subjects), so plain LOO with a label-shuffle null matches the
generalization target.

Does NOT create findings and does NOT modify metadata -- analyze and report only.
Writes results to ``results/qc_sex_age_check/``.

Run::

    ./.venv/bin/python scripts/scratch/qc_sex_age_annotation_check.py
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import warnings
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
from numpy.random import Generator, default_rng
from scipy import stats
from sklearn.base import BaseEstimator
from sklearn.linear_model import LogisticRegression, RidgeCV
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from project_loader import load_analysis_subset, load_paired_dataset

logger = logging.getLogger(__name__)

__script_meta__: dict[str, object] = {
    "task": "QC: verify SEX/AGE metadata annotations against the proteome",
    "kind": "analysis",
    "provides": [
        "results/qc_sex_age_check/sex_per_sample.csv",
        "results/qc_sex_age_check/age_per_sample.csv",
        "results/qc_sex_age_check/summary.json",
    ],
    "uses": ["project_loader", "common.data_loading"],
    "seeded_from": {
        "template": "leakage-safe-linear-classifier (local adaptation; no lib/)",
        "version": None,
    },
    "description": (
        "Independent QC of Sex/isMale (primary) and Age (secondary) metadata vs the "
        "40x2340 log2 plasma proteome. No Y-linked proteins are present, so sex uses a "
        "leakage-safe LOO logistic classifier with a label-shuffle null; age uses LOO "
        "RidgeCV predicted-age correlation with gross-residual outlier flagging. "
        "Reports per-sample concordance + explicit mismatch lists; creates no findings "
        "and never edits metadata."
    ),
}

# UniProt ENTRY-NAME tokens (the NAME in db|ACC|NAME_HUMAN) for the requested
# Y-chromosome / sex-specific proteins. Matched case-insensitively as substrings.
Y_LINKED_ENTRY_TOKENS: tuple[str, ...] = (
    "DDX3Y",
    "RS4Y1",  # RPS4Y1
    "RS4Y2",  # RPS4Y2
    "UTY",
    "KDM5D",
    "USP9Y",
    "NLGNY",  # NLGN4Y
    "IF1AY",  # EIF1AY
    "ZFY",
    "TYB4Y",  # TMSB4Y
    "TMSB4Y",
    "PRKY",
    "TXLNGY",
    "PCDH11Y",
    "AMELY",
)

SexModelName = Literal["ridge_logistic", "elasticnet_logistic"]


@dataclass(frozen=True)
class Measurement:
    """A single reported statistic with its uncertainty and test provenance."""

    metric: str
    value: float
    ci_low: float | None
    ci_high: float | None
    p_value: float | None
    correction: str | None
    test: str
    n: int


def _sex_pipeline(name: SexModelName) -> Pipeline:
    """Build a leakage-safe scaler+logistic pipeline (all learning is fit in-fold)."""
    if name == "ridge_logistic":
        clf: BaseEstimator = LogisticRegression(C=0.01, max_iter=5000)  # l2 (ridge) end
    elif name == "elasticnet_logistic":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            clf = LogisticRegression(
                penalty="elasticnet",
                solver="saga",
                l1_ratio=0.5,
                C=0.05,
                max_iter=5000,
                tol=1e-3,
            )
    else:  # pragma: no cover - exhaustive
        raise ValueError(f"unknown sex model {name!r}")
    return Pipeline([("scaler", StandardScaler()), ("clf", clf)])


def _loo_oof_proba(pipe: Pipeline, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Leave-one-out out-of-fold P(class=1); scaler+model refit per fold (no leakage)."""  # noqa: E501
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        proba = cross_val_predict(
            pipe, x, y, cv=LeaveOneOut(), method="predict_proba"
        )
    return np.asarray(proba)[:, 1]


def _auc(y: np.ndarray, score: np.ndarray) -> float:
    from sklearn.metrics import roc_auc_score

    return float(roc_auc_score(y, score))


def _balanced_accuracy(y: np.ndarray, pred: np.ndarray) -> float:
    from sklearn.metrics import balanced_accuracy_score

    return float(balanced_accuracy_score(y, pred))


def _bootstrap_auc_ci(
    y: np.ndarray, score: np.ndarray, rng: Generator, n_boot: int
) -> tuple[float, float]:
    """Percentile bootstrap CI for LOO AUC (resample samples with replacement)."""
    n = len(y)
    aucs: list[float] = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yb, sb = y[idx], score[idx]
        if len(np.unique(yb)) < 2:
            continue
        aucs.append(_auc(yb, sb))
    if not aucs:
        return (float("nan"), float("nan"))
    lo, hi = np.percentile(aucs, [2.5, 97.5])
    return (float(lo), float(hi))


def _permutation_null_auc(
    pipe: Pipeline, x: np.ndarray, y: np.ndarray, rng: Generator, n_perm: int
) -> tuple[float, np.ndarray]:
    """Label-shuffle null: recompute LOO AUC on permuted labels. Returns (p, null)."""
    obs = _auc(y, _loo_oof_proba(pipe, x, y))
    null = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        yp = rng.permutation(y)
        null[i] = _auc(yp, _loo_oof_proba(pipe, x, yp))
    p = float((1 + int(np.sum(null >= obs))) / (n_perm + 1))
    return p, null


@dataclass(frozen=True)
class SexResult:
    model: str
    auc: Measurement
    balanced_accuracy: Measurement
    per_sample: pd.DataFrame
    mismatches: pd.DataFrame
    null_auc_mean: float


def run_sex_check(
    x: np.ndarray,
    meta: pd.DataFrame,
    model: SexModelName,
    *,
    rng: Generator,
    n_perm: int,
    n_boot: int,
    hc_prob: float,
) -> SexResult:
    """Leakage-safe LOO sex classifier + label-shuffle null + per-sample mismatch call.

    ``hc_prob`` is the high-confidence threshold: a sample is a HIGH-CONFIDENCE mismatch
    only if its out-of-fold probability for the OPPOSITE sex exceeds ``hc_prob``.
    """
    y = meta["isMale"].astype(int).to_numpy()
    n = len(y)
    pipe = _sex_pipeline(model)
    proba_male = _loo_oof_proba(pipe, x, y)
    pred = (proba_male > 0.5).astype(int)

    auc = _auc(y, proba_male)
    lo, hi = _bootstrap_auc_ci(y, proba_male, rng, n_boot)
    if n_perm > 0:
        p_auc, null = _permutation_null_auc(pipe, x, y, rng, n_perm)
        null_mean = float(np.mean(null))
    else:  # skip the (expensive) null for a secondary comparison model
        p_auc, null_mean = None, float("nan")
    bacc = _balanced_accuracy(y, pred)

    per = pd.DataFrame(
        {
            "SampleID": meta.index.to_numpy(),
            "Condition": meta["Condition"].to_numpy(),
            "annotated_Sex": meta["Sex"].to_numpy(),
            "annotated_isMale": y,
            "oof_P_male": np.round(proba_male, 4),
            "predicted_isMale": pred,
            "predicted_Sex": np.where(pred == 1, "Male", "Female"),
            "concordant": pred == y,
            # confidence toward the OPPOSITE of the annotation:
            "opp_prob": np.where(y == 1, proba_male < 0.5, proba_male > 0.5)
            * np.abs(proba_male - (1 - y)),
        }
    )
    # opp confidence = P assigned to the opposite class
    per["P_opposite_class"] = np.round(
        np.where(y == 1, 1 - proba_male, proba_male), 4
    )
    per = per.drop(columns=["opp_prob"])

    mism = per[(~per["concordant"]) & (per["P_opposite_class"] > hc_prob)].copy()

    return SexResult(
        model=model,
        auc=Measurement(
            metric="LOO_ROC_AUC",
            value=round(auc, 4),
            ci_low=round(lo, 4),
            ci_high=round(hi, 4),
            p_value=p_auc,
            correction=(
                "label-shuffle permutation null (one-sided)"
                if n_perm > 0
                else "no null (secondary comparison model)"
            ),
            test=(
                f"LOO CV {model}; {n_perm}-shuffle null; {n_boot}-boot CI"
                if n_perm > 0
                else f"LOO CV {model}; {n_boot}-boot CI (null skipped)"
            ),
            n=n,
        ),
        balanced_accuracy=Measurement(
            metric="LOO_balanced_accuracy",
            value=round(bacc, 4),
            ci_low=None,
            ci_high=None,
            p_value=None,
            correction=None,
            test=f"LOO CV {model}",
            n=n,
        ),
        per_sample=per,
        mismatches=mism,
        null_auc_mean=null_mean,
    )


@dataclass(frozen=True)
class AgeResult:
    per_sample: pd.DataFrame
    pearson: Measurement
    spearman: Measurement
    rmse_years: float
    r2: float
    outliers: pd.DataFrame
    null_r_mean: float


def run_age_check(
    x: np.ndarray,
    meta: pd.DataFrame,
    *,
    rng: Generator,
    n_perm: int,
    resid_z: float,
) -> AgeResult:
    """Leakage-safe LOO RidgeCV predicted age vs annotated age + gross-outlier flag."""
    age = meta["Age"].astype(float).to_numpy()
    n = len(age)
    alphas = np.logspace(1, 6, 20)
    pipe = Pipeline(
        [("scaler", StandardScaler()), ("ridge", RidgeCV(alphas=alphas))]
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pred = cross_val_predict(pipe, x, age, cv=LeaveOneOut())
    pred = np.asarray(pred, dtype=float)

    r, _ = stats.pearsonr(age, pred)
    rho, _ = stats.spearmanr(age, pred)
    # Fisher-z CI for Pearson r
    if abs(r) < 1.0:
        z = np.arctanh(r)
        se = 1.0 / np.sqrt(n - 3)
        lo, hi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
    else:  # pragma: no cover
        lo, hi = r, r

    # permutation null on Pearson r (shuffle annotated age, refit LOO)
    null = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        ap = rng.permutation(age)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pp = cross_val_predict(pipe, x, ap, cv=LeaveOneOut())
        null[i], _ = stats.pearsonr(ap, np.asarray(pp, dtype=float))
    p_r = float((1 + int(np.sum(null >= r))) / (n_perm + 1))

    resid = age - pred
    z_resid = (resid - resid.mean()) / resid.std(ddof=1)
    ss_res = float(np.sum(resid**2))
    ss_tot = float(np.sum((age - age.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot
    rmse = float(np.sqrt(np.mean(resid**2)))

    per = pd.DataFrame(
        {
            "SampleID": meta.index.to_numpy(),
            "Condition": meta["Condition"].to_numpy(),
            "annotated_Age": np.round(age, 1),
            "predicted_Age": np.round(pred, 1),
            "residual_years": np.round(resid, 1),
            "residual_z": np.round(z_resid, 2),
        }
    )
    outliers = per[np.abs(per["residual_z"]) > resid_z].copy()

    return AgeResult(
        per_sample=per,
        pearson=Measurement(
            metric="pearson_r(annotated_age, LOO_predicted_age)",
            value=round(float(r), 4),
            ci_low=round(float(lo), 4),
            ci_high=round(float(hi), 4),
            p_value=p_r,
            correction="label-shuffle permutation null (one-sided)",
            test=f"LOO RidgeCV predicted age; {n_perm}-shuffle null; Fisher-z CI",
            n=n,
        ),
        spearman=Measurement(
            metric="spearman_rho(annotated_age, LOO_predicted_age)",
            value=round(float(rho), 4),
            ci_low=None,
            ci_high=None,
            p_value=None,
            correction=None,
            test="LOO RidgeCV predicted age",
            n=n,
        ),
        rmse_years=round(rmse, 2),
        r2=round(r2, 4),
        outliers=outliers,
        null_r_mean=float(np.mean(null)),
    )


def search_sex_proteins(feature_names: np.ndarray) -> pd.DataFrame:
    """Report which requested Y-linked / sex-specific proteins are present."""
    names = [str(n) for n in feature_names]
    upper = [n.upper() for n in names]
    rows: list[dict[str, object]] = []
    for tok in Y_LINKED_ENTRY_TOKENS:
        hits = [names[i] for i, u in enumerate(upper) if tok in u]
        rows.append({"entry_token": tok, "present": bool(hits), "matches": hits})
    return pd.DataFrame(rows)


def _measurement_json(m: Measurement) -> dict[str, object]:
    return asdict(m)


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("data/protein-data.csv"))
    parser.add_argument("--metadata", type=Path, default=Path("data/metadata.csv"))
    parser.add_argument(
        "--out-dir", type=Path, default=Path("results/qc_sex_age_check")
    )
    parser.add_argument("--seed", type=int, default=20260915)
    parser.add_argument("--n-perm", type=int, default=1000)
    parser.add_argument("--n-boot", type=int, default=2000)
    parser.add_argument("--n-perm-age", type=int, default=1000)
    parser.add_argument(
        "--hc-prob",
        type=float,
        default=0.90,
        help="high-confidence mismatch: P(opposite sex) must exceed this",
    )
    parser.add_argument("--resid-z", type=float, default=3.0)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    ds = load_analysis_subset(load_paired_dataset(args.data, args.metadata))
    x = ds.abundances
    meta = ds.metadata
    logger.info("loaded %d samples x %d features", *x.shape)

    # (a) sex-informative protein search
    protein_search = search_sex_proteins(ds.feature_names)
    n_present = int(protein_search["present"].sum())
    logger.info("Y-linked / sex-specific proteins present: %d", n_present)
    protein_search.to_csv(args.out_dir / "sex_protein_search.csv", index=False)

    # (b) sex classifier(s) -- fallback path since no Y-linked proteins present
    sex_primary = run_sex_check(
        x,
        meta,
        "ridge_logistic",
        rng=default_rng(args.seed),
        n_perm=args.n_perm,
        n_boot=args.n_boot,
        hc_prob=args.hc_prob,
    )
    # Secondary comparison model: skip the (saga, expensive) null -- the primary ridge
    # model carries the label-shuffle null; this is reported only as a sensitivity read.
    sex_enet = run_sex_check(
        x,
        meta,
        "elasticnet_logistic",
        rng=default_rng(args.seed + 1),
        n_perm=0,
        n_boot=args.n_boot,
        hc_prob=args.hc_prob,
    )
    sex_primary.per_sample.to_csv(args.out_dir / "sex_per_sample.csv", index=False)
    sex_primary.mismatches.to_csv(args.out_dir / "sex_mismatches.csv", index=False)
    logger.info(
        "SEX primary(ridge): AUC=%.3f [%.3f,%.3f] perm-p=%.3g bacc=%.3f null_mean=%.3f",
        sex_primary.auc.value,
        sex_primary.auc.ci_low,
        sex_primary.auc.ci_high,
        sex_primary.auc.p_value,
        sex_primary.balanced_accuracy.value,
        sex_primary.null_auc_mean,
    )
    logger.info(
        "SEX enet(l1=0.5) comparison: AUC=%.3f (null skipped)", sex_enet.auc.value
    )
    logger.info("SEX high-confidence mismatches: %d", len(sex_primary.mismatches))

    # (c) age check
    age = run_age_check(
        x,
        meta,
        rng=default_rng(args.seed + 2),
        n_perm=args.n_perm_age,
        resid_z=args.resid_z,
    )
    age.per_sample.to_csv(args.out_dir / "age_per_sample.csv", index=False)
    age.outliers.to_csv(args.out_dir / "age_outliers.csv", index=False)
    logger.info(
        "AGE: pearson r=%.3f [%.3f,%.3f] perm-p=%.3g rmse=%.1fy r2=%.3f outliers=%d",
        age.pearson.value,
        age.pearson.ci_low,
        age.pearson.ci_high,
        age.pearson.p_value,
        age.rmse_years,
        age.r2,
        len(age.outliers),
    )

    summary: dict[str, object] = {
        "data_version": (
            "protein-data.csv:90a0462eb03a02bd metadata.csv:62f10c7072fe613e"
        ),
        "n_samples": int(x.shape[0]),
        "n_features": int(x.shape[1]),
        "seed": args.seed,
        "params": {
            "n_perm_sex": args.n_perm,
            "n_boot": args.n_boot,
            "n_perm_age": args.n_perm_age,
            "hc_prob": args.hc_prob,
            "resid_z": args.resid_z,
        },
        "sex_proteins_present": n_present,
        "sex_protein_search": protein_search.to_dict(orient="records"),
        "sex_primary_model": sex_primary.model,
        "sex_primary_auc": _measurement_json(sex_primary.auc),
        "sex_primary_balanced_accuracy": _measurement_json(
            sex_primary.balanced_accuracy
        ),
        "sex_primary_null_auc_mean": sex_primary.null_auc_mean,
        "sex_elasticnet_auc": _measurement_json(sex_enet.auc),
        "sex_high_confidence_mismatches": sex_primary.mismatches.to_dict(
            orient="records"
        ),
        "age_pearson": _measurement_json(age.pearson),
        "age_spearman": _measurement_json(age.spearman),
        "age_rmse_years": age.rmse_years,
        "age_r2": age.r2,
        "age_null_r_mean": age.null_r_mean,
        "age_gross_outliers": age.outliers.to_dict(orient="records"),
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    logger.info("wrote results to %s", args.out_dir)


if __name__ == "__main__":
    main()
