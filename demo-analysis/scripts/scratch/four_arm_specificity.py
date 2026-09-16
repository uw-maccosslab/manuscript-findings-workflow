"""FOUR-ARM ADD-SPECIFICITY check for the Boruta-confirmed ADD-vs-Other signature.

The primary contrast pooled {PDD, PDCN, HCN} as "Other", discarding the four-arm
confound-control structure the study was designed around. This script asks the design's
core question: is the ADD-discriminative signature **ADD-SPECIFIC** or
**dementia-GENERIC**? Concretely -- does the down-in-ADD protein set ALSO drop in PDD
(Parkinson's disease dementia, the dementia control), or is the shift specific to ADD?

Two analyses, both EXPLORATORY (n=10/arm; confound-bounded; findings 0001/0003/0004):

  1. **Per-arm profile of the signature.** Split the 40 Boruta-Confirmed proteins into
     the up-in-ADD and down-in-ADD subsets by the sign of the ADD-vs-Other moderated DE
     log2FC. z-score each protein across the 40 samples, then for each SAMPLE take the
     mean z over the up-set and over the down-set -> two per-sample summary scores.
     Report each arm's (ADD/PDD/PDCN/HCN) mean score with a t 95% CI, test ADD vs each
     arm (Welch t) with Hedges' g + 95% CI, BH-corrected across the 6 comparisons.

  2. **ADD-vs-PDD moderated DE (dementia-vs-dementia)** on the full 2340 features (same
     limma-style eBayes as finding 0006). How many of the Boruta-40 (and of each DE-hit
     direction) still separate ADD from PDD? CRITICAL: ADD is 100% Kerchner and PDD is
     ~70% Poston, so ADD-vs-PDD is a near-pure CROSS-COHORT contrast -- even MORE
     cohort-confounded than the pooled contrast. Every ADD-vs-PDD number is flagged as
     inseparable from a site effect.

Deterministic; seed 20260915 (no RNG -- recorded for provenance/fingerprint). Fail-loud.
Caches result tables + meta.json under results/four_arm/<fingerprint>/.

Run::

    ./.venv/bin/python scripts/scratch/four_arm_specificity.py
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import logging
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from common.data_loading import Dataset
from common.result_io import save_cached_result
from differential_abundance import differential_abundance
from project_loader import load_analysis_subset, load_paired_dataset

logger = logging.getLogger(__name__)

__script_meta__: dict[str, object] = {
    "task": "four-arm-add-specificity",
    "kind": "analysis",
    "provides": [
        "FourArmResult",
        "confirmed_features",
        "signature_direction",
        "per_arm_scores",
        "arm_comparisons",
        "add_vs_pdd_de",
        "run_four_arm",
        "main",
    ],
    "uses": ["project_loader", "differential_abundance", "common.data_loading",
             "common.result_io"],
    "seeded_from": {"template": "differential-abundance", "version": "0.1"},
    "description": (
        "Four-arm ADD-specificity check for the Boruta-confirmed ADD-vs-Other "
        "signature: per-arm (ADD/PDD/PDCN/HCN) mean z summary scores of the up/down "
        "signature subsets with Hedges' g + 95% CI (ADD vs each arm, BH-corrected) and "
        "a moderated ADD-vs-PDD DE contrast (cross-cohort-confounded). Exploratory; "
        "confound-bounded (findings 0001/0003/0004)."
    ),
}

CONDITION_COLUMN = "Condition"
ARMS = ("ADD", "PDD", "PDCN", "HCN")
SEED = 20260915
DATA_VERSION = "protein-data.csv sha256:90a0462eb03a02bd"
SIGNATURE_SET = "boruta:762e16f294e3 (40 Confirmed, ADD-vs-Other)"

DEFAULT_BORUTA_DIR = Path("results/boruta/762e16f294e3")
DEFAULT_DE_TABLE = Path(
    "results/differential_abundance/add-vs-other-moderated/contrast_table.csv"
)
DEFAULT_CACHE_ROOT = Path("results")
ANALYSIS_NAME = "four_arm"

# Per-arm cohort composition (state/METADATA.md, Table 1): Kerchner/Poston/ADRC-PUC/PUC.
ARM_COHORT = {
    "ADD": {"Kerchner": 10, "Poston": 0, "ADRC/PUC": 0, "PUC only": 0},
    "PDD": {"Kerchner": 0, "Poston": 7, "ADRC/PUC": 2, "PUC only": 1},
    "PDCN": {"Kerchner": 0, "Poston": 10, "ADRC/PUC": 0, "PUC only": 0},
    "HCN": {"Kerchner": 7, "Poston": 3, "ADRC/PUC": 0, "PUC only": 0},
}


@dataclass(frozen=True)
class FourArmResult:
    """Cacheable result bundle for the four-arm ADD-specificity check."""

    signature_direction: pd.DataFrame  # per confirmed feature: effect, direction
    sample_scores: pd.DataFrame  # per sample: arm, up_score, down_score
    per_arm_scores: pd.DataFrame  # per arm x set: mean, ci_low, ci_high, sd, n
    arm_comparisons: pd.DataFrame  # ADD vs each arm x set: g, ci, p, q
    add_vs_pdd_table: pd.DataFrame  # full 2340-feature ADD-vs-PDD contrast
    n_up: int
    n_down: int
    n_confirmed: int
    seed: int
    signature_set: str


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #
def confirmed_features(boruta_dir: Path) -> list[str]:
    """Return the Boruta-Confirmed feature ids (decision.json/feature_names.json)."""
    decision = json.loads((boruta_dir / "decision.json").read_text())["values"]
    names = json.loads((boruta_dir / "feature_names.json").read_text())["values"]
    if len(decision) != len(names):
        raise ValueError(
            f"decision ({len(decision)}) / names ({len(names)}) length mismatch "
            f"in {boruta_dir}."
        )
    conf = [n for n, d in zip(names, decision, strict=True) if d == "Confirmed"]
    if not conf:
        raise ValueError(f"No Confirmed features in {boruta_dir}.")
    logger.info("Boruta Confirmed: %d features", len(conf))
    return conf


def signature_direction(confirmed: list[str], de_table: Path) -> pd.DataFrame:
    """Attach the ADD-vs-Other DE log2FC + up/down direction to each confirmed feature.

    effect = log2FC(ADD - Other); effect > 0 -> up-in-ADD, effect < 0 -> down-in-ADD.
    Fails loud if any confirmed feature is absent from the DE table.
    """
    de = pd.read_csv(de_table)[["feature", "effect", "q"]]
    de = de.set_index("feature")
    missing = [f for f in confirmed if f not in de.index]
    if missing:
        raise ValueError(
            f"{len(missing)} Boruta-confirmed feature(s) absent from the DE table "
            f"{de_table}: {missing[:5]}"
        )
    rows = []
    for f in confirmed:
        eff = float(de.loc[f, "effect"])
        rows.append(
            {
                "feature": f,
                "effect_add_minus_other": eff,
                "add_vs_other_q": float(de.loc[f, "q"]),
                "direction": "up_in_ADD" if eff > 0 else "down_in_ADD",
            }
        )
    out = pd.DataFrame(rows)
    n_up = int((out["direction"] == "up_in_ADD").sum())
    n_down = int((out["direction"] == "down_in_ADD").sum())
    logger.info("signature direction: %d up-in-ADD, %d down-in-ADD", n_up, n_down)
    return out


# --------------------------------------------------------------------------- #
# Part 1: per-arm summary scores
# --------------------------------------------------------------------------- #
def _sample_scores(
    dataset: Dataset, direction: pd.DataFrame
) -> pd.DataFrame:
    """Per-sample mean-z summary score of the up-set and the down-set.

    z-score each signature protein across the 40 samples, then average the z over the
    up-in-ADD proteins and over the down-in-ADD proteins for each sample.
    """
    feat_index = {name: i for i, name in enumerate(dataset.feature_names.tolist())}
    X = dataset.abundances  # (40, n_features) log2
    z_full = (X - X.mean(axis=0, keepdims=True)) / X.std(axis=0, ddof=0, keepdims=True)

    up = direction.loc[direction["direction"] == "up_in_ADD", "feature"].tolist()
    down = direction.loc[direction["direction"] == "down_in_ADD", "feature"].tolist()
    up_cols = [feat_index[f] for f in up]
    down_cols = [feat_index[f] for f in down]

    up_score = z_full[:, up_cols].mean(axis=1)
    down_score = z_full[:, down_cols].mean(axis=1)
    arms = dataset.metadata[CONDITION_COLUMN].astype(str).to_numpy()
    return pd.DataFrame(
        {
            "sample": dataset.metadata["SampleID"].astype(str).to_numpy(),
            "arm": arms,
            "up_score": up_score,
            "down_score": down_score,
        }
    )


def _per_arm_scores(sample_scores: pd.DataFrame) -> pd.DataFrame:
    """Per-arm mean summary score + t-based 95% CI, for the up-set and the down-set."""
    rows = []
    for set_name in ("up_score", "down_score"):
        for arm in ARMS:
            vals = sample_scores.loc[sample_scores["arm"] == arm, set_name].to_numpy()
            n = vals.size
            mean = float(vals.mean())
            sd = float(vals.std(ddof=1))
            tcrit = float(stats.t.ppf(0.975, df=n - 1))
            half = tcrit * sd / np.sqrt(n)
            rows.append(
                {
                    "set": "up_in_ADD" if set_name == "up_score" else "down_in_ADD",
                    "arm": arm,
                    "mean_z": mean,
                    "ci_low": mean - half,
                    "ci_high": mean + half,
                    "sd": sd,
                    "n": n,
                }
            )
    return pd.DataFrame(rows)


def _hedges_g(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float]:
    """Hedges' g (a vs b) + approximate 95% CI. Positive g -> a > b."""
    n1, n2 = a.size, b.size
    s1, s2 = a.std(ddof=1), b.std(ddof=1)
    sp = np.sqrt(((n1 - 1) * s1**2 + (n2 - 1) * s2**2) / (n1 + n2 - 2))
    d = (a.mean() - b.mean()) / sp
    J = 1.0 - 3.0 / (4.0 * (n1 + n2) - 9.0)
    g = J * d
    se = np.sqrt((n1 + n2) / (n1 * n2) + g**2 / (2.0 * (n1 + n2 - 2)))
    return float(g), float(g - 1.96 * se), float(g + 1.96 * se)


def _bh(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(pvals, dtype=float)
    m = p.size
    order = np.argsort(p)
    ranked = p[order] * m / (np.arange(m) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(ranked, 0, 1)
    return out


def _arm_comparisons(sample_scores: pd.DataFrame) -> pd.DataFrame:
    """ADD vs each other arm, up-set and down-set: Hedges' g + CI + Welch p + BH-q.

    BH is applied over the full family of 6 comparisons (2 sets x 3 other arms).
    """
    rows = []
    for set_name in ("up_score", "down_score"):
        add = sample_scores.loc[sample_scores["arm"] == "ADD", set_name].to_numpy()
        for arm in ("PDD", "PDCN", "HCN"):
            other = sample_scores.loc[sample_scores["arm"] == arm, set_name].to_numpy()
            g, lo, hi = _hedges_g(add, other)
            t, p = stats.ttest_ind(add, other, equal_var=False)
            rows.append(
                {
                    "set": "up_in_ADD" if set_name == "up_score" else "down_in_ADD",
                    "comparison": f"ADD_vs_{arm}",
                    "mean_diff_add_minus_arm": float(add.mean() - other.mean()),
                    "hedges_g": g,
                    "g_ci_low": lo,
                    "g_ci_high": hi,
                    "welch_t": float(t),
                    "p": float(p),
                }
            )
    out = pd.DataFrame(rows)
    out["q_bh"] = _bh(out["p"].to_numpy())
    return out


# --------------------------------------------------------------------------- #
# Part 2: ADD vs PDD moderated DE (cross-cohort-confounded)
# --------------------------------------------------------------------------- #
def add_vs_pdd_de(dataset: Dataset) -> pd.DataFrame:
    """Moderated (limma eBayes) ADD-vs-PDD contrast on the ADD+PDD subset (20 samples).

    Subsets to the two dementia arms, builds a binary factor (reference=PDD so effect =
    log2FC ADD - PDD), and runs the moderated model. Returns the full 2340-feature
    contrast table. CROSS-COHORT-CONFOUNDED (ADD 100% Kerchner vs PDD ~70% Poston).
    """
    meta = dataset.metadata
    mask = meta[CONDITION_COLUMN].astype(str).isin(["ADD", "PDD"]).to_numpy()
    n_add = int((meta[CONDITION_COLUMN].astype(str) == "ADD").sum())
    n_pdd = int((meta[CONDITION_COLUMN].astype(str) == "PDD").sum())
    if (n_add, n_pdd) != (10, 10):
        raise ValueError(f"expected 10 ADD / 10 PDD; got {n_add}/{n_pdd}.")
    sub_meta = meta.loc[mask].reset_index(drop=True).copy()
    sub_meta["ADDvsPDD"] = sub_meta[CONDITION_COLUMN].astype(str)
    sub = dataclasses.replace(
        dataset,
        abundances=dataset.abundances[mask, :],
        metadata=sub_meta,
    )
    result = differential_abundance(
        sub,
        contrast="ADDvsPDD",
        method="moderated",
        reference={"ADDvsPDD": "PDD"},  # effect = ADD - PDD
        covariates=(),
    )
    logger.info(
        "ADD-vs-PDD moderated DE: %d features, %d samples",
        len(result.contrast_table),
        result.n_samples,
    )
    return result.contrast_table


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def run_four_arm(
    dataset: Dataset, boruta_dir: Path, de_table: Path
) -> FourArmResult:
    confirmed = confirmed_features(boruta_dir)
    direction = signature_direction(confirmed, de_table)
    sample_scores = _sample_scores(dataset, direction)
    per_arm = _per_arm_scores(sample_scores)
    comparisons = _arm_comparisons(sample_scores)
    add_pdd = add_vs_pdd_de(dataset)

    n_up = int((direction["direction"] == "up_in_ADD").sum())
    n_down = int((direction["direction"] == "down_in_ADD").sum())
    return FourArmResult(
        signature_direction=direction,
        sample_scores=sample_scores,
        per_arm_scores=per_arm,
        arm_comparisons=comparisons,
        add_vs_pdd_table=add_pdd,
        n_up=n_up,
        n_down=n_down,
        n_confirmed=len(confirmed),
        seed=SEED,
        signature_set=SIGNATURE_SET,
    )


def _overlap_report(result: FourArmResult) -> dict[str, object]:
    """ADD-vs-PDD hit counts + overlap with the Boruta-40 + direction concordance."""
    de = result.add_vs_pdd_table.set_index("feature")
    q = de["q"].to_numpy(dtype=float)
    n_q05 = int(np.nansum(q < 0.05))
    n_q10 = int(np.nansum(q < 0.10))

    conf = result.signature_direction.copy()
    conf = conf.merge(
        de[["effect", "q"]].rename(
            columns={"effect": "add_vs_pdd_effect", "q": "add_vs_pdd_q"}
        ),
        left_on="feature",
        right_index=True,
        how="left",
    )
    conf_q10 = int((conf["add_vs_pdd_q"] < 0.10).sum())
    conf_q05 = int((conf["add_vs_pdd_q"] < 0.05).sum())
    # Direction concordance among the Boruta-40: does ADD-vs-PDD point the same way as
    # ADD-vs-Other?
    same_sign = int(
        (
            np.sign(conf["effect_add_minus_other"])
            == np.sign(conf["add_vs_pdd_effect"])
        ).sum()
    )
    return {
        "add_vs_pdd_n_bh_q_lt_0.05": n_q05,
        "add_vs_pdd_n_bh_q_lt_0.10": n_q10,
        "boruta40_add_vs_pdd_q_lt_0.05": conf_q05,
        "boruta40_add_vs_pdd_q_lt_0.10": conf_q10,
        "boruta40_direction_concordant_add_other_vs_add_pdd": same_sign,
        "boruta40_n": len(conf),
    }


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("data/protein-data.csv"))
    parser.add_argument("--metadata", type=Path, default=Path("data/metadata.csv"))
    parser.add_argument("--boruta-dir", type=Path, default=DEFAULT_BORUTA_DIR)
    parser.add_argument("--de-table", type=Path, default=DEFAULT_DE_TABLE)
    parser.add_argument("--cache-root", type=Path, default=DEFAULT_CACHE_ROOT)
    args = parser.parse_args(argv)

    dataset = load_analysis_subset(load_paired_dataset(args.data, args.metadata))
    logger.info(
        "loaded analysis subset: %d samples x %d features, scale=%s",
        dataset.abundances.shape[0],
        dataset.abundances.shape[1],
        dataset.scale,
    )

    result = run_four_arm(dataset, args.boruta_dir, args.de_table)
    overlap = _overlap_report(result)

    params = {
        "signature_set": SIGNATURE_SET,
        "arms": list(ARMS),
        "score": "per-protein z across 40 samples, per-arm mean of up/down subset",
        "effect_size": "Hedges g (ADD vs arm), 95% CI",
        "comparison_test": "Welch t",
        "correction": "BH over 6 comparisons (2 sets x 3 arms)",
        "add_vs_pdd": "moderated eBayes DE, effect=ADD-PDD, BH-q, CROSS-COHORT",
        "method": "z-summary-score + moderated-DE",
    }
    meta = save_cached_result(
        result,
        cache_root=args.cache_root,
        analysis=ANALYSIS_NAME,
        data_version=DATA_VERSION,
        params=params,
        seed=SEED,
        label="four-arm ADD-specificity (Boruta-40 up/down profile + ADD-vs-PDD DE)",
    )
    cache_dir = args.cache_root / ANALYSIS_NAME / meta.fingerprint
    (cache_dir / "overlap_report.json").write_text(json.dumps(overlap, indent=2))

    logger.info("CACHE: %s (fingerprint %s)", cache_dir, meta.fingerprint)
    logger.info("PER-ARM SCORES:\n%s", result.per_arm_scores.to_string(index=False))
    logger.info(
        "ARM COMPARISONS (ADD vs arm):\n%s",
        result.arm_comparisons.to_string(index=False),
    )
    logger.info("OVERLAP / ADD-vs-PDD: %s", json.dumps(overlap, indent=2))


if __name__ == "__main__":
    main()
