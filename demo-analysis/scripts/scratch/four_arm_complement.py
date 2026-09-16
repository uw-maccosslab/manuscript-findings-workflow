"""FOUR-ARM COMPLEMENT/COAGULATION: ADD-specific program, or Kerchner site artifact?

Finding 0011's headline up-in-ADD GSEA signal is the KEGG "Complement and coagulation
cascades" program -- which is ALSO a textbook plasma pre-analytical / collection
artifact. ADD is 100% Kerchner cohort, so an up-in-ADD complement/coagulation shift is
confounded with the collection site. The four-arm design tests this directly, and HCN is
the discriminating arm: 70% Kerchner (7/10) but cognitively normal.

Primary set: the 27 leading-edge genes of KEGG "Complement and coagulation cascades"
(the Lead_genes of that row in results/gsea/fed7a349c550/gsea_kegg.csv -- the
finding-0011 up-in-ADD program). Mapped to matrix feature-ids via ranked_symbols.csv.

Readouts (all EXPLORATORY -- n=10/arm; confound-bounded; findings 0001/0003/0004/0006):

  1. Per-arm profile. z-score each mapped complement/coagulation protein across the 40
     samples, then per SAMPLE take the mean z over the set -> one summary score/sample.
     Report each arm's (ADD/PDD/PDCN/HCN) mean + t 95% CI; test ADD vs each arm (Welch)
     with Hedges' g + 95% CI, BH over the 3 comparisons. KEY non-circular reads: where
     HCN (majority-Kerchner) and PDD sit relative to ADD and to BASELINE (z=0).
       - HCN at baseline despite majority-Kerchner  => weak evidence AGAINST a pure
         Kerchner pre-analytical artifact (a pure site artifact should elevate HCN too).
       - HCN elevated toward ADD                     => consistent with a site artifact.

  2. Within-HCN site test. Within HCN only (7 Kerchner vs 3 Poston), is the program
     higher in Kerchner? Hedges' g + CI. A direct (tiny-n) probe of site vs disease.
     Same within PDD (7 Poston / 2 ADRC-PUC / 1 PUC-only) if the cells allow.

  Secondary/robustness: the GO:BP "regulation of complement activation" + "regulation of
  humoral immune response" leading-edge UNION (15 genes each, identical here) as an
  alternative program definition.

The leading-edge set is defined FROM the ADD-vs-Other ranking, so ADD is extreme BY
CONSTRUCTION; the non-circular readouts are (a) the PDD/PDCN/HCN positions and (b) the
within-HCN site test.

Deterministic; seed 20260915 (no RNG -- recorded for provenance/fingerprint). Fail-loud.
Caches result tables + meta.json under results/four_arm_complement/<fingerprint>/.

Run::

    ./.venv/bin/python scripts/scratch/four_arm_complement.py
"""

from __future__ import annotations

import argparse
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
from project_loader import load_analysis_subset, load_paired_dataset

logger = logging.getLogger(__name__)

__script_meta__: dict[str, object] = {
    "task": "four-arm-complement-coagulation-site-vs-add",
    "kind": "analysis",
    "provides": [
        "ComplementResult",
        "map_symbols_to_features",
        "run_complement",
        "main",
    ],
    "uses": ["project_loader", "common.data_loading", "common.result_io"],
    "seeded_from": {"template": "four_arm_specificity.py", "version": "0.1"},
    "description": (
        "Four-arm complement/coagulation program check: is the KEGG complement "
        "leading-edge (up-in-ADD, finding 0011) ADD-specific or does it track the "
        "Kerchner collection site? Per-arm z-summary score + ADD-vs-arm Hedges g (BH), "
        "a within-HCN Kerchner-vs-Poston probe. Exploratory; confound-bounded."
    ),
}

CONDITION_COLUMN = "Condition"
COHORT_COLUMN = "Cohort"
ARMS = ("ADD", "PDD", "PDCN", "HCN")
SEED = 20260915
DATA_VERSION = "protein-data.csv sha256:90a0462eb03a02bd"
SOURCE_GSEA = "gsea:fed7a349c550 (KEGG Complement/coagulation cascades leading-edge)"
ANALYSIS_NAME = "four_arm_complement"

DEFAULT_GSEA_DIR = Path("results/gsea/fed7a349c550")
DEFAULT_CACHE_ROOT = Path("results")

# KEGG "Complement and coagulation cascades" leading-edge (27), finding 0011 up-in-ADD.
KEGG_COMPLEMENT_LEAD = [
    "SERPIND1", "CPB2", "KNG1", "PROS1", "PROC", "C4A", "FGA", "F2", "SERPING1",
    "C5", "C3", "CFB", "CFH", "C8A", "SERPINC1", "C8B", "CFI", "VTN", "SERPINA1",
    "A2M", "C1QB", "F9", "ITGB2", "SERPINE1", "FGB", "F5", "MASP1",
]
# GO:BP humoral-immune union (secondary): reg. of complement activation + reg. of
# humoral immune response leading-edge (15 each; identical membership here).
GOBP_HUMORAL_LEAD = [
    "CPB2", "PROS1", "C4A", "F2", "SERPING1", "C5", "C3", "CFB", "CFH", "C8A",
    "C8B", "CFI", "VTN", "C1QB", "CPN1",
]


@dataclass(frozen=True)
class ComplementResult:
    """Cacheable bundle for the complement/coagulation four-arm site check."""

    mapping: pd.DataFrame  # symbol -> accession -> feature, in_matrix (primary set)
    mapping_secondary: pd.DataFrame  # secondary GO:BP humoral set mapping
    sample_scores: pd.DataFrame  # per sample: arm, cohort, primary + secondary score
    per_arm_scores: pd.DataFrame  # per arm: mean, ci, sd, n  (primary + secondary)
    arm_comparisons: pd.DataFrame  # ADD vs arm: g, ci, welch p, BH-q (both sets)
    within_site_tests: pd.DataFrame  # within-arm Kerchner-vs-other site probes
    cohort_orthogonality: pd.DataFrame  # per-arm max/mean |per-protein cohort diff|
    n_primary_mapped: int
    n_primary_requested: int
    n_secondary_mapped: int
    seed: int
    source_gsea: str


# --------------------------------------------------------------------------- #
# Symbol -> matrix-feature mapping
# --------------------------------------------------------------------------- #
def map_symbols_to_features(
    symbols: list[str], ranked_symbols_csv: Path, matrix_features: set[str]
) -> pd.DataFrame:
    """Map gene symbols to matrix feature-ids via the GSEA ranked_symbols table.

    ranked_symbols.csv columns: symbol, accession, feature, statistic, effect. The
    ``feature`` field is the matrix feature-id (e.g. ``sp|P01023|A2MG_HUMAN``). Returns
    one row per requested symbol with the resolved feature and an ``in_matrix`` flag.
    """
    ranked = pd.read_csv(ranked_symbols_csv)
    by_symbol: dict[str, tuple[str, str]] = {}
    for _, r in ranked.iterrows():
        by_symbol.setdefault(str(r["symbol"]), (str(r["accession"]), str(r["feature"])))
    rows = []
    for sym in symbols:
        if sym in by_symbol:
            acc, feat = by_symbol[sym]
            in_mat = feat in matrix_features
            rows.append({"symbol": sym, "accession": acc, "feature": feat,
                         "in_ranking": True, "in_matrix": in_mat})
        else:
            rows.append({"symbol": sym, "accession": "", "feature": "",
                         "in_ranking": False, "in_matrix": False})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Summary score + stats
# --------------------------------------------------------------------------- #
def _score_from_features(dataset: Dataset, features: list[str]) -> np.ndarray:
    """Per-sample mean-z summary score over ``features`` (z across the 40 samples)."""
    feat_index = {name: i for i, name in enumerate(dataset.feature_names.tolist())}
    cols = [feat_index[f] for f in features]
    X = dataset.abundances[:, cols]  # (40, k) log2
    z = (X - X.mean(axis=0, keepdims=True)) / X.std(axis=0, ddof=0, keepdims=True)
    return np.asarray(z.mean(axis=1), dtype=float)


def _per_arm(sample_scores: pd.DataFrame, col: str, label: str) -> pd.DataFrame:
    rows = []
    for arm in ARMS:
        vals = sample_scores.loc[sample_scores["arm"] == arm, col].to_numpy()
        n = vals.size
        mean = float(vals.mean())
        sd = float(vals.std(ddof=1))
        tcrit = float(stats.t.ppf(0.975, df=n - 1))
        half = tcrit * sd / np.sqrt(n)
        rows.append({"set": label, "arm": arm, "mean_z": mean,
                     "ci_low": mean - half, "ci_high": mean + half, "sd": sd, "n": n})
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
    p = np.asarray(pvals, dtype=float)
    m = p.size
    order = np.argsort(p)
    ranked = p[order] * m / (np.arange(m) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(ranked, 0, 1)
    return out


def _arm_comparisons(sample_scores: pd.DataFrame, col: str, label: str) -> pd.DataFrame:
    """ADD vs each other arm on ``col``: Hedges' g + CI + Welch p + BH-q (3 tests)."""
    rows = []
    add = sample_scores.loc[sample_scores["arm"] == "ADD", col].to_numpy()
    for arm in ("PDD", "PDCN", "HCN"):
        other = sample_scores.loc[sample_scores["arm"] == arm, col].to_numpy()
        g, lo, hi = _hedges_g(add, other)
        t, p = stats.ttest_ind(add, other, equal_var=False)
        rows.append({"set": label, "comparison": f"ADD_vs_{arm}",
                     "mean_diff_add_minus_arm": float(add.mean() - other.mean()),
                     "hedges_g": g, "g_ci_low": lo, "g_ci_high": hi,
                     "welch_t": float(t), "p": float(p)})
    out = pd.DataFrame(rows)
    out["q_bh"] = _bh(out["p"].to_numpy())
    return out


def _within_site_tests(
    sample_scores: pd.DataFrame, col: str, label: str
) -> pd.DataFrame:
    """Within-arm site probes: majority cohort vs the rest, per arm where estimable.

    HCN: Kerchner (7) vs Poston (3). PDD: Poston (7) vs non-Poston (2 ADRC/PUC + 1 PUC).
    Reports Hedges' g + CI + Welch p; flags tiny-n. Skips a test if either group < 2.
    """
    rows = []
    plans = [
        ("HCN", "Kerchner", ["Poston"]),
        ("PDD", "Poston", ["ADRC/PUC", "PUC only"]),
    ]
    for arm, site_a, sites_b in plans:
        sub = sample_scores[sample_scores["arm"] == arm]
        a = sub.loc[sub["cohort"] == site_a, col].to_numpy()
        b = sub.loc[sub["cohort"].isin(sites_b), col].to_numpy()
        if a.size < 2 or b.size < 2:
            rows.append({"set": label, "arm": arm,
                         "site_a": site_a, "n_a": int(a.size),
                         "site_b": "+".join(sites_b), "n_b": int(b.size),
                         "mean_a": float(a.mean()) if a.size else float("nan"),
                         "mean_b": float(b.mean()) if b.size else float("nan"),
                         "mean_diff_a_minus_b": float("nan"),
                         "hedges_g": float("nan"), "g_ci_low": float("nan"),
                         "g_ci_high": float("nan"), "welch_p": float("nan"),
                         "estimable": False})
            continue
        g, lo, hi = _hedges_g(a, b)
        _, p = stats.ttest_ind(a, b, equal_var=False)
        rows.append({"set": label, "arm": arm,
                     "site_a": site_a, "n_a": int(a.size),
                     "site_b": "+".join(sites_b), "n_b": int(b.size),
                     "mean_a": float(a.mean()), "mean_b": float(b.mean()),
                     "mean_diff_a_minus_b": float(a.mean() - b.mean()),
                     "hedges_g": g, "g_ci_low": lo, "g_ci_high": hi,
                     "welch_p": float(p), "estimable": True})
    return pd.DataFrame(rows)


def _cohort_orthogonality(dataset: Dataset) -> pd.DataFrame:
    """Per-arm max/mean |per-protein cohort-group mean difference| over ALL features.

    A within-arm site test can only detect a site effect if cohort actually imprints on
    the abundances within that arm. This diagnostic measures that imprint directly on
    the full 2340-feature matrix (raw log2, no z), so the cache records WHETHER the
    within-arm site probe is informative here or degenerate (~0 imprint => useless).
    """
    meta = dataset.metadata
    arm = meta[CONDITION_COLUMN].astype(str).to_numpy()
    coh = meta[COHORT_COLUMN].astype(str).to_numpy()
    X = dataset.abundances
    plans = [("HCN", "Kerchner", ["Poston"]),
             ("PDD", "Poston", ["ADRC/PUC", "PUC only"])]
    rows = []
    for a, site_a, sites_b in plans:
        ma = (arm == a) & (coh == site_a)
        mb = (arm == a) & np.isin(coh, sites_b)
        if ma.sum() < 2 or mb.sum() < 2:
            rows.append({"arm": a, "site_a": site_a, "site_b": "+".join(sites_b),
                         "n_a": int(ma.sum()), "n_b": int(mb.sum()),
                         "max_abs_protein_mean_diff": float("nan"),
                         "mean_abs_protein_mean_diff": float("nan"),
                         "informative": False})
            continue
        d = np.abs(X[ma].mean(0) - X[mb].mean(0))
        rows.append({"arm": a, "site_a": site_a, "site_b": "+".join(sites_b),
                     "n_a": int(ma.sum()), "n_b": int(mb.sum()),
                     "max_abs_protein_mean_diff": float(d.max()),
                     "mean_abs_protein_mean_diff": float(d.mean()),
                     "informative": bool(d.max() > 1e-6)})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def run_complement(dataset: Dataset, gsea_dir: Path) -> ComplementResult:
    matrix_features = set(dataset.feature_names.tolist())
    ranked_csv = gsea_dir / "ranked_symbols.csv"

    map_primary = map_symbols_to_features(
        KEGG_COMPLEMENT_LEAD, ranked_csv, matrix_features
    )
    map_secondary = map_symbols_to_features(
        GOBP_HUMORAL_LEAD, ranked_csv, matrix_features
    )

    feats_primary = map_primary.loc[map_primary["in_matrix"], "feature"].tolist()
    feats_secondary = map_secondary.loc[map_secondary["in_matrix"], "feature"].tolist()
    if not feats_primary:
        raise ValueError("No primary complement/coagulation feature mapped in matrix.")
    logger.info("primary mapped %d/%d; secondary mapped %d/%d",
                len(feats_primary), len(KEGG_COMPLEMENT_LEAD),
                len(feats_secondary), len(GOBP_HUMORAL_LEAD))

    score_primary = _score_from_features(dataset, feats_primary)
    score_secondary = _score_from_features(dataset, feats_secondary)
    sample_scores = pd.DataFrame({
        "sample": dataset.metadata["SampleID"].astype(str).to_numpy(),
        "arm": dataset.metadata[CONDITION_COLUMN].astype(str).to_numpy(),
        "cohort": dataset.metadata[COHORT_COLUMN].astype(str).to_numpy(),
        "complement_score": score_primary,
        "humoral_score": score_secondary,
    })

    per_arm = pd.concat([
        _per_arm(sample_scores, "complement_score", "KEGG_complement_coag"),
        _per_arm(sample_scores, "humoral_score", "GOBP_humoral_union"),
    ], ignore_index=True)
    comparisons = pd.concat([
        _arm_comparisons(sample_scores, "complement_score", "KEGG_complement_coag"),
        _arm_comparisons(sample_scores, "humoral_score", "GOBP_humoral_union"),
    ], ignore_index=True)
    within = pd.concat([
        _within_site_tests(sample_scores, "complement_score", "KEGG_complement_coag"),
        _within_site_tests(sample_scores, "humoral_score", "GOBP_humoral_union"),
    ], ignore_index=True)
    orthogonality = _cohort_orthogonality(dataset)

    return ComplementResult(
        mapping=map_primary,
        mapping_secondary=map_secondary,
        sample_scores=sample_scores,
        per_arm_scores=per_arm,
        arm_comparisons=comparisons,
        within_site_tests=within,
        cohort_orthogonality=orthogonality,
        n_primary_mapped=len(feats_primary),
        n_primary_requested=len(KEGG_COMPLEMENT_LEAD),
        n_secondary_mapped=len(feats_secondary),
        seed=SEED,
        source_gsea=SOURCE_GSEA,
    )


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("data/protein-data.csv"))
    parser.add_argument("--metadata", type=Path, default=Path("data/metadata.csv"))
    parser.add_argument("--gsea-dir", type=Path, default=DEFAULT_GSEA_DIR)
    parser.add_argument("--cache-root", type=Path, default=DEFAULT_CACHE_ROOT)
    args = parser.parse_args(argv)

    dataset = load_analysis_subset(load_paired_dataset(args.data, args.metadata))
    logger.info("loaded analysis subset: %d x %d, scale=%s",
                dataset.abundances.shape[0], dataset.abundances.shape[1], dataset.scale)

    result = run_complement(dataset, args.gsea_dir)

    params = {
        "primary_set": SOURCE_GSEA,
        "primary_symbols": KEGG_COMPLEMENT_LEAD,
        "secondary_set": "gsea:fed7a349c550 GO:BP reg-complement + reg-humoral union",
        "secondary_symbols": GOBP_HUMORAL_LEAD,
        "arms": list(ARMS),
        "score": "per-protein z across 40 samples, per-sample mean z over mapped set",
        "effect_size": "Hedges g (ADD vs arm), 95% CI",
        "comparison_test": "Welch t",
        "correction": "BH over 3 ADD-vs-arm comparisons (per set)",
        "within_site_probe": "HCN Kerchner-vs-Poston; PDD Poston-vs-rest (tiny-n)",
        "method": "z-summary-score four-arm profile + within-arm site test",
    }
    meta = save_cached_result(
        result,
        cache_root=args.cache_root,
        analysis=ANALYSIS_NAME,
        data_version=DATA_VERSION,
        params=params,
        seed=SEED,
        label="four-arm complement/coagulation: ADD-specific vs Kerchner site artifact",
    )
    cache_dir = args.cache_root / ANALYSIS_NAME / meta.fingerprint
    summary = {
        "n_primary_mapped": result.n_primary_mapped,
        "n_primary_requested": result.n_primary_requested,
        "n_secondary_mapped": result.n_secondary_mapped,
        "unmapped_primary_symbols": result.mapping.loc[
            ~result.mapping["in_matrix"], "symbol"].tolist(),
    }
    (cache_dir / "map_summary.json").write_text(json.dumps(summary, indent=2))

    logger.info("CACHE: %s (fingerprint %s)", cache_dir, meta.fingerprint)
    logger.info("MAP SUMMARY: %s", json.dumps(summary))
    logger.info("PER-ARM:\n%s", result.per_arm_scores.to_string(index=False))
    logger.info("ADD-vs-ARM:\n%s", result.arm_comparisons.to_string(index=False))
    logger.info("WITHIN-SITE:\n%s", result.within_site_tests.to_string(index=False))
    logger.info("COHORT-ORTHOGONALITY (within-arm site imprint on full matrix):\n%s",
                result.cohort_orthogonality.to_string(index=False))


if __name__ == "__main__":
    main()
