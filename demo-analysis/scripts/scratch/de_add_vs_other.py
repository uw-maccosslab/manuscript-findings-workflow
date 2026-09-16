"""PRIMARY differential-abundance analysis: ADD vs all-others (pooled), plasma-EV study.

Call site for the project copy of the ``differential-abundance`` lib/ template
(``scripts/scratch/differential_abundance.py`` — seeded from
``lib/analysis/differential-abundance`` @ 0.1). It:

  1. loads the 40-sample x 2340-protein log2 analysis subset via ``project_loader``
     (controls already excluded; yeast spike-in + two floored-sentinel features dropped
     upstream — state/DATA_DESCRIPTION.md);
  2. builds the binary primary contrast **ADD (n=10) vs Other (n=30)**, where
     ``Other = {PDD, PDCN, HCN}`` pooled, as a new metadata factor;
  3. runs ``differential_abundance(..., method="moderated")`` (limma-style empirical
     Bayes) with ``reference="Other"`` so the effect is **log2FC = ADD - Other**;
  4. writes the full per-feature contrast table (effect + 95% CI + p + BH-q +
     mean abundance) to ``results/differential_abundance/<label>/`` and a JSON summary
     carrying the p-value-histogram calibration read + hit counts.

CRITICAL CAVEATS baked into the reporting (state/DATA_DESCRIPTION.md; findings 0001 +
0004):

  * **No Cohort covariate.** ADD is 100% Kerchner (single-cohort), so Cohort is aliased
    with the contrast (rank-deficient) — adding it would make the design singular. The
    contrast is therefore reported WITHOUT cohort adjustment; the ADD effect is
    inseparable from a Kerchner-site effect. This is a hard design confound, not a knob.
  * **Matrix is already limma-batch-corrected on Cohort** (removeBatchEffect protecting
    Condition). The pre-correction matrix is unavailable, so the ideal
    batch-as-covariate-on-uncorrected-data robustness check is IMPOSSIBLE. DE on
    removeBatchEffect output DEFLATES within-feature variance -> p-values are
    ANTICONSERVATIVE. Emphasize effect sizes + CIs; treat significance cautiously.
  * **Exploratory.** Small n (10 vs 30), imbalanced, confounded. Not confirmatory.

No label-shuffle null is run here (that is the classifier/regressor discipline, not the
univariate DE family); the moderated model's own eBayes moderation is the variance
control. The honesty check available without a null is the **p-value histogram
calibration read**, computed and reported below.

Non-interactive, parameterized (argparse), seeded where relevant (the analysis is
deterministic — no RNG — but the run stamp is recorded). Fail-loud at every boundary.

Run::

    ./.venv/bin/python scripts/scratch/de_add_vs_other.py
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Make `scripts/` and `scripts/scratch/` importable regardless of invocation (direct
# script / -m / pytest), mirroring project_loader.py's pattern. (ruff exempts bare
# `sys.path` mutations from E402.)
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from common.data_loading import Dataset
from differential_abundance import (
    DifferentialAbundanceResult,
    differential_abundance,
)
from project_loader import load_analysis_subset, load_paired_dataset

logger = logging.getLogger(__name__)

__script_meta__: dict[str, object] = {
    "task": "de-add-vs-other",
    "kind": "analysis",
    "provides": [
        "CONTRAST_COLUMN",
        "ADD_LEVEL",
        "OTHER_LEVEL",
        "build_binary_contrast",
        "run_primary_de",
        "pvalue_calibration_read",
        "annotate_known_markers",
        "main",
    ],
    "uses": ["project_loader", "differential_abundance", "common.data_loading"],
    "seeded_from": {"template": "differential-abundance", "version": "0.1"},
    "description": (
        "PRIMARY plasma-EV differential abundance: ADD (10) vs pooled Other "
        "{PDD,PDCN,HCN} (30), moderated (limma eBayes) linear model on the delivered "
        "log2 40x2340 analysis subset, effect = log2FC(ADD-Other) + 95% CI + p + BH-q. "
        "NO cohort covariate (ADD is single-cohort -> aliased with contrast, findings "
        "0001/0004); matrix already batch-corrected -> p-values anticonservative; "
        "exploratory. Writes the full contrast CSV + a JSON calibration/summary read."
    ),
}

CONTRAST_COLUMN = "ADDvsOther"
ADD_LEVEL = "ADD"
OTHER_LEVEL = "Other"
CONDITION_COLUMN = "Condition"
OTHER_CONDITIONS: frozenset[str] = frozenset({"PDD", "PDCN", "HCN"})

DEFAULT_OUTPUT_DIR = Path("results/differential_abundance")
DEFAULT_LABEL = "add-vs-other-moderated"

# Curated, OUTCOME-INDEPENDENT reference markers (research/ad-plasma-ev-proteomic-
# markers.md + research/ad-plasma-*-biomarkers.md), keyed by UniProt accession. Used to
# ANNOTATE hits for interpretation — never to select or filter features.
KNOWN_MARKERS: dict[str, tuple[str, str]] = {
    # MISEV2018 canonical EV-quality-control markers (identity/recovery, NOT disease).
    "P21926": ("CD9", "EV-QC (MISEV2018 tetraspanin)"),
    "P08962": ("CD63", "EV-QC (MISEV2018 tetraspanin)"),
    "P60033": ("CD81", "EV-QC (MISEV2018 tetraspanin)"),
    "O00560": ("SDCBP/syntenin-1", "EV-QC (MISEV2018 cytosolic)"),
    "Q8WUM4": ("PDCD6IP/ALIX", "EV-QC (MISEV2018 cytosolic)"),
    "Q99816": ("TSG101", "EV-QC (MISEV2018 cytosolic)"),
    "O75955": ("FLOT1", "EV-QC (MISEV2018 flotillin)"),
    # Contested neuron-derived-EV capture markers.
    "P32004": ("L1CAM", "contested NDEV capture marker"),
    "P13591": ("NCAM1", "contested NDEV capture marker"),
    # Unbiased plasma-EV MS AD candidates (single small discovery study; hypothesis).
    "P00488": ("F13A1", "unbiased plasma-EV AD candidate (Nielsen 2021, unreplicated)"),
    "P05160": ("F13B", "unbiased plasma-EV AD candidate (Nielsen 2021, unreplicated)"),
    # AD-marker parent proteins detected in QC.
    "P14136": ("GFAP", "AD glial marker (parent protein)"),
    "P05067": ("APP", "AD amyloid parent protein"),
}


def build_binary_contrast(dataset: Dataset) -> Dataset:
    """Return a copy of ``dataset`` with the binary ``ADDvsOther`` factor added.

    ``ADD`` stays ``ADD``; ``{PDD, PDCN, HCN}`` collapse to ``Other``. Fails loud if
    the ``Condition`` column is absent, carries an unexpected level, or the resulting
    split is not the expected 10 ADD / 30 Other (a data change must not slip through).
    """
    metadata = dataset.metadata
    if CONDITION_COLUMN not in metadata.columns:
        raise ValueError(
            f"metadata is missing required column {CONDITION_COLUMN!r}; cannot build "
            f"the ADD-vs-Other contrast."
        )
    conditions = metadata[CONDITION_COLUMN].astype(str)
    expected_levels = {ADD_LEVEL, *OTHER_CONDITIONS}
    seen = set(conditions.unique())
    unexpected = seen - expected_levels
    if unexpected:
        raise ValueError(
            f"Unexpected {CONDITION_COLUMN!r} level(s) {sorted(unexpected)} in the "
            f"analysis subset; expected exactly {sorted(expected_levels)} "
            f"(controls should already be excluded upstream)."
        )
    labels = np.where(conditions.to_numpy() == ADD_LEVEL, ADD_LEVEL, OTHER_LEVEL)
    new_metadata = metadata.copy()
    new_metadata[CONTRAST_COLUMN] = labels

    n_add = int((labels == ADD_LEVEL).sum())
    n_other = int((labels == OTHER_LEVEL).sum())
    if (n_add, n_other) != (10, 30):
        raise ValueError(
            f"ADD-vs-Other split is {n_add} ADD / {n_other} Other; expected 10 / 30. "
            f"Confirm the analysis subset and Condition labels (state/METADATA.md)."
        )
    logger.info(
        "built %r factor: %d %s vs %d %s",
        CONTRAST_COLUMN,
        n_add,
        ADD_LEVEL,
        n_other,
        OTHER_LEVEL,
    )
    return dataclasses.replace(dataset, metadata=new_metadata)


def run_primary_de(dataset: Dataset) -> DifferentialAbundanceResult:
    """Run the moderated ADD-vs-Other contrast (effect sign = ADD - Other).

    NO covariates: Cohort is aliased with the contrast (ADD is single-cohort; findings
    0001/0004), so a Cohort term would make the design singular; every other covariate
    is left out deliberately to keep the primary contrast the plain group difference.
    """
    with_contrast = build_binary_contrast(dataset)
    result = differential_abundance(
        with_contrast,
        contrast=CONTRAST_COLUMN,
        method="moderated",
        reference={
            CONTRAST_COLUMN: OTHER_LEVEL
        },  # effect = non-ref - ref = ADD - Other
        covariates=(),
    )
    logger.info(
        "moderated DE done: %d features, prior_var=%.4g prior_df=%.4g, term=%s",
        len(result.contrast_table),
        result.prior_variance if result.prior_variance is not None else float("nan"),
        result.prior_df if result.prior_df is not None else float("nan"),
        result.contrast_terms,
    )
    return result


def pvalue_calibration_read(pvalues: np.ndarray) -> dict[str, object]:
    """Summarize the contrast p-value distribution as a calibration honesty check.

    Reports the 20-bin histogram, the null-region density (a Storey pi0 estimate from
    the p>0.5 mass), the near-zero spike ratio (bin[0] density / null density), and a
    verdict string: well-calibrated (uniform + modest spike), anti-conservative
    (null-region density > 1, i.e. inflated small-p mass beyond the top bin), or a
    strong near-zero spike.
    """
    p = np.asarray(pvalues, dtype=float)
    p = p[np.isfinite(p)]
    n = int(p.size)
    if n == 0:
        raise ValueError("no finite p-values to calibrate.")
    edges = np.linspace(0.0, 1.0, 21)
    counts, _ = np.histogram(p, bins=edges)
    density = counts / (n / 20.0)  # expected count per bin under the uniform null

    # Storey-style pi0 from the upper half (p > 0.5), where true nulls dominate.
    tail_mass = float((p > 0.5).sum())
    pi0 = min(1.0, (tail_mass / (0.5 * n)))
    null_density = float(density[10:].mean())  # mean density over p in (0.5, 1]
    spike_ratio = float(density[0] / null_density) if null_density > 0 else float("inf")
    frac_lt_05 = float((p < 0.05).mean())
    frac_lt_01 = float((p < 0.01).mean())

    # Interpretation. The healthy shape is a monotone spike at 0 on a ~uniform (density
    # ~= pi0, typically 0.7-0.95) tail. Two pathologies: a tail ABOVE uniform that RISES
    # toward p=1 (a U-shape -> unmodeled correlation/confounding structure), and a tail
    # SUPPRESSED well BELOW uniform (density << pi0-plausible). Variance deflation from
    # DE on removeBatchEffect output (finding 0004) compresses the WHOLE distribution
    # toward 0, which depletes the tail below uniform and inflates the spike -- i.e. the
    # anticonservative tell here is a SUPPRESSED tail, not an elevated one.
    right_end = float(density[-1])
    u_shaped = right_end > 1.20 and right_end > null_density
    if u_shaped:
        verdict = (
            f"U-SHAPED — the tail RISES toward p=1 (right-end density "
            f"{right_end:.2f}x), a sign of unmodeled correlation/confounding "
            "structure; treat with caution"
        )
    elif spike_ratio >= 2.0 and null_density < 0.75:
        verdict = (
            f"SIGNAL PRESENT but TAIL SUPPRESSED — a strong near-zero spike "
            f"({spike_ratio:.1f}x) sits on a null region DEPLETED to "
            f"{null_density:.2f}x uniform (a mostly-null proteome should sit "
            f"~0.7-0.95). A tail this low implies either ~{null_density * 100:.0f}% of "
            "proteins are truly differential (implausibly high for this contrast) OR "
            "the p-values are ANTICONSERVATIVE from variance deflation (finding 0004) "
            "/ the ADD-Kerchner cohort confound inflating the small-p mass. Real "
            "enrichment is present but the SIGNIFICANT-HIT COUNT is likely inflated -- "
            "trust effect sizes + CIs over q-values."
        )
    elif spike_ratio >= 2.0 and null_density <= 1.20:
        verdict = (
            f"WELL-CALIBRATED with real signal — near-uniform null (density "
            f"{null_density:.2f}x) plus a near-zero spike ({spike_ratio:.1f}x the null)"
        )
    else:
        verdict = (
            f"FLAT / no clear enrichment — near-zero spike only {spike_ratio:.1f}x the "
            f"null (density {null_density:.2f}x)"
        )
    return {
        "n_features_tested": n,
        "hist_bin_edges": edges.tolist(),
        "hist_counts": counts.astype(int).tolist(),
        "hist_density_vs_uniform": [round(float(d), 3) for d in density],
        "pi0_storey_upper_half": round(pi0, 3),
        "null_region_density_p_gt_0p5": round(null_density, 3),
        "near_zero_spike_ratio": round(spike_ratio, 3),
        "frac_p_lt_0.05": round(frac_lt_05, 4),
        "frac_p_lt_0.01": round(frac_lt_01, 4),
        "verdict": verdict,
    }


def _parse_accession(feature_id: str) -> str:
    """Return the UniProt accession from a ``db|ACC|NAME_ORGANISM`` id (else the id)."""
    parts = feature_id.split("|")
    return parts[1] if len(parts) >= 3 else feature_id


def annotate_known_markers(contrast_table: pd.DataFrame) -> pd.DataFrame:
    """Add ``gene``/``marker_class`` from the curated marker map (keyed on accession).

    Annotation only — outcome-independent and never used to filter/select features.
    """
    table = contrast_table.copy()
    accessions = table["feature"].astype(str).map(_parse_accession)
    table["gene"] = accessions.map(lambda a: KNOWN_MARKERS.get(a, ("", ""))[0])
    table["marker_class"] = accessions.map(lambda a: KNOWN_MARKERS.get(a, ("", ""))[1])
    return table


def _build_summary(
    result: DifferentialAbundanceResult,
    annotated: pd.DataFrame,
    calibration: dict[str, object],
) -> dict[str, object]:
    """Assemble the JSON summary (hit counts, top hits, effect extremes, markers)."""
    contrast = annotated  # already the contrast table, annotated + q-sorted
    q = contrast["q"].to_numpy(dtype=float)
    n_q05 = int(np.nansum(q < 0.05))
    n_q10 = int(np.nansum(q < 0.10))

    def _direction(effect: float) -> str:
        return "up in ADD" if effect > 0 else "down in ADD"

    top_cols = [
        "feature",
        "gene",
        "marker_class",
        "effect",
        "ci_low",
        "ci_high",
        "p",
        "q",
    ]
    top15 = contrast.head(15).copy()
    top15["direction"] = top15["effect"].map(_direction)
    top_hits = [
        {
            "feature": str(r["feature"]),
            "gene": str(r["gene"]),
            "marker_class": str(r["marker_class"]),
            "log2fc_add_minus_other": round(float(r["effect"]), 3),
            "ci95": [round(float(r["ci_low"]), 3), round(float(r["ci_high"]), 3)],
            "p": float(r["p"]),
            "bh_q": float(r["q"]),
            "direction": str(r["direction"]),
        }
        for _, r in top15[[*top_cols, "direction"]].iterrows()
    ]

    by_effect = contrast.reindex(
        contrast["effect"].abs().sort_values(ascending=False).index
    ).head(10)
    strongest_effects = [
        {
            "feature": str(r["feature"]),
            "gene": str(r["gene"]),
            "log2fc_add_minus_other": round(float(r["effect"]), 3),
            "ci95": [round(float(r["ci_low"]), 3), round(float(r["ci_high"]), 3)],
            "bh_q": float(r["q"]),
            "direction": _direction(float(r["effect"])),
        }
        for _, r in by_effect.iterrows()
    ]

    markers_present = contrast[contrast["marker_class"] != ""]
    marker_hits = [
        {
            "feature": str(r["feature"]),
            "gene": str(r["gene"]),
            "marker_class": str(r["marker_class"]),
            "log2fc_add_minus_other": round(float(r["effect"]), 3),
            "bh_q": float(r["q"]),
            "direction": _direction(float(r["effect"])),
        }
        for _, r in markers_present.iterrows()
    ]

    return {
        "contrast": "ADD vs Other (pooled PDD+PDCN+HCN)",
        "method": "moderated (limma-style empirical-Bayes linear model)",
        "effect": "log2 fold change, ADD minus Other",
        "n_samples": int(result.n_samples),
        "n_add": 10,
        "n_other": 30,
        "n_features_tested": len(result.contrast_table),
        "covariates": list(result.covariates),
        "reference_level": result.reference,
        "prior_variance": (
            None if result.prior_variance is None else float(result.prior_variance)
        ),
        "prior_df": None if result.prior_df is None else float(result.prior_df),
        "n_bh_q_lt_0.05": n_q05,
        "n_bh_q_lt_0.10": n_q10,
        "pvalue_calibration": calibration,
        "top_hits_by_bh_q": top_hits,
        "strongest_effects_by_abs_log2fc": strongest_effects,
        "known_marker_hits": marker_hits,
        "caveats": [
            "NO cohort covariate: ADD is 100% Kerchner (single-cohort), aliased with "
            "the contrast (rank-deficient); the ADD effect is inseparable from a "
            "Kerchner-site effect (findings 0001/0004).",
            "Matrix is already limma removeBatchEffect-corrected on Cohort; the "
            "pre-correction matrix is unavailable, so "
            "batch-as-covariate-on-uncorrected robustness is impossible. DE on "
            "corrected output deflates variance -> "
            "p-values/BH-q are ANTICONSERVATIVE. Weight effect sizes + CIs over "
            "significance.",
            "Exploratory: small n (10 vs 30), imbalanced, confounded.",
        ],
    }


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("data/protein-data.csv"))
    parser.add_argument("--metadata", type=Path, default=Path("data/metadata.csv"))
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--label", type=str, default=DEFAULT_LABEL)
    args = parser.parse_args(argv)

    dataset = load_analysis_subset(load_paired_dataset(args.data, args.metadata))
    logger.info(
        "loaded analysis subset: %d samples x %d features, scale=%s",
        dataset.abundances.shape[0],
        dataset.abundances.shape[1],
        dataset.scale,
    )

    result = run_primary_de(dataset)
    annotated = annotate_known_markers(result.contrast_table)
    calibration = pvalue_calibration_read(result.contrast_table["p"].to_numpy())
    summary = _build_summary(result, annotated, calibration)

    out_dir = args.output_dir / args.label
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "contrast_table.csv"
    # Ship the full per-feature contrast table (every protein), q-sorted, annotated.
    csv_cols = [
        "feature",
        "gene",
        "marker_class",
        "term",
        "effect",
        "ci_low",
        "ci_high",
        "statistic",
        "p",
        "q",
        "mean_abundance",
        "n",
    ]
    annotated[csv_cols].to_csv(csv_path, index=False)
    logger.info("wrote full contrast table (%d rows) -> %s", len(annotated), csv_path)

    summary_path = out_dir / "summary.json"
    with summary_path.open("w") as fh:
        json.dump(summary, fh, indent=2)
    logger.info("wrote summary -> %s", summary_path)

    logger.info(
        "SUMMARY: BH-q<0.05 = %d, BH-q<0.10 = %d | calibration: %s",
        summary["n_bh_q_lt_0.05"],
        summary["n_bh_q_lt_0.10"],
        calibration["verdict"],
    )


if __name__ == "__main__":
    main()
