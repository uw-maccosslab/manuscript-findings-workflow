"""BORUTA ADD-vs-Other — SUMMARY / direction attachment / overlaps (reads the cache).

The *analysis* half of the compute-once / render-many split: it LOADS the cached
``BorutaResult`` (written by ``boruta_add_vs_other.py``) — it never re-runs Boruta — and
derives the finding-facing tables:

  1. Confirmed / Tentative / Rejected counts.
  2. The full Confirmed list with a readable protein short-name and a **DE-direction
     label** attached from the moderated differential-abundance contrast table
     (``results/differential_abundance/add-vs-other-moderated/contrast_table.csv``). The
     sign of ``effect`` = log2FC(ADD - Other): positive => up in ADD, negative => down
     in ADD. A Confirmed feature with ``|log2FC| < --ambiguous-threshold`` (default 0.1)
     is flagged **direction-ambiguous** rather than forced into an up/down bin.
  3. The up / down / ambiguous split of the Confirmed set (feeds direction enrichment).
  4. Overlap of the Confirmed set with (a) the elastic-net minimal-optimal set (the
     non-zero features in the cached EN coefficients.parquet) and (b) the DE hits
     (BH q < --q-threshold). The **redundant-but-real** count = Confirmed features the
     EN classifier zeroed (all-relevant minus minimal-optimal).

Writes a JSON summary to ``results/boruta_add_vs_other/summary.json`` and logs the
headline numbers. Pure post-processing (seconds); no model is fit here.

CAVEATS carried in the JSON: EXPLORATORY; small n (10 positives); the ADD==Kerchner
cohort confound (finding 0001) — "all-relevant to the ADD-vs-Other partition" is aliased
with collection site, so it is not "all-relevant to disease"; the shadow null addresses
relevance-vs-noise, not the confound.

Run::

    ./.venv/bin/python scripts/scratch/boruta_add_vs_other_summary.py \\
        --fingerprint 762e16f294e3
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from boruta_selection import BorutaResult
from common.result_io import load_cached_result

logger = logging.getLogger(__name__)

__script_meta__: dict[str, object] = {
    "task": "boruta-add-vs-other-summary",
    "kind": "analysis",
    "provides": [
        "parse_short_name",
        "attach_direction",
        "confirmed_overlaps",
        "main",
    ],
    "uses": ["boruta_selection", "common.result_io"],
    "seeded_from": {"cache": {"template": "result-io", "version": "0.2"}},
    "description": (
        "Post-processing over a cached Boruta ADD-vs-Other result: Confirmed/Tentative/"
        "Rejected counts, the Confirmed list with DE-direction labels (sign of "
        "log2FC(ADD-Other), |log2FC|<threshold flagged ambiguous), the up/down/"
        "ambiguous split, and overlap of the Confirmed set with the elastic-net "
        "minimal-optimal set and the DE hits (BH q<threshold) -> the "
        "redundant-but-real count. Reads the cache; fits nothing. EXPLORATORY; "
        "ADD==Kerchner confound."
    ),
}

DEFAULT_CACHE_ROOT = Path("results")
DEFAULT_DE_TABLE = Path(
    "results/differential_abundance/add-vs-other-moderated/contrast_table.csv"
)
# The cached elastic-net NULL-pass result whose non-zero coefficients are the
# minimal-optimal 257 (results/manifest.md; task brief).
DEFAULT_EN_COEF = Path("results/classification/902c1e98c018/coefficients.parquet")
DEFAULT_OUT_DIR = Path("results/boruta_add_vs_other")


def parse_short_name(feature_id: str) -> str:
    """Return the protein short-name (``NAME`` of ``db|ACC|NAME_ORGANISM``), else id.

    The DE table's own ``gene`` column is essentially empty for this dataset, so the
    readable identifier is the UniProt entry mnemonic (e.g. ``MP2K2`` from
    ``sp|P36507|MP2K2_HUMAN``). This is a mnemonic, not a curated HGNC symbol.
    """
    parts = feature_id.split("|")
    if len(parts) >= 3:
        name = parts[2]
        return name.rsplit("_", 1)[0] if "_" in name else name
    return feature_id


def attach_direction(
    confirmed: np.ndarray,
    de_table: pd.DataFrame,
    *,
    ambiguous_threshold: float,
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    """Attach the DE log2FC(ADD-Other) direction to each Confirmed feature.

    Returns the per-feature rows (importance-descending order from ``confirmed``) and
    the up/down/ambiguous/not-in-DE split counts.
    """
    de_by_feature = de_table.set_index("feature")
    rows: list[dict[str, Any]] = []
    split = {"up_in_add": 0, "down_in_add": 0, "ambiguous": 0, "not_in_de_table": 0}
    for feature in confirmed:
        feature = str(feature)
        in_de = feature in de_by_feature.index
        effect = float(de_by_feature.loc[feature, "effect"]) if in_de else float("nan")
        q = float(de_by_feature.loc[feature, "q"]) if in_de else float("nan")
        marker = de_by_feature.loc[feature, "marker_class"] if in_de else np.nan
        marker_str = None if (isinstance(marker, float) and pd.isna(marker)) else marker
        if not in_de:
            direction = "not_in_de_table"
            split["not_in_de_table"] += 1
        elif not np.isfinite(effect) or abs(effect) < ambiguous_threshold:
            direction = "ambiguous"
            split["ambiguous"] += 1
        elif effect > 0:
            direction = "up_in_add"
            split["up_in_add"] += 1
        else:
            direction = "down_in_add"
            split["down_in_add"] += 1
        rows.append(
            {
                "feature": feature,
                "short_name": parse_short_name(feature),
                "de_log2fc_add_minus_other": (
                    None if not np.isfinite(effect) else round(effect, 4)
                ),
                "de_bh_q": None if not np.isfinite(q) else q,
                "de_bh_q_lt_thresh": bool(np.isfinite(q) and q < 0.05),
                "direction": direction,
                "marker_class": None if marker_str is None else str(marker_str),
            }
        )
    return rows, split


def confirmed_overlaps(
    confirmed: set[str],
    en_set: set[str],
    de_hits: set[str],
) -> dict[str, Any]:
    """Overlap of the Confirmed set with the EN minimal-optimal set and the DE hits."""
    in_en = confirmed & en_set
    not_in_en = confirmed - en_set  # redundant-but-real: EN zeroed, Boruta confirmed
    in_de = confirmed & de_hits
    not_in_de = confirmed - de_hits
    return {
        "n_confirmed": len(confirmed),
        "n_en_minimal_optimal": len(en_set),
        "n_de_hits": len(de_hits),
        "confirmed_in_en": len(in_en),
        "confirmed_not_in_en_redundant_but_real": len(not_in_en),
        "confirmed_in_de_hits": len(in_de),
        "confirmed_not_in_de_hits": len(not_in_de),
        "redundant_but_real_features": sorted(
            {parse_short_name(f) for f in not_in_en}
        ),
    }


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fingerprint", type=str, required=True)
    parser.add_argument("--cache-root", type=Path, default=DEFAULT_CACHE_ROOT)
    parser.add_argument("--de-table", type=Path, default=DEFAULT_DE_TABLE)
    parser.add_argument("--en-coef", type=Path, default=DEFAULT_EN_COEF)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--ambiguous-threshold", type=float, default=0.1)
    parser.add_argument("--q-threshold", type=float, default=0.05)
    args = parser.parse_args(argv)

    result = load_cached_result(
        BorutaResult,
        cache_root=args.cache_root,
        analysis="boruta",
        fingerprint=args.fingerprint,
    )
    logger.info(
        "loaded boruta result %s: %d Confirmed / %d Tentative / %d Rejected",
        args.fingerprint,
        result.n_confirmed,
        result.n_tentative,
        result.n_rejected,
    )

    de_table = pd.read_csv(args.de_table)
    de_hits = set(de_table.loc[de_table["q"] < args.q_threshold, "feature"].tolist())
    en_set = set(pd.read_parquet(args.en_coef)["feature"].tolist())

    confirmed = result.confirmed_features  # importance-descending
    rows, split = attach_direction(
        confirmed, de_table, ambiguous_threshold=args.ambiguous_threshold
    )
    overlaps = confirmed_overlaps(
        {str(f) for f in confirmed}, en_set, de_hits
    )

    tentative = sorted(
        parse_short_name(str(f))
        for f in result.feature_names[result.decision == "Tentative"]
    )

    summary: dict[str, Any] = {
        "analysis": "BORUTA all-relevant selection — ADD vs Other (summary)",
        "boruta_fingerprint": args.fingerprint,
        "boruta_params": {
            "max_iter": result.max_iter,
            "alpha": result.alpha,
            "perc": result.perc,
            "seed": result.random_state,
            "n_iter_run": result.n_iter,
            "task": result.task,
            "classes": list(result.classes) if result.classes else None,
            "n_samples": result.n_samples,
            "n_features": result.n_features,
            "shadow_threshold": round(result.shadow_threshold, 5),
        },
        "counts": {
            "confirmed": result.n_confirmed,
            "tentative": result.n_tentative,
            "rejected": result.n_rejected,
        },
        "direction_split_of_confirmed": split,
        "ambiguous_threshold_abs_log2fc": args.ambiguous_threshold,
        "overlaps": overlaps,
        "confirmed_features": rows,
        "tentative_features_short_names": tentative,
        "caveats": [
            "EXPLORATORY. Small n (10 positives / 30 negatives); Boruta's shadow null "
            "addresses relevance-vs-noise, not sampling stability at this n.",
            "COHORT CONFOUND (hard, finding 0001): ADD is 100% Kerchner cohort. Among "
            "these 40 samples 'ADD' and 'Kerchner' are the SAME partition, so the "
            "Confirmed set is all-relevant to the ADD-vs-Other partition, which is "
            "ALIASED with collection site — NOT necessarily all-relevant to disease. "
            "The shadow null does NOT resolve this confound.",
            "Matrix is already median-normalized + limma/ComBat batch-corrected on "
            "Cohort (finding 0004); feeding Boruta the corrected matrix is legitimate "
            "for selection (not a per-feature significance test).",
            "Direction labels are taken from the moderated DE contrast (log2FC "
            "ADD-Other); the DE 'gene' column is empty for this dataset, so short_name "
            "is the UniProt entry mnemonic, not a curated HGNC symbol.",
        ],
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.out_dir / "summary.json"
    with out_path.open("w") as fh:
        json.dump(summary, fh, indent=2)
    logger.info("wrote summary -> %s", out_path)
    logger.info(
        "direction split: %d up / %d down / %d ambiguous / %d not-in-DE",
        split["up_in_add"],
        split["down_in_add"],
        split["ambiguous"],
        split["not_in_de_table"],
    )
    logger.info(
        "overlaps: %d/%d Confirmed in EN-257; %d redundant-but-real (EN-zeroed); "
        "%d/%d Confirmed in DE hits (q<%.2f)",
        overlaps["confirmed_in_en"],
        overlaps["n_confirmed"],
        overlaps["confirmed_not_in_en_redundant_but_real"],
        overlaps["confirmed_in_de_hits"],
        overlaps["n_confirmed"],
        args.q_threshold,
    )


if __name__ == "__main__":
    main()
