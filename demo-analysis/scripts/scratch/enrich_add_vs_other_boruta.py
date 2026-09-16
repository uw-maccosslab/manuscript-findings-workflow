"""GO/KEGG over-representation of the Boruta Confirmed set: ADD vs Other — COMPUTE.

Call site for the project copy of the vetted ``lib/analysis/enrichment`` template
(``scripts/scratch/enrichment.py``, ``enrich``). Runs **directional** GO/KEGG
over-representation (ORA) via g:Profiler on the Boruta *all-relevant* Confirmed set for
the ADD-vs-Other contrast, **separately for the up-in-ADD and down-in-ADD** proteins,
each against the experiment's **detected-proteome background** (the 2340-feature Stage-4
analysis subset — ``domain_scope="custom"``, NOT the whole genome).

What it does:

  1. Loads the 40-feature **Confirmed** set by round-tripping the cached
     ``BorutaResult`` (``results/boruta/762e16f294e3/`` via
     ``result_io.load_cached_result``).
  2. Assigns each Confirmed feature a **direction** from the moderated DE contrast
     (``results/differential_abundance/add-vs-other-moderated/contrast_table.csv``):
     ``effect`` is log2FC ADD-minus-Other, so ``effect > 0`` -> up-in-ADD, ``< 0`` ->
     down-in-ADD. A |log2FC| below ``AMBIGUOUS_EPS`` is flagged **ambiguous** and
     dropped from both queries (expected: 0). The split is asserted 18 up / 22 down /
     0 ambiguous (a data-change guard).
  3. Decodes UniProt accessions from the feature ids with the template's signed-off
     ``parse_uniprot_accession`` (``sp|ACC|NAME`` -> ``ACC``; a `` / `` protein group ->
     the first member's accession). Loads the **background** = the 2340 analysis-subset
     feature accessions via ``project_loader.load_analysis_subset``.
  4. Records the decode audit (``build_mapping_audit``: per-id ``original ->
     accession``, source, in-query/in-background) plus per-source counts and mappable %.
  5. Runs ``enrich`` for the UP query and the DOWN query separately (organism
     ``hsapiens``; sources GO:BP/MF/CC/KEGG; correction **g:SCS** — the g:Profiler
     default the scientist chose; ``all_results=True`` so every tested term is kept).
  6. Writes per-term CSVs (full + significant-only), the decode audit CSV, and a
     provenance JSON that pins the **g:Profiler data version** to
     ``results/enrichment/``.

CAVEATS (state honestly; findings 0001, 0007; state/DATA_DESCRIPTION.md):

  * **EXPLORATORY / limited power.** The queries are small (18 up, 22 down); at this
    size few terms may survive g:SCS. Report "none significant" honestly if so.
  * **CONFOUND CEILING (hard, finding 0001).** The Confirmed set is all-relevant to
    the ADD-vs-Other partition, which is ALIASED with collection site (ADD is 100%
    Kerchner cohort). This ORA therefore characterizes *what compartment/pathway the
    discriminative proteins belong to* — it does NOT establish them as disease biology
    (could be site/handling). Enrichment describes the set; it does not resolve
    disease-vs-site.
  * **Network + reproducibility.** g:Profiler is a live, versioned service; the returned
    data version is captured to provenance. A network/service failure RAISES
    (``EnrichmentServiceError``) — never a silent empty result.

Non-interactive, parameterized (argparse), seeded-from the enrichment template.

Run::

    ./.venv/bin/python scripts/scratch/enrich_add_vs_other_boruta.py
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

# Make `scripts/` and `scripts/scratch/` importable regardless of invocation, mirroring
# boruta_add_vs_other.py / de_add_vs_other.py.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from boruta_selection import BorutaResult
from common.result_io import load_cached_result
from enrichment import (
    DEFAULT_GPROFILER_URL,
    DEFAULT_SOURCES,
    Correction,
    EnrichmentResult,
    build_mapping_audit,
    enrich,
    parse_uniprot_accession,
)
from project_loader import load_analysis_subset, load_paired_dataset

logger = logging.getLogger(__name__)

__script_meta__: dict[str, object] = {
    "task": "enrich-add-vs-other-boruta-confirmed",
    "kind": "analysis",
    "provides": [
        "ORGANISM",
        "CORRECTION",
        "DE_TERM",
        "BORUTA_FINGERPRINT",
        "DATA_VERSION",
        "GProfilerQuery",
        "confirmed_directions",
        "run_directional_enrichment",
        "main",
    ],
    "uses": [
        "enrichment",
        "boruta_selection",
        "project_loader",
        "common.result_io",
    ],
    "seeded_from": {"enrichment": {"template": "enrichment", "version": "0.1"}},
    "description": (
        "Directional GO(BP/MF/CC)+KEGG over-representation (g:Profiler, g:SCS) of the "
        "Boruta all-relevant Confirmed set (40) for ADD-vs-Other, run SEPARATELY for "
        "the 18 up-in-ADD and 22 down-in-ADD proteins against the 2340-protein "
        "detected-proteome background (domain_scope=custom). Directions taken from the "
        "moderated DE contrast (log2FC ADD-Other). UniProt accessions decoded with the "
        "template's signed-off parser (protein group -> first member). Writes per-term "
        "CSVs (full + significant), the decode audit, and a provenance JSON pinning "
        "the g:Profiler data version. EXPLORATORY (small queries); ADD==Kerchner site "
        "confound (finding 0001) caps this at describing the set, not disease biology."
    ),
}

ORGANISM = "hsapiens"
CORRECTION: Correction = "g_SCS"
DE_TERM = "ADDvsOther[ADD vs Other]"
BORUTA_FINGERPRINT = "762e16f294e3"
DATA_VERSION = "protein-data.csv sha256:90a0462eb03a02bd"

# |log2FC| below this is treated as directionally ambiguous and dropped from both
# queries (matches the summary's ambiguous_threshold_abs_log2fc = 0.1). Expect 0.
AMBIGUOUS_EPS = 0.1

EXPECTED_N_CONFIRMED = 40
EXPECTED_N_UP = 18
EXPECTED_N_DOWN = 22
EXPECTED_N_BACKGROUND = 2340

DEFAULT_CACHE_ROOT = Path("results")
DEFAULT_OUT_DIR = Path("results/enrichment/add-vs-other-boruta-confirmed")

_TERM_COLUMNS = [
    "source",
    "term_id",
    "term_name",
    "p_adj_gSCS",
    "significant",
    "intersection_size",
    "term_size",
    "query_size",
    "gene_ratio",
    "recall",
    "query_genes",
]


@dataclass(frozen=True)
class GProfilerQuery:
    """One directional query: the raw feature ids and their decoded accessions."""

    label: str
    feature_ids: tuple[str, ...]
    accessions: tuple[str, ...]


def _load_confirmed(cache_root: Path, fingerprint: str) -> list[str]:
    """Round-trip the cached BorutaResult and return the Confirmed feature ids."""
    result = load_cached_result(
        BorutaResult, cache_root=cache_root, analysis="boruta", fingerprint=fingerprint
    )
    confirmed = [str(f) for f in result.confirmed_features]
    if len(confirmed) != EXPECTED_N_CONFIRMED:
        raise ValueError(
            f"Boruta cache {fingerprint} has {len(confirmed)} Confirmed features; "
            f"expected {EXPECTED_N_CONFIRMED}."
        )
    logger.info(
        "loaded %d Confirmed features from cache %s", len(confirmed), fingerprint
    )
    return confirmed


def _load_de_effects(de_table: Path) -> dict[str, float]:
    """Map feature id -> log2FC (ADD - Other) from the moderated DE contrast table."""
    table = pd.read_csv(de_table)
    for col in ("feature", "effect", "term"):
        if col not in table.columns:
            raise ValueError(f"DE table {de_table} is missing required column {col!r}.")
    terms = set(table["term"].astype(str).unique())
    if terms != {DE_TERM}:
        raise ValueError(
            f"DE table {de_table} terms {sorted(terms)} != expected {{{DE_TERM!r}}}."
        )
    effects = {
        str(f): float(e)
        for f, e in zip(table["feature"], table["effect"], strict=True)
    }
    return effects


def confirmed_directions(
    confirmed: list[str], effects: dict[str, float]
) -> tuple[GProfilerQuery, GProfilerQuery, list[str]]:
    """Split the Confirmed set into up/down queries by the sign of log2FC (ADD-Other).

    Returns ``(up_query, down_query, ambiguous_feature_ids)``. A feature whose |log2FC|
    is below :data:`AMBIGUOUS_EPS`, or that is absent from the DE table, is ambiguous
    and excluded from both queries. Fails loud unless the split is exactly 18 up /
    22 down / 0 ambiguous (a data-change guard aligned with the Boruta summary).
    """
    up_ids: list[str] = []
    down_ids: list[str] = []
    ambiguous: list[str] = []
    for feature in confirmed:
        if feature not in effects:
            ambiguous.append(feature)
            continue
        log2fc = effects[feature]
        if abs(log2fc) < AMBIGUOUS_EPS:
            ambiguous.append(feature)
        elif log2fc > 0:
            up_ids.append(feature)
        else:
            down_ids.append(feature)
    if (len(up_ids), len(down_ids), len(ambiguous)) != (
        EXPECTED_N_UP,
        EXPECTED_N_DOWN,
        0,
    ):
        raise ValueError(
            f"Direction split is {len(up_ids)} up / {len(down_ids)} down / "
            f"{len(ambiguous)} ambiguous; want {EXPECTED_N_UP}/{EXPECTED_N_DOWN}/0. "
            f"Ambiguous: {ambiguous}."
        )
    up = GProfilerQuery(
        label="up in ADD",
        feature_ids=tuple(up_ids),
        accessions=tuple(_decode(up_ids)),
    )
    down = GProfilerQuery(
        label="down in ADD",
        feature_ids=tuple(down_ids),
        accessions=tuple(_decode(down_ids)),
    )
    logger.info(
        "direction split: %d up / %d down / %d ambiguous",
        len(up_ids),
        len(down_ids),
        len(ambiguous),
    )
    return up, down, ambiguous


def _decode(feature_ids: list[str]) -> list[str]:
    """Decode UniProt accessions; fail loud on any miss (decode is verified clean)."""
    accessions: list[str] = []
    misses: list[str] = []
    for fid in feature_ids:
        acc = parse_uniprot_accession(fid)
        if acc is None:
            misses.append(fid)
        else:
            accessions.append(acc)
    if misses:
        raise ValueError(
            f"{len(misses)} feature id(s) failed to decode to a UniProt accession "
            f"(decode was verified 100% clean): {misses}."
        )
    return accessions


def _term_table(result: EnrichmentResult) -> pd.DataFrame:
    """Project the enrichment table to the recorded output columns (all terms)."""
    table = result.table.rename(
        columns={"p_value": "p_adj_gSCS", "intersecting_genes": "query_genes"}
    )
    return table[_TERM_COLUMNS].copy()


def _source_counts(table: pd.DataFrame) -> dict[str, dict[str, int]]:
    """Per-source tested-term and significant-term counts for the summary."""
    counts: dict[str, dict[str, int]] = {}
    for source, group in table.groupby("source"):
        counts[str(source)] = {
            "tested": len(group),
            "significant": int(group["significant"].sum()),
        }
    return counts


def _mapping_summary(mapping: pd.DataFrame) -> dict[str, object]:
    """Per-source decode counts + mappable % + any unmappable ids for the audit."""
    by_source = {
        str(s): int(n) for s, n in mapping["source"].value_counts().to_dict().items()
    }
    n_total = len(mapping)
    n_mappable = int((mapping["uniprot_accession"] != "").sum())
    unmappable = mapping.loc[
        mapping["uniprot_accession"] == "", "original_id"
    ].tolist()
    mappable = mapping.loc[
        mapping["uniprot_accession"] != "", ["original_id", "uniprot_accession"]
    ]
    examples = mappable.head(5).to_dict("records")
    return {
        "n_distinct_ids": n_total,
        "n_mappable": n_mappable,
        "mappable_pct": round(100.0 * n_mappable / n_total, 4) if n_total else 0.0,
        "by_source": by_source,
        "n_unmappable": len(unmappable),
        "unmappable_ids": unmappable,
        "examples": examples,
    }


def run_directional_enrichment(
    query: GProfilerQuery,
    background: list[str],
    *,
    sources: tuple[str, ...],
    correction: Correction,
    base_url: str,
    out_dir: Path,
) -> tuple[EnrichmentResult, dict[str, object]]:
    """Run g:Profiler ORA for one directional query and write its per-term CSVs."""
    result = enrich(
        list(query.feature_ids),
        background,
        ORGANISM,
        sources=list(sources),
        correction=correction,
        query_label=query.label,
        base_url=base_url,
    )
    table = _term_table(result)
    slug = query.label.replace(" ", "_")
    full_path = out_dir / f"{slug}_terms.csv"
    sig_path = out_dir / f"{slug}_significant.csv"
    table.to_csv(full_path, index=False)
    table[table["significant"]].to_csv(sig_path, index=False)
    n_sig = int(table["significant"].sum())
    logger.info(
        "%s: %d accessions, %d recognized; %d/%d terms significant (g:SCS<=%.3f); "
        "gprofiler_version=%s",
        query.label,
        result.n_query_mapped,
        result.n_query_recognized,
        n_sig,
        len(table),
        result.user_threshold,
        result.gprofiler_version,
    )
    summary: dict[str, object] = {
        "label": query.label,
        "n_query_features": len(query.feature_ids),
        "n_query_mapped": result.n_query_mapped,
        "n_query_recognized": result.n_query_recognized,
        "query_accessions": list(query.accessions),
        "n_terms_tested": len(table),
        "n_terms_significant": n_sig,
        "per_source": _source_counts(table),
        "gprofiler_version": result.gprofiler_version,
        "effective_domain_size": result.effective_domain_size,
        "csv_full": str(full_path),
        "csv_significant": str(sig_path),
    }
    return result, summary


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("data/protein-data.csv"))
    parser.add_argument("--metadata", type=Path, default=Path("data/metadata.csv"))
    parser.add_argument("--cache-root", type=Path, default=DEFAULT_CACHE_ROOT)
    parser.add_argument("--boruta-fingerprint", type=str, default=BORUTA_FINGERPRINT)
    parser.add_argument(
        "--de-table",
        type=Path,
        default=Path(
            "results/differential_abundance/add-vs-other-moderated/contrast_table.csv"
        ),
    )
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument(
        "--base-url",
        type=str,
        default=None,
        help="Pin g:Profiler to an archived release for exact reproducibility.",
    )
    args = parser.parse_args(argv)

    base_url: str = (
        args.base_url if args.base_url is not None else DEFAULT_GPROFILER_URL
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)

    confirmed = _load_confirmed(args.cache_root, args.boruta_fingerprint)
    effects = _load_de_effects(args.de_table)
    up, down, ambiguous = confirmed_directions(confirmed, effects)

    dataset = load_analysis_subset(load_paired_dataset(args.data, args.metadata))
    background = [str(x) for x in dataset.feature_names]
    if len(background) != EXPECTED_N_BACKGROUND:
        raise ValueError(
            f"Background has {len(background)} features; expected "
            f"{EXPECTED_N_BACKGROUND} (the analysis subset)."
        )

    # One decode audit over the whole Confirmed set + background (queries are subsets).
    mapping = build_mapping_audit(confirmed, background)
    audit_path = args.out_dir / "decode_audit.csv"
    mapping.to_csv(audit_path, index=False)
    mapping_summary = _mapping_summary(mapping)
    logger.info(
        "decode audit: %d/%d distinct ids mappable (%.2f%%); %d unmappable",
        mapping_summary["n_mappable"],
        mapping_summary["n_distinct_ids"],
        mapping_summary["mappable_pct"],
        mapping_summary["n_unmappable"],
    )

    sources = tuple(str(s) for s in DEFAULT_SOURCES)
    up_result, up_summary = run_directional_enrichment(
        up,
        background,
        sources=sources,
        correction=CORRECTION,
        base_url=base_url,
        out_dir=args.out_dir,
    )
    down_result, down_summary = run_directional_enrichment(
        down,
        background,
        sources=sources,
        correction=CORRECTION,
        base_url=base_url,
        out_dir=args.out_dir,
    )

    versions = {up_result.gprofiler_version, down_result.gprofiler_version}
    if len(versions) != 1:
        logger.warning("g:Profiler versions differ between queries: %s", versions)

    provenance: dict[str, object] = {
        "analysis": "GO/KEGG over-representation (g:Profiler) of the Boruta Confirmed "
        "set — ADD vs Other, directional (up/down separately)",
        "seeded_from": {"template": "enrichment", "version": "0.1"},
        "organism": ORGANISM,
        "sources": list(sources),
        "correction": CORRECTION,
        "user_threshold": up_result.user_threshold,
        "domain_scope": "custom (detected-proteome background)",
        "gprofiler_version": sorted(versions),
        "base_url": base_url,
        "base_url_pinned": args.base_url is not None,
        "data_version": DATA_VERSION,
        "boruta_fingerprint": args.boruta_fingerprint,
        "de_table": str(args.de_table),
        "de_effect_orientation": "log2FC ADD minus Other",
        "ambiguous_eps_abs_log2fc": AMBIGUOUS_EPS,
        "n_confirmed": len(confirmed),
        "n_up": len(up.feature_ids),
        "n_down": len(down.feature_ids),
        "n_ambiguous": len(ambiguous),
        "n_background": len(background),
        "decode_audit": mapping_summary,
        "decode_audit_csv": str(audit_path),
        "up": up_summary,
        "down": down_summary,
        "caveats": [
            "EXPLORATORY: small queries (18 up, 22 down) -> limited power; few terms "
            "may survive g:SCS.",
            "CONFOUND CEILING (finding 0001): the Confirmed set is all-relevant to the "
            "ADD-vs-Other partition, which is ALIASED with collection site (ADD is "
            "100% Kerchner cohort). This ORA characterizes WHAT compartment/pathway "
            "discriminative proteins belong to; it does NOT establish them as disease "
            "biology (could be site/handling). Enrichment describes the set; it does "
            "not resolve disease-vs-site.",
            "g:Profiler is a live versioned service; the data version is pinned above. "
            "A network/service failure raises EnrichmentServiceError (never empty).",
        ],
    }
    prov_path = args.out_dir / "enrichment_provenance.json"
    prov_path.write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    logger.info("wrote provenance -> %s", prov_path)
    logger.info(
        "DONE: up %d sig / down %d sig terms; gprofiler_version=%s",
        up_summary["n_terms_significant"],
        down_summary["n_terms_significant"],
        sorted(versions),
    )


if __name__ == "__main__":
    main()
