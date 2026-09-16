"""Preranked GSEA of the ADD-vs-Other moderated contrast — COMPUTE.

Rank-based (abundance-weighted) enrichment follow-up to the count-based ORA of
finding 0009, which found NOTHING significant in GO:BP/CC/MF or KEGG on the 40-protein
discriminative hit-set. ORA is unweighted and thresholded; **preranked GSEA** tests
whether a gene set's members are *coordinately* shifted across the FULL effect-ranked
proteome, weighting every protein by the magnitude of its moderated t. It can surface
coordinated signal where a thresholded ORA is null. That is the hypothesis under test.

Pipeline:

  1. Load the pre-computed moderated DE contrast table (2340 proteins, one per row;
     `statistic` = moderated t = ADD - Other, signed). DE is NOT recomputed here.
  2. Decode a UniProt accession from each `feature` (`sp|ACC|NAME_HUMAN`; protein group
     `A / B` -> FIRST/representative member per the Stage-2 protein-group policy) using
     the project's signed-off `enrichment.parse_uniprot_accession`.
  3. Map accession -> HGNC symbol via mygene (`scopes=uniprot, fields=symbol,
     species=human`) and VERIFY: the UniProt entry name (`MP2K2`) is NOT the HGNC symbol
     (`MAP2K2`), so the entry name is never used. Report a decode audit (mapped/unmapped
     of 2340). Unmapped are dropped. Multiple accessions collapsing to one symbol: keep
     the entry with max |t|.
  4. Rank ALL mapped, deduped symbols by the signed moderated t, descending.
  5. gseapy `prerank` over GO_Biological_Process_2021, GO_Molecular_Function_2021,
     GO_Cellular_Component_2021, KEGG_2021_Human (Enrichr, symbol-keyed). Fixed seed.
     GSEA's gene-set-permutation null is intrinsic (no separate label-shuffle needed).
  6. Cache result tables + ranked list + mapping audit + meta.json under
     results/gsea/<fingerprint>/ (fingerprint over DE-table sha / data-version / ranking
     metric / libraries / seed / gseapy version / mapping method).

CAVEATS (carry into every interpretation):

  * CONFOUND (findings 0001, 0004): ADD is 100% Kerchner cohort; the matrix is limma
    batch-corrected on cohort with EXTRAPOLATED correction for the single-cohort
    ADD/PDCN groups. Any GSEA hit is inseparable from an ADD/Kerchner collection-SITE
    effect -- confound-bounded, NOT disease-proven. The ranking metric derives from a
    confounded contrast.
  * ANTICONSERVATIVE INPUT (finding 0006): the DE p/q are anticonservative (variance
    deflated by upstream batch correction; pi0 ~ 0.53). The t-statistics inherit this;
    the RELATIVE ranking stays meaningful for GSEA, but absolute significance is
    inflated.
  * SMALL n / EXPLORATORY (finding 0003): n=10 ADD vs 30 Other. Everything is
    exploratory.

Deterministic; seeded. Not promoted (scratch).

Run::

    ./.venv/bin/python scripts/scratch/gsea_add_vs_other.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import sys
from pathlib import Path

import gseapy as gp
import mygene
import pandas as pd

# Make `scripts/` and `scripts/scratch/` importable regardless of invocation.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from common.result_io import result_fingerprint  # noqa: E402
from enrichment import parse_uniprot_accession  # noqa: E402

logger = logging.getLogger(__name__)

__script_meta__: dict[str, object] = {
    "task": "gsea-preranked-add-vs-other",
    "kind": "analysis",
    "provides": ["decode_symbols", "build_ranked_list", "run_prerank", "main"],
    "uses": ["enrichment", "common.result_io", "gseapy", "mygene"],
    "seeded_from": {"enrichment": {"template": "enrichment", "version": "0.1"}},
    "description": (
        "Preranked (abundance-weighted) GSEA over GO:BP/MF/CC + KEGG of the ADD-vs-Other "
        "moderated contrast, ranked by signed moderated t. UniProt->HGNC via mygene, "
        "verified; protein group -> first member; collisions keep max|t|. GSEA's "
        "gene-set-permutation null is intrinsic. Follow-up to null ORA finding 0009. "
        "Confound-bounded (ADD==Kerchner site), anticonservative input, small-n -> "
        "exploratory."
    ),
}

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DE = (
    REPO
    / "results/differential_abundance/add-vs-other-moderated/contrast_table.csv"
)
DEFAULT_CACHE_ROOT = REPO / "results"
DATA_VERSION = "protein-data.csv sha256:90a0462eb03a02bd"
DE_TERM = "ADDvsOther[ADD vs Other]"
RANKING_METRIC = "moderated_t_signed_desc"
MAPPING_METHOD = "mygene.querymany(scopes=uniprot,fields=symbol,species=human)"
LIBRARIES = [
    "GO_Biological_Process_2021",
    "GO_Molecular_Function_2021",
    "GO_Cellular_Component_2021",
    "KEGG_2021_Human",
]
SEED = 20260915
MIN_SIZE = 15
MAX_SIZE = 500
PERMUTATION_NUM = 1000


def sha16(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16]


def decode_symbols(de: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Attach a HGNC symbol to each DE row; return (ranked-input frame, audit).

    Steps: parse first accession (protein group -> representative) -> mygene map ->
    drop unmapped -> collapse duplicate symbols keeping the max-|t| entry.
    """
    df = de.copy()
    df["accession"] = df["feature"].map(parse_uniprot_accession)
    n_total = len(df)
    n_no_acc = int(df["accession"].isna().sum())

    accessions = sorted(df["accession"].dropna().unique().tolist())
    mg = mygene.MyGeneInfo()
    hits = mg.querymany(
        accessions,
        scopes="uniprot",
        fields="symbol",
        species="human",
        verbose=False,
    )
    # Keep the best-scoring symbol per accession (querymany may return multiple hits).
    acc2sym: dict[str, str] = {}
    acc2score: dict[str, float] = {}
    for h in hits:
        q = h.get("query")
        sym = h.get("symbol")
        if sym is None or h.get("notfound"):
            continue
        score = float(h.get("_score", 0.0) or 0.0)
        if q not in acc2score or score > acc2score[q]:
            acc2sym[q] = sym
            acc2score[q] = score

    df["symbol"] = df["accession"].map(acc2sym)
    n_unmapped = int(df["symbol"].isna().sum())
    n_mapped_rows = n_total - n_unmapped

    mapped = df.dropna(subset=["symbol"]).copy()
    mapped["abs_t"] = mapped["statistic"].abs()
    n_symbols_before = mapped["symbol"].nunique()

    # Collision handling: multiple accessions -> same symbol -> keep max |t|.
    mapped = mapped.sort_values("abs_t", ascending=False)
    n_collisions = len(mapped) - n_symbols_before
    deduped = mapped.drop_duplicates(subset="symbol", keep="first")

    ranked = deduped[["symbol", "accession", "feature", "statistic", "effect"]].copy()
    ranked = ranked.sort_values("statistic", ascending=False).reset_index(drop=True)

    audit = {
        "n_total_de_rows": n_total,
        "n_no_accession_parsed": n_no_acc,
        "n_accessions_queried": len(accessions),
        "n_accessions_mapped_to_symbol": len(acc2sym),
        "n_rows_mapped": n_mapped_rows,
        "n_rows_unmapped_dropped": n_unmapped,
        "n_unique_symbols_before_dedup": int(n_symbols_before),
        "n_symbol_collisions_dropped_keep_max_abs_t": int(n_collisions),
        "n_symbols_in_ranked_list": int(len(ranked)),
    }
    return ranked, audit


def run_prerank(ranked: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Run gseapy prerank per library; return {library: res2d frame}."""
    rnk = ranked[["symbol", "statistic"]].copy()
    results: dict[str, pd.DataFrame] = {}
    for lib in LIBRARIES:
        logger.info("prerank over %s", lib)
        pre = gp.prerank(
            rnk=rnk,
            gene_sets=lib,
            min_size=MIN_SIZE,
            max_size=MAX_SIZE,
            permutation_num=PERMUTATION_NUM,
            seed=SEED,
            threads=4,
            outdir=None,
            verbose=False,
        )
        res = pre.res2d.copy()
        res["Library"] = lib
        res["leading_edge_size"] = res["Lead_genes"].map(
            lambda s: len([g for g in str(s).split(";") if g])
        )
        results[lib] = res
    return results


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--de", type=Path, default=DEFAULT_DE)
    ap.add_argument("--cache-root", type=Path, default=DEFAULT_CACHE_ROOT)
    args = ap.parse_args()

    de = pd.read_csv(args.de)
    if len(de) != 2340:
        raise ValueError(f"expected 2340 DE rows, got {len(de)} (data changed?)")
    if de["statistic"].isna().any():
        raise ValueError("NaN in ranking statistic column")

    de_sha = sha16(args.de)
    ranked, audit = decode_symbols(de)
    logger.info("decode audit: %s", json.dumps(audit, indent=2))

    results = run_prerank(ranked)

    params = {
        "contrast": DE_TERM,
        "de_table_sha256_16": de_sha,
        "ranking_metric": RANKING_METRIC,
        "libraries": LIBRARIES,
        "min_size": MIN_SIZE,
        "max_size": MAX_SIZE,
        "permutation_num": PERMUTATION_NUM,
        "gseapy_version": gp.__version__,
        "mapping_method": MAPPING_METHOD,
        "sample_set": "experimental (n=40; 10 ADD vs 30 Other)",
    }
    fp = result_fingerprint(
        analysis="gsea", data_version=DATA_VERSION, params=params, seed=SEED
    )
    outdir = args.cache_root / "gsea" / fp
    outdir.mkdir(parents=True, exist_ok=True)

    ranked.to_csv(outdir / "ranked_symbols.csv", index=False)
    pd.DataFrame([audit]).to_csv(outdir / "decode_audit.csv", index=False)
    all_res = []
    for lib, res in results.items():
        short = lib.replace("_2021", "").replace("_Human", "").lower()
        res.to_csv(outdir / f"gsea_{short}.csv", index=False)
        all_res.append(res)
    combined = pd.concat(all_res, ignore_index=True)
    combined.to_csv(outdir / "gsea_all_libraries.csv", index=False)

    meta = {
        "analysis": "gsea",
        "data_version": DATA_VERSION,
        "params": params,
        "seed": SEED,
        "fingerprint": fp,
        "label": (
            "ADD-vs-Other preranked GSEA (signed moderated t; GO:BP/MF/CC + KEGG 2021; "
            "mygene UniProt->symbol). Follow-up to null ORA finding 0009."
        ),
        "created": None,
        "decode_audit": audit,
    }
    (outdir / "meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    # Console summary.
    def q(x: pd.DataFrame) -> pd.Series:
        return pd.to_numeric(x["FDR q-val"], errors="coerce")

    def nes(x: pd.DataFrame) -> pd.Series:
        return pd.to_numeric(x["NES"], errors="coerce")

    print("\n=== GSEA summary ===")
    print(f"cache: {outdir}")
    for lib, res in results.items():
        qv, nv = q(res), nes(res)
        for thr in (0.25, 0.05):
            pos = int(((qv < thr) & (nv > 0)).sum())
            neg = int(((qv < thr) & (nv < 0)).sum())
            print(f"{lib}: FDR<{thr}: +{pos} / -{neg}  (n tested={len(res)})")
    combined["_nes"] = pd.to_numeric(combined["NES"], errors="coerce")
    combined["_q"] = pd.to_numeric(combined["FDR q-val"], errors="coerce")
    top = combined.reindex(
        combined["_nes"].abs().sort_values(ascending=False).index
    ).head(12)
    print("\nTop by |NES|:")
    for _, r in top.iterrows():
        print(
            f"  {r['Term'][:55]:55s} {r['Library'][:28]:28s} "
            f"NES={r['_nes']:+.2f} q={r['_q']:.3g} LE={r['leading_edge_size']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
