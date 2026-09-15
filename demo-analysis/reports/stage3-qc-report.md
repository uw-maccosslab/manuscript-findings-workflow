# Stage 3 — Loaders, pairing & QC report  [INTEGRITY GATE]

**Study:** Plasma-EV proteomics for Alzheimer's disease dementia markers (see `state/PROJECT.md`).
**Date:** 2026-09-15 · **Status:** ✅ **INTEGRITY GATE PASSED** — signed off by Michael Riffle, 2026-09-15.

## Data version certified by this gate

| File | sha256 (first 16) | Shape |
|------|-------------------|-------|
| `data/protein-data.csv` | `90a0462eb03a02bd` | 2343 features × 40 samples |
| `data/metadata.csv` | `62f10c7072fe613e` | 48 rows |

If either file changes, this gate is invalidated and must be re-run.

## 1. Loader — built, tested, independently reviewed

Loader: `scripts/scratch/project_loader.py` (seeds the vetted `wide-data-loader` template
`scripts/common/data_loading.py`). Reconciles the 48-row metadata into 40 experimental /
8 control by an explicit `Condition` allow-list, reconciles **both ways** against the data
columns (fail-loud on any orphan), preserves feature-id strings verbatim, records
`scale=log2`, pairs 40↔40 exactly, and exposes an analysis-subset function.

- **Obligation A — tests pass.** `test_project_loader.py`: **23 tests pass** (ruff + mypy --strict clean).
  Unit fixtures (orientation, pairing, scale, id-fidelity incl. pipe/date-like ids),
  Hypothesis property/invariant tests (every loaded value == source cell; Run-Order ordering),
  a planted-truth mean-difference fixture, and fail-loud edge cases (true orphan both
  directions, duplicate feature id, NaN under missing="error", unknown Condition, missing
  drop target).
- **Obligation B — verified against the source (two independent ways).** An independent
  pandas transpose read reconciles **byte-identically** to the loaded matrix (full-matrix
  `array_equal` = True; 12/12 random spot cells exact). Counts reconcile (40×2343);
  orientation confirmed; dtype float64, no NaN; range 0.373–38.74; samples ordered by
  Run Order (17→85); IDs verbatim.
- **Independent code-review:** APPROVED — no correctness/data-handling defects (two
  registry-hygiene nits fixed). Reconciliation judged airtight both ways.

## 2. Sample↔metadata pairing & experimental/control classification

- **Pairing complete and exact:** all 40 data columns map 1:1 to their experimental metadata
  rows; no orphans on either side.
- **Analysis set = 40 experimental samples** (ADD/PDD/PDCN/HCN, 10 each). The **8 pooled
  references (CP01–04, KP01–04) are metadata-only** — absent from the data matrix
  (reconciled both ways). Consequence: no pooled-reference data exists to assess run drift
  from (caveat **0003**).

## 3. Preprocessing applied (decisions from Stage 2)

| Decision | Applied |
|----------|---------|
| Yeast spike-in `ENO1_YEAST` | dropped (QC control, not human biology) |
| Floored features `IMPA3`, `CD048` (0.373 = not-detected) | dropped (both >30% floored) |
| Normalization | none — already median-normalized (verified: sample medians span 0.15 log2) |
| Batch correction | already applied upstream (limma removeBatchEffect on Cohort, protecting Condition); not re-applied |

**Analysis matrix = 2340 human features × 40 samples**, log2, pre-normalized,
pre-batch-corrected. Canonical QC states materialized once under `results/qc_states/`
(`delivered_log`, `delivered_linear`).

## 4. QC figures (all reviewed & accepted)

Single delivered processing state (log2, normalized, batch-corrected on cohort); the
raw→normalized→batch-corrected series collapses to one state because pre-processing matrices
are unavailable (caveat **0004**). Controls absent (**0003**) → all QC is experimental-only.

| Family | File (`figures/qc/…`) | Reading |
|--------|------------------------|---------|
| Abundance box plots | `abundance-boxplot/` | Per-sample medians flat across the run (drift −0.03 log2) — normalization confirmed. |
| CV | `cv/` | Median per-protein CV 0.397 (de-logged linear); single-state snapshot. |
| Dynamic range | `dynamic-range/` | Spans 6.8 orders of magnitude; smooth, no hyper-abundant cliff (normalized data — albumin at rank 206, not the head). |
| Sample correlation | `sample-correlation/` | Samples cohere (off-diag r 0.80–0.99, median 0.94). **TPAD1012 is a global correlation outlier** (mean r 0.88, ~5 MADs below pack) — candidate QC concern. |
| PCA (by cohort) | `pca/pca-delivered-by-cohort.*` | PC1 (33%): cohorts overlap (p=0.40) — correction worked on the main axis. **PC2 (8.7%): residual cohort split** (p<0.0001) — correction incomplete on the minor axis. |
| PCA (by group) | `pca/pca-delivered-by-group.*` | Groups intermix on PC1 (p=0.26), split on PC2 (p<0.0001) — but that split coincides with the PC2 cohort split, so it **cannot be attributed to biology** (confound 0001/0004). |

## 5. Batch-confounding re-confirmation

`assess_batch_confounding(Cohort, Condition)` on the delivered matrix: bias-corrected
Cramér's V = **0.515**, "strongly confounded" warning fired — reproduces caveat **0001**
exactly. Corroborates 0001/0004; not recorded as a new finding.

## 6. Metadata annotation check (sex & age vs the proteome)

Requested spot-check for misannotated samples. Script `scripts/scratch/qc_sex_age_annotation_check.py`
(leakage-safe LOO CV, label-shuffle null). **Verdict: no sample confidently misannotated or
swapped — but a power limitation, not a clean bill.**
- **Sex:** no Y-chromosome proteins present in this plasma matrix; multivariate signal weak
  (AUC 0.71 [0.53–0.87], perm p=0.017, ~60% accuracy). Two nominal flags (TPAD1010, TPAD1009)
  are ordinary classifier errors, **not** credible swaps. Cannot adjudicate individuals.
- **Age:** no usable signal at n=40 (apparent r is a LOO regularization artifact, sits in the
  null p=0.63); no gross swaps.
- Capped `exploratory`; not recorded as a finding (a null/underpowered check, not a claim).

## 7. EV-enrichment & AD-marker detectability QC

Cross-referenced verified AD/EV markers (research corpus: `ad-plasma-atn-core-biomarkers`,
`ad-plasma-glial-axonal-biomarkers`, `ad-plasma-ev-proteomic-markers`) against the 2340
detected features. Figure: `figures/qc/dynamic-range/dynamic-range-delivered-linear-ad-markers.*`.
- **EV enrichment supported:** 5/7 canonical EV-QC anchors detected — CD9, syntenin (SDCB1),
  ALIX (PDC6I), TSG101, FLOT1 (CD63, CD81 not detected).
- **Classic AD plasma markers mostly undetectable** (ultralow-abundance, below untargeted
  plasma DIA-MS floor): Tau, NfL, NRGN, CHI3L1, GAP43, SNAP25 absent. Only GFAP and APP
  *parent proteins* are present — **a parent-protein peptide, NOT the diagnostic Aβ42/40
  fragment or p-tau phosphopeptide**; do not read as a biomarker hit.
- EV-cargo candidates present: L1CAM (contested), NCAM1, F13A1, F13B.

## 8. Caveats carried into analysis (from Stages 1–2, unchanged)

- **0001** — group aliased with cohort (V=0.62 ADD-vs-others, 0.83 PD-vs-others).
- **0002** — collection-year skew (storage-time confound).
- **0003** — class imbalance (10-vs-30) + pooled references absent from matrix.
- **0004** — delivered matrix limma-batch-corrected on cohort; correction extrapolated for
  single-cohort ADD/PDCN; inference on removeBatchEffect output may understate variance.

At sign-off these Stage-1 caveats' `integrity_signoff` flips to `true` for this data version.

## 9. Sign-off decisions (2026-09-15)

1. **TPAD1012** — global correlation outlier: **KEEP in the analysis set, FLAGGED** (scientist
   decision). Recorded as caveat finding **0005**; Stage-4 analyses that could be influenced by
   a single sample must run a with/without-TPAD1012 sensitivity check. (It is an ADD sample, so
   dropping it would worsen the 10-vs-30 imbalance — relates to 0003.)
2. Annotation check acknowledged as **clean-but-underpowered** (not a positive verification);
   not recorded as a finding.
3. **Integrity gate PASSED** for data version `protein-data.csv sha256:90a0462eb03a02bd;
   metadata.csv sha256:62f10c7072fe613e`, signed off by Michael Riffle. Caveats 0001–0004 (and
   0005) carry `integrity_signoff: true` for this data version. **Stage-4 analysis is unlocked.**
