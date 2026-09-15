# Data description — verified understanding

> Canonical reference for *what the data matrix is* and how it will be preprocessed.
> Generated from evidence (code-tested), with preprocessing decisions confirmed by the
> scientist at the Stage 2 checkpoint (2026-09-15). Regenerate from verified understanding;
> never edit into inconsistency with the data.

## Data-version stamp

| File | sha256 (first 16) | Size | mtime |
|------|-------------------|------|-------|
| `data/protein-data.csv` | `90a0462eb03a02bd` | 1,643,587 B | 2026-09-14 |
| `data/metadata.csv` | `62f10c7072fe613e` | 5504 B | 2026-09-14 |

## Orientation, shape, identifiers

- **Orientation:** rows = **features (proteins)**, columns = **samples**. (Determined, not
  assumed: row index holds UniProt protein ids; column headers are `TPAD####` sample ids.)
- **Shape:** **2343 features × 40 samples.** All 40 columns are experimental samples; **no
  control columns** are present (the 8 pooled references CP/KP are metadata-only — confirms
  Stage-1 caveat **0003**, reconciled both ways: no Pool/QC/Ref/Bridge/Std/Blank-named column
  exists in the data, and no experimental sample is missing from it).
- **Sample-id scheme:** `TPAD####`, all match `TPAD\d+`; these are exactly the metadata
  `SampleID` join key (Stage 1). 0 duplicate sample columns.
- **Feature-id scheme:** UniProt `db|ACCESSION|NAME_ORGANISM`. 2340 `sp|` (SwissProt) + 3
  `tr|` (TrEMBL). Organisms: 2342 `_HUMAN` + **1 `_YEAST`** (`sp|P00924|ENO1_YEAST`). 0
  duplicate feature ids.
- **Protein groups:** **83 features are ambiguous protein groups**, delimited by `" / "`
  (not `;`) — e.g. histone, tubulin, HLA, immunoglobulin-variable families. Policy (below).

## Transformation / normalization / batch state

- **Scale: `log2`.** Values distribute median 22.7, IQR 20.7–25.0, p1–p99 ≈ 16.2–31.9, max
  38.7, no negatives — characteristic log2 intensities. **The loader must record
  `scale = log2`** (so downstream steps refuse scale-incorrect operations, e.g. re-logging).
- **Already normalized (median/quantile).** Per-sample medians span only **0.15 log2 units**
  (22.63–22.78, sd 0.033); per-sample sums CV 0.003. The data arrived normalized.
- **Already batch-corrected on Cohort** (scientist-reported provenance). Applied upstream in
  R via limma:
  ```r
  limma::removeBatchEffect(nobatchProt,
                           design = model.matrix(~Condition, data = nobatchMeta),
                           batch  = nobatchMeta$Cohort)
  ```
  i.e. the **Cohort** effect was removed while **protecting the Condition (diagnostic-group)
  effect** via `design`. **Verified in-data:** within the only groups that span cohorts, the
  residual Kerchner-vs-Poston mean difference is ~0 (HCN mean Δ = −0.002, median |Δ| = 0.000;
  PDD mean Δ = 0.003) — consistent with the stated correction. **Critical caveat →
  Finding 0004:** ADD (100% Kerchner) and PDCN (100% Poston) are **single-cohort**, so the
  cohort effect is **not estimable within those groups**; their correction is **extrapolated**
  from HCN/PDD and cannot resolve the design confound (0001). The pre-correction matrix
  (`nobatchProt`) is **not available in the project** — only the corrected matrix is.

## Missing-value semantics

- The matrix has **no `NaN`, no empty cells, and no true `0`** — it looks 0% missing.
- **BUT a disguised sentinel exists:** exactly **2 features** carry an identical floor value
  **`0.373`** standing for *not-detected / left-censored*, not a real measurement:
  `sp|Q9NX62|IMPA3_HUMAN` (14/40 cells) and `sp|Q5BLP8|CD048_HUMAN` (12/40 cells). No other
  feature contains any value < 13.3. This is the classic "tool token disguised as a number"
  trap; `0.373` is **not** interchangeable with a genuine low intensity.
- Aside from those two features, the matrix is genuinely **complete** (upstream DIA processing
  filtered to features quantified across the cohort). No general imputation is required.

## Contaminants / decoys

- **None present.** 0 `CON__`/contaminant ids, 0 `REV__`/reverse/decoy ids. No handling
  decision needed beyond the yeast spike-in below.

## Preprocessing decisions (confirmed with scientist, 2026-09-15)

Applied by the Stage-3 loader / Stage-4 analysis on the **experimental subset**:

1. **Yeast spike-in → EXCLUDE.** Drop `sp|P00924|ENO1_YEAST` (enolase; stable ~32.8 log2, a
   loading/normalization standard, not human biology) from the analysis matrix. It may be
   viewed in QC as a loading control.
2. **Floored features → DROP BOTH.** Treat `0.373` as not-detected; remove `IMPA3_HUMAN` and
   `CD048_HUMAN` (both >30% floored — low information, and the placeholder would distort
   statistics). Handle via the `lib/common/missing-values` template on the raw matrix, or as
   an explicit feature drop; record dropped features in `provenance.params`.
3. **Normalization → NONE (treat as pre-normalized).** Data arrived median-normalized; apply
   no further normalization. Record `already-normalized`. (Re-applying median would be a
   near no-op; MAD/VSN are inapplicable — they expect a linear matrix, and this is log2.)
4. **Batch correction → ALREADY APPLIED upstream (do not re-apply).** limma
   `removeBatchEffect` on Cohort, protecting Condition (see above). No additional batch
   correction in this pipeline. **Do not** re-run ComBat/removeBatchEffect (double-correction;
   and the confound makes cohort non-identifiable for ADD/PDCN anyway).

**Resulting analysis matrix:** 2343 − 1 (yeast) − 2 (floored) = **2340 human protein features
× 40 experimental samples**, log2, pre-normalized, pre-batch-corrected.

## Protein-group policy

The 83 `" / "`-delimited groups are quantified as a single unit (proteins indistinguishable
by the evidence). Policy: **keep the group as one feature**, annotate/report by its **first
(representative) accession**, and **flag member-ambiguity** on any finding that lands on a
group (the specific member is not resolved). Do not split groups into members.

## Known data-quality issues / caveats

- **0004 (new, this stage):** delivered matrix is limma-batch-corrected on Cohort; correction
  is extrapolated for the single-cohort groups ADD/PDCN and cannot resolve the 0001 confound;
  downstream inference on `removeBatchEffect` output tends to **understate variance** (limma
  authors intend it for visualization/clustering; including batch in the model is the usual
  inferential route, but that is impossible here because Cohort is aliased with the contrast).
  Related to **0001**. Pre-correction matrix unavailable.
- **0001 / 0002 / 0003** (Stage 1) carry forward unchanged at the design level.
- `0.373` sentinel handled by dropping the 2 affected features (decision 2).

## Statistical implications for Stage 4 (carry forward)

- Group means are preserved by the batch correction, so differential abundance on Condition is
  meaningful — but p-values on the corrected matrix may be **anticonservative**; prefer effect
  sizes + CIs, and validate with the label-shuffle null (caveat 0003) and blinded replication.
- Any ADD-specific hit remains **inseparable** from an ADD-specific (Kerchner-specific) site
  effect beyond what HCN/PDD could estimate (0001 + 0004) — qualify accordingly, and use the
  four-arm structure (is a feature dementia-generic? does it match the HCN cohort pattern?).
