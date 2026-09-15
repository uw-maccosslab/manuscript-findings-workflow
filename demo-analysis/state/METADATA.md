# Metadata — verified understanding

> Canonical reference for *what the experiment is*. Generated from evidence (code-tested,
> not hand-asserted) and confirmed with the scientist at the Stage 1 human checkpoint
> (2026-09-15). Regenerate from verified understanding; never edit into inconsistency with
> the data.

## Data-version stamp

Describes these exact files (drift guard):

| File | sha256 (first 16) | Size | mtime |
|------|-------------------|------|-------|
| `data/metadata.csv` | `62f10c7072fe613e` | 5504 B | 2026-09-14 |
| `data/protein-data.csv` | `90a0462eb03a02bd` | 1,643,587 B | 2026-09-14 |

## Overview

`data/metadata.csv` holds **48 rows** (samples) × **19 columns**. The samples split into
**40 experimental** and **8 control references** (see classification below). The protein
data matrix (`data/protein-data.csv`) contains **only the 40 experimental samples** as
columns; the 8 controls are metadata-only.

## Experimental vs. control sample classification

**Rule:** the split is derived from the **`Condition`** column.

| Class | Condition value(s) | n | Sample IDs | Control subtype |
|-------|--------------------|---|-----------|-----------------|
| **Experimental** (analysis set) | ADD, PDD, PDCN, HCN | 40 | `TPAD10xx` | — |
| Control | `Plasma Reference` | 4 | `CP01`–`CP04` | pooled reference (common/commercial plasma pool) |
| Control | `TPAD plasma pooled reference` | 4 | `KP01`–`KP04` | pooled reference (pool built from the TPAD cohort) |

- The **experimental subset (40 samples) is the analysis set.** The 8 references are
  QC/normalization controls, **excluded from all Stage-4 biological analysis** and viewed
  **separately in QC**. (Confirmed by the scientist.)
- The 8 control references are **present in metadata but absent from the data matrix**
  (`data/protein-data.csv` has 40 columns, all experimental). Consequence: this data file
  carries **no pooled-reference quantitative data**, so technical reproducibility / run
  drift **cannot be assessed from the pooled references in this matrix** — flagged as a
  caveat to confirm at Stage 2/3. See finding **0003** below.

## Experimental design encoded

Four diagnostic groups, **10 samples each** (balanced), banked human EDTA plasma:

| Group | Meaning | Role in design |
|-------|---------|----------------|
| **ADD** | Alzheimer's disease dementia | group of interest |
| **PDD** | Parkinson's disease dementia | controls for dementia-generic effects |
| **PDCN** | Parkinson's disease, cognitively normal | controls for PD / movement-disorder treatment without dementia |
| **HCN** | Healthy cognitively normal | neurodegeneration-free reference |

**No technical replicates.** **Single plate/batch — no acquisition-batch axis** (confirmed
by the scientist); `KF Plate Position` is a well location, not a batch factor.

## Join key to the data matrix

**`SampleID`** (values like `TPAD1019`). The `SampleID` values of the 40 experimental rows
are **exactly** the 40 column headers of `data/protein-data.csv`. Join verified complete:
40/40 map, no orphans on either side; the 8 control `SampleID`s (`CP0x`/`KP0x`) do not
appear in the matrix.

## Columns — inferred meaning and validated domain

| Column | Meaning | Type / domain | Notes |
|--------|---------|---------------|-------|
| `Chris ID` | sample identifier (e.g. `AD01`) | string, unique (48/48) | group-prefixed local ID |
| `SampleID` | sample identifier — **join key** | string, unique (48/48) | `TPAD10xx`; = data-matrix columns |
| `Poston_Plasma_ID` | plasma-bank ID | int-like; NA for 8 controls | |
| `Linkcode` | de-identification linking code | numeric string; NA for 8 controls | |
| `Condition` | diagnostic group + control role | 6 levels (4 experimental + 2 control) | drives the exp/control split |
| `Age` | age at draw, years | float 60.0–87.9; complete for all 40 experimental | |
| `Sex` | sex | {Female, Male}; NA for controls | |
| `isMale` | sex, redundant encoding | {0,1}; **0 inconsistencies** vs `Sex` | |
| `Education` | years of education | int 12–20; **1 experimental value missing** (HCN) | |
| `Freezer Box Position` | storage location | string; NA for controls | handling, not a design factor |
| `KF Plate Position` | plate well location | string, unique (48/48) | **not** a batch axis (single plate) |
| `Run Order` | MS acquisition order | int, unique; range **17–85** | interleaved among a larger run; **not** confounded with group (p=0.98) |
| `Tube Type` | collection tube | {Colored, 0.5 mL Colored, Clear} | handling |
| `Label` | tube label style | 4 levels | handling |
| `Volume (mL)` | plasma volume | {0.5, 1} | handling |
| `Cohort` | source cohort/site | {Kerchner, Poston, ADRC/PUC, PUC only} | **confounded with group — see 0001** |
| `Collection date` | plasma draw date | M/D/YYYY; 2011–2021 | **skewed across groups — see 0002** |
| `alpha-syn` | α-synuclein seed-amplification assay (SAA) | {SSA+, SSA−, not completed} | orthogonal biomarker, not a design factor |
| `plasma ptau` | plasma phospho-tau (p-tau181) | float; NA for controls + 3 experimental | orthogonal biomarker |

## Validity checks (code-tested)

- `SampleID`, `Chris ID`, `KF Plate Position`, `Run Order` all **unique across 48 rows**.
- `Sex` vs `isMale`: **0 inconsistencies**.
- `Age`: numeric, all 40 experimental present, range 60.0–87.9 (plausible).
- Missing values are structural: the 8 controls carry `na` for demographic/handling fields;
  additionally 1 experimental `Education` and 3 experimental `plasma ptau` are missing.

## Relationships tested — imbalance, skew, confounding

Bias-corrected **Cramér's V** for categorical associations; **Kruskal–Wallis** for continuous
covariates across the four groups. Numbers in `results/confounding_stats.json`,
`results/crosstab_group_x_cohort.csv`, `results/table1_cohort_summary.csv`.

**Balanced / not concerning:**
- **Sex** ⟂ group (V = 0.00; 5/5 in ADD, PDCN, HCN; 4/6 PDD).
- **Age** not materially different (Kruskal p = 0.24). **Education** balanced (p = 0.83).
- **Run Order** evenly spread across groups (Kruskal p = 0.98) — acquisition order is **not**
  confounded with group.

**Material gotchas (recorded as caveat findings):**
- 🔴 **Cohort ↔ group confound** — Cramér's V(group×cohort) = **0.515**; for the primary
  contrast **ADD vs others V = 0.62** (ADD is **10/10 Kerchner**), and for the secondary
  contrast **PD vs others V = 0.83** (PD is 17/20 Poston; Other is 17/20 Kerchner). The
  disease contrasts are **aliased with source cohort/collection site**. → **Finding 0001**
  (`group-aliased-with-cohort`, `kind: caveat`). Figure:
  `figures/metadata/crosstabs/group-x-cohort-confound.*`.
- 🟠 **Collection-year skew** — Kruskal p = 0.026; ADD collected earlier (median ~2013) than
  others (~2014), PDD spread wide (to 2021). Travels with the cohort split ⇒ differential
  plasma storage/freezer time. → **Finding 0002** (`collection-year-skews-across-groups`,
  `kind: caveat`; `relates_to` 0001). Figure:
  `figures/metadata/distributions/collection-year-by-group.*`.
- 🟠 **Class imbalance / no pooled-reference data in the matrix** — the primary contrast is
  10 vs 30 (secondary 20 vs 20; only 10 positives), demanding imbalance-aware metrics,
  stratified CV, and a label-shuffle null; and the pooled references (CP/KP) are absent from
  the data matrix so run drift can't be assessed from them here. → **Finding 0003**
  (`contrasts-class-imbalanced-small-n`, `kind: caveat`; `relates_to` 0001; confirm
  reference-data availability at Stage 2/3). Figure:
  `figures/metadata/distributions/cohort-characterization.*`.

*(All three are `status: candidate`, `phase: exploratory`, `integrity_signoff: false` —
their signoff is set at the Stage 3 integrity gate for the certified data version. See
`findings/manifest.md`.)*

## Cohort summary ("Table 1")

Full table: `results/table1_cohort_summary.csv`. Visual: `figures/metadata/distributions/`.

| Group | n | Age mean(sd) | F / M | Educ mean(sd) | Cohort (Kerch/Post/ADRC-PUC/PUC-only) | Collect yr (med) |
|-------|---|--------------|-------|---------------|----------------------------------------|------------------|
| ADD | 10 | 75.7 (7.4) | 5 / 5 | 16.5 (2.5) | 10/0/0/0 | 2013 |
| PDD | 10 | 72.8 (7.1) | 4 / 6 | 16.3 (2.6) | 0/7/2/1 | 2014 |
| PDCN | 10 | 69.2 (4.6) | 5 / 5 | 16.5 (3.2) | 0/10/0/0 | 2014 |
| HCN | 10 | 74.0 (7.2) | 5 / 5 | 17.4 (1.8) | 7/3/0/0 | 2014 |
