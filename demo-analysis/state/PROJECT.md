# Project: Plasma EV proteomics for Alzheimer's disease dementia markers

> Living document — captures the scientific framing established at Stage 0. It will be
> updated as understanding deepens through later stages.

## Domain

Biomarker discovery in **neurodegenerative disease** using **plasma extracellular vesicle
(EV) proteomics**. The aim is to assess the molecular/diagnostic potential of plasma EVs
for **Alzheimer's disease dementia (ADD)**.

- **Sample type:** banked human **EDTA plasma** — *not* CSF, brain tissue, or
  neuron-specific immunocaptured EVs. This matters for interpretation: the EV proteome
  here is whole-plasma-derived, and markers should be read as accessible-in-plasma signals,
  not necessarily brain- or neuron-origin.
- **EV enrichment:** **Mag-Net** method (magnetic bead–based EV enrichment).
- **Measurement:** quantitative **proteomics by data-independent acquisition (DIA)** mass
  spectrometry. Measured features are protein-level quantities.
- **Cohort source:** Stanford longitudinal neurodegenerative-disease cohorts.

## What is being examined and why

A **pilot study** evaluating whether the Mag-Net plasma-EV method can surface molecular
markers that distinguish **ADD** from other clinical groups. The scientific motivation is to
find protein markers associated **specifically with ADD**, as opposed to nonspecific effects
of dementia in general or of Parkinson's disease and its treatments.

The four-group design is deliberately built to isolate an ADD-specific signature from
confounds:

- **ADD** — Alzheimer's disease dementia (the group of interest).
- **PDD** — Parkinson's disease dementia. Controls for **dementia-associated** lifestyle
  changes and symptomatic treatments that could confound ADD-specific signatures.
- **PDCN** — Parkinson's disease, cognitively normal. Controls for **Parkinson's disease
  and movement-disorder treatments** *without* dementia.
- **HCN** — healthy cognitively normal. Neurodegenerative-disease-free reference group.

Reading the arms together lets an ADD signal be separated from "dementia-generic" and
"PD/treatment-generic" effects.

## Experimental design

- **40 individuals, 10 per group**, balanced across ADD / PDD / PDCN / HCN.
- Groups were **approximately age-matched and balanced by sex** (to be verified against the
  metadata in Stage 1).
- **No technical replicates.**
- **Physical sample prep and MS acquisition details** (instrument, timing) are considered
  out of scope for this analysis by the scientist.
- **Batch structure, QC/pooled-reference samples, and acquisition order** were explicitly
  **deferred to the metadata stage (Stage 1)** for examination — do not assume single-batch
  or randomized-order until verified there.

## Scientific goals

Both a **primary** and a **secondary** contrast, each pursued as **discovery + binary
classification** (not pre-specified hypothesis-testing of named candidates):

### Primary — ADD vs. all others
- **Contrast:** ADD (n=10) vs. all other samples pooled (PDD + PDCN + HCN, n=30).
- **Aims:**
  1. Build a **binary classifier** predicting ADD vs. all others.
  2. Identify **features (proteins) that distinguish ADD** — from the classifier, from
     differential abundance analysis, or by other means.

### Secondary — Parkinson's vs. all others
- **Contrast:** Parkinson's (PDD + PDCN, n=20) vs. all others (ADD + HCN, n=20).
- **Aims:** the same — a binary classifier and distinguishing features for the PD grouping.

### Design tension to keep visible
The primary contrast pools PDD/PDCN/HCN as "others," which sets aside the carefully
constructed confound-control structure of the four arms. Where a distinguishing feature
emerges, it will be worth checking whether it is **ADD-specific** or merely
**dementia-generic** (i.e. also elevated in PDD) — the four-group structure supports that
follow-up and should be used to qualify claims of ADD specificity.

### Practical notes for later stages
- **Class imbalance:** primary contrast is 10 vs. 30; classification and validation must
  account for this (appropriate metrics, CV matched to the target, label-shuffle null).
- **Small n:** 10 per group is a pilot — power is limited; effect sizes with CIs and honest
  power statements matter more than bare p-values.

## Guard against motivated reasoning

Stating goals up front sharpens relevance but **invites motivated reasoning** — the hope of
finding an ADD-specific marker can bias what gets reported. These goals are recorded
faithfully as *aims*, not as expected conclusions. The workflow's countermeasures stay in
force: skepticism is generous during exploration but **ruthless at promotion**, and every
promoted finding passes **blind independent validation**. Keep this tension visible when
interpreting results.
