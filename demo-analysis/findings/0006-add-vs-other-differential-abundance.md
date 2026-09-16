---
id: 6
title: "ADD-vs-others differential abundance: a large but anticonservative, cohort-confounded hit list dominated by intracellular proteins — and no known plasma/EV AD marker is differential"
status: candidate
phase: exploratory
kind: discovery
created: 2026-09-15
updated: 2026-09-15

summary: "A moderated (limma eBayes) linear model contrasting ADD (n=10) vs pooled Other {PDD,PDCN,HCN} (n=30) on the delivered 40×2340 log2 plasma-EV matrix returns 363 proteins at BH-q<0.05 (133 up in ADD, 230 down) and 594 at q<0.10, effect range log2FC −2.90…+3.52. The strongest hits are intracellular/organellar proteins (kinases, proteasome, mitochondrial, vesicle-trafficking), NOT canonical secreted or EV-cargo proteins, and none of the known plasma/EV AD reference markers reaches significance."
verdict: "A large ADD-specific hit list, but the p-value histogram shows the anticonservative signature (spike at p≈0 with a suppressed null tail, π0≈0.53), so the 363/594 hit COUNT is inflated — trust the effect sizes + CIs, not the q-values. The effect is inseparable from an ADD/Kerchner site effect (single-cohort, no cohort covariate possible) and runs on already-batch-corrected data (variance deflated), so this is hypothesis-generating only. Two honesty checks are load-bearing: the hits are intracellular proteins rather than EV cargo, and every known plasma/EV AD marker is non-significant. A label-shuffle-nulled classifier and blinded/held-out replication are the honest next tests."

entities:
  - { db: uniprot, id: "P36507", label: "MP2K2/MEK2" }
  - { db: uniprot, id: "O60763", label: "USO1" }
  - { db: uniprot, id: "Q9P253", label: "VPS18" }
  - { db: uniprot, id: "P29350", label: "PTN6/SHP1" }
  - { db: uniprot, id: "P19971", label: "TYPH/TP" }
  - { db: uniprot, id: "Q9Y262", label: "EIF3L" }
  - { db: uniprot, id: "P26440", label: "IVD" }
  - { db: uniprot, id: "P51812", label: "KS6A3/RSK2" }
  - { db: uniprot, id: "P19404", label: "NDUV2/NDUFV2" }
  - { db: uniprot, id: "Q16186", label: "ADRM1" }
  - { db: uniprot, id: "O00560", label: "SDCB1/SDCBP/syntenin-1" }
  - { db: uniprot, id: "P21926", label: "CD9" }
  - { db: uniprot, id: "O75955", label: "FLOT1" }
  - { db: uniprot, id: "Q99816", label: "TS101/TSG101" }
  - { db: uniprot, id: "Q8WUM4", label: "PDC6I/ALIX" }
  - { db: uniprot, id: "P14136", label: "GFAP" }
  - { db: uniprot, id: "P05067", label: "A4/APP" }
  - { db: uniprot, id: "P32004", label: "L1CAM" }
  - { db: uniprot, id: "P13591", label: "NCAM1" }
  - { db: uniprot, id: "P00488", label: "F13A/F13A1" }
  - { db: uniprot, id: "P05160", label: "F13B" }

relationships:
  - { type: relates_to, target: 1, note: "ADD is 100% Kerchner, so the ADD-vs-Other effect is inseparable from a collection-site effect; no cohort covariate is possible (single-cohort → rank-deficient design), exactly the group×cohort aliasing recorded in 0001." }
  - { type: relates_to, target: 4, note: "DE is run on the already-removeBatchEffect-corrected matrix, which deflates within-feature variance → anticonservative p/q; the ideal both-ways robustness check (uncorrected + batch-as-covariate) is impossible because the pre-correction matrix is unavailable (0004)." }

provenance:
  data_version: "sha256:90a0462eb03a02bd"
  script: { path: "scripts/scratch/de_add_vs_other.py", commit: "uncommitted" }
  params:
    sample_set: "experimental"
    n_experimental: 40
    n_controls_excluded: 0
    n_add: 10
    n_other: 30
    contrast: "ADD_vs_pooled_Other"
    other_groups: "PDD, PDCN, HCN (pooled)"
    reference_level: "Other"
    effect: "log2FC = ADD - Other"
    method: "moderated (limma-style empirical-Bayes linear model)"
    covariates: "none"
    cohort_covariate: "impossible — ADD single-cohort, aliased with contrast (rank-deficient); see finding 0001"
    fdr: 0.05
    correction: "BH"
    n_features_tested: 2340
    normalization: "median (upstream, on delivered matrix)"
    log2: true
    batch_correct: { method: "limma::removeBatchEffect", batch_column: "Cohort", covariate_preserved: "Condition", applied: "upstream on delivered matrix; not re-applied", precorrection_matrix_available: false }
    null_model: "none — no shuffle null (univariate DE family; eBayes moderation is the variance control); honesty check is the p-value-histogram calibration read"
  environment: "uv.lock@uncommitted"
  seeded_from: { template: "differential-abundance", version: "0.1" }
  seed: null
  result_id: null

evidence:
  - metric: "hit count, BH-q<0.05"
    value: 363
    ci: null
    p_value: null
    p_adjusted: 0.05
    correction: "BH"
    test: "limma-style moderated (eBayes) linear model, ADD vs pooled Other"
    n: 40
    note: "363 proteins at q<0.05 (133 up in ADD, 230 down); 594 at q<0.10. COUNT is inflated — the p-value null tail is suppressed (π0≈0.53), so read effect sizes + CIs, not this count."
  - metric: "log2FC"
    value: -1.256
    ci: [-1.638, -0.874]
    p_value: 5.18e-08
    p_adjusted: 1.21e-04
    correction: "BH"
    test: "limma moderated t"
    n: 40
    note: "MP2K2/MEK2 (P36507) — top hit by q; down in ADD. Intracellular kinase, not an EV-cargo protein."
  - metric: "log2FC"
    value: -0.776
    ci: [-1.022, -0.529]
    p_value: 1.31e-07
    p_adjusted: 1.54e-04
    correction: "BH"
    test: "limma moderated t"
    n: 40
    note: "USO1 (O60763) — vesicle-trafficking; down in ADD."
  - metric: "log2FC"
    value: 1.700
    ci: [1.097, 2.303]
    p_value: 1.16e-06
    p_adjusted: 9.01e-04
    correction: "BH"
    test: "limma moderated t"
    n: 40
    note: "VPS18 (Q9P253) — vesicle-trafficking; up in ADD."
  - metric: "log2FC"
    value: -0.885
    ci: [-1.209, -0.561]
    p_value: 2.11e-06
    p_adjusted: 1.23e-03
    correction: "BH"
    test: "limma moderated t"
    n: 40
    note: "PTN6/SHP1 (P29350) — intracellular phosphatase; down in ADD."
  - metric: "log2FC"
    value: -1.323
    ci: [-1.850, -0.797]
    p_value: 8.61e-06
    p_adjusted: 1.37e-03
    correction: "BH"
    test: "limma moderated t"
    n: 40
    note: "TYPH/thymidine phosphorylase (P19971); down in ADD."
  - metric: "log2FC"
    value: 1.740
    ci: [1.046, 2.434]
    p_value: 9.04e-06
    p_adjusted: 1.37e-03
    correction: "BH"
    test: "limma moderated t"
    n: 40
    note: "EIF3L (Q9Y262) — translation-initiation; up in ADD."
  - metric: "log2FC"
    value: 1.431
    ci: [0.850, 2.012]
    p_value: 1.21e-05
    p_adjusted: 1.37e-03
    correction: "BH"
    test: "limma moderated t"
    n: 40
    note: "IVD (P26440) — mitochondrial; up in ADD."
  - metric: "log2FC"
    value: -2.905
    ci: [-4.054, -1.755]
    p_value: 7.99e-06
    p_adjusted: 1.37e-03
    correction: "BH"
    test: "limma moderated t"
    n: 40
    note: "KS6A3/RSK2 (P51812) — intracellular kinase; strongest negative effect among top hits; down in ADD."
  - metric: "log2FC"
    value: -1.844
    ci: [-2.598, -1.091]
    p_value: 1.35e-05
    p_adjusted: 1.37e-03
    correction: "BH"
    test: "limma moderated t"
    n: 40
    note: "NDUV2/NDUFV2 (P19404) — mitochondrial complex I; down in ADD."
  - metric: "log2FC"
    value: 3.517
    ci: [2.002, 5.031]
    p_value: 3.01e-05
    p_adjusted: 2.27e-03
    correction: "BH"
    test: "limma moderated t"
    n: 40
    note: "ADRM1 (Q16186) — proteasome; STRONGEST effect in the whole contrast (+3.52); up in ADD."
  - metric: "log2FC"
    value: -0.510
    ci: [-0.883, -0.137]
    p_value: null
    p_adjusted: 0.052
    correction: "BH"
    test: "limma moderated t"
    n: 40
    note: "SDCBP/syntenin-1 (O00560) — closest known plasma/EV AD marker to significance, still q>0.05."
  - metric: "log2FC"
    value: -0.465
    ci: [-0.979, 0.048]
    p_value: null
    p_adjusted: 0.192
    correction: "BH"
    test: "limma moderated t"
    n: 40
    note: "CD9 (P21926) — canonical EV tetraspanin marker; non-significant. Reference markers CD9 0.19, FLOT1 0.29, TSG101 0.38, ALIX 0.70, GFAP 0.61, APP 0.95, L1CAM 0.82, NCAM1 0.94, F13A1 0.56, F13B 0.76 — none significant."
  - metric: "p-value histogram calibration (bin[0], p<0.05)"
    value: 6.504
    ci: null
    p_value: null
    p_adjusted: null
    correction: null
    test: "density vs uniform, 2340 features"
    n: 40
    note: "Spike at p≈0 is 6.5× uniform BUT the null tail (p>0.5) is suppressed to ~0.53× uniform (π0≈0.53) — the anticonservative signature, NOT a clean U-shape. Consequence: hit count inflated; trust effect sizes + CIs, not q."

figures:
  - png:        "figures/analysis/differential-abundance/add-vs-other/volcano-add-vs-other-moderated-q05.png"
    svg:        "figures/analysis/differential-abundance/add-vs-other/volcano-add-vs-other-moderated-q05.svg"
    legend_png: "figures/analysis/differential-abundance/add-vs-other/volcano-add-vs-other-moderated-q05.legend.png"
    legend_svg: "figures/analysis/differential-abundance/add-vs-other/volcano-add-vs-other-moderated-q05.legend.svg"
    caption:    "Volcano, ADD vs pooled Other, moderated (eBayes) model: x = log2FC (ADD − Other), y = −log10(BH q); n=40 (ADD 10, Other 30). Points above q=0.05 colored by direction (higher in ADD vs higher in Other); the 11 known plasma/EV AD reference markers overlaid in a distinct style, all at/below the q=0.05 line. Exploratory — ADD=100% Kerchner site confound (findings 0001/0004); q-values anticonservative and hit count inflated (read effect sizes, not q)."
    script:       { path: "scripts/figures/plot_volcano_add_vs_other.py", commit: "uncommitted" }
    data_version: "sha256:90a0462eb03a02bd"
    result_id:    null
    params:       { input: "results/differential_abundance/add-vs-other-moderated/contrast_table.csv", q_threshold: 0.05 }
  - png:        "figures/analysis/differential-abundance/add-vs-other/pvalue-hist-add-vs-other.png"
    svg:        "figures/analysis/differential-abundance/add-vs-other/pvalue-hist-add-vs-other.svg"
    legend_png: "figures/analysis/differential-abundance/add-vs-other/pvalue-hist-add-vs-other.legend.png"
    legend_svg: "figures/analysis/differential-abundance/add-vs-other/pvalue-hist-add-vs-other.legend.svg"
    caption:    "Calibration check — p-value histogram of the 2340-feature moderated contrast: spike at p≈0 (bin[0] 6.5× uniform) with a suppressed null tail (p>0.5 at ~0.53× uniform, π0≈0.53). This is the anticonservative signature (variance deflation from batch correction 0004 + confound 0001), NOT a clean U-shape; the BH hit count is therefore inflated."
    script:       { path: "scripts/figures/plot_pvalue_hist_add_vs_other.py", commit: "uncommitted" }
    data_version: "sha256:90a0462eb03a02bd"
    result_id:    null
    params:       { input: "results/differential_abundance/add-vs-other-moderated/contrast_table.csv" }

references: []
  # FLAG (see report): the "known plasma/EV AD reference marker" and MISEV2018 EV-QC
  # designations, and the Nielsen-2021 unbiased-candidate class, are background/interpretive
  # claims (conventions/findings.md §2.5) drawn from the research corpus
  # `ad-plasma-ev-proteomic-markers`. They must be added here with proper ids and
  # fact-checked by the research-reviewer before promotion. Left empty rather than
  # fabricating identifiers.

validation:
  computational_reproduction: { status: not_attempted }
  analytic_replication:       { status: not_attempted }
  data_replication:           { status: not_attempted }

integrity_signoff: true
---

# ADD-vs-others differential abundance: a large but anticonservative, cohort-confounded hit list dominated by intracellular proteins — and no known plasma/EV AD marker is differential

## Summary
The primary Stage-4 contrast is a moderated (limma-style empirical-Bayes) linear model of **ADD (n=10)** vs a **pooled Other = {PDD, PDCN, HCN} (n=30)** group, run on the delivered 40×2340 log2 plasma-EV matrix (already median-normalized and batch-corrected on cohort). The effect is `log2FC = ADD − Other` with a 95% CI, raw p, and BH-q per protein. It returns **363 proteins at BH-q<0.05** (133 up in ADD, 230 down) and **594 at q<0.10**, over an effect range of log2FC **−2.90 … +3.52**. Two honesty observations dominate the reading: the strongest hits are **intracellular / organellar** proteins (kinases, proteasome, mitochondrial, vesicle-trafficking), not canonical secreted plasma or EV-cargo proteins; and **none of the known plasma/EV AD reference markers reaches significance**.

## Verdict
This is a large but fragile hit list. The p-value histogram shows the **anticonservative signature** — a spike at p≈0 with a *suppressed* null tail (π0≈0.53), not a clean U-shape — so the 363/594 hit **count is inflated**: read the effect sizes and CIs, not the q-values. The effect is also inseparable from an **ADD/Kerchner collection-site effect** (ADD is single-cohort, so no cohort covariate is possible; [finding 0001](0001-group-aliased-with-cohort.md)) and is computed on **already-batch-corrected data** whose variance is deflated ([finding 0004](0004-batch-correction-extrapolated-single-cohort.md)). It is strictly **hypothesis-generating**. The honest next tests are a **label-shuffle-nulled classifier** and **blinded / held-out replication**.

## Evidence
The contrast returns a large hit list with a wide effect range. 363 proteins clear BH-q<0.05 (133 up in ADD, 230 down) and 594 clear q<0.10; the top hits by q are all intracellular/organellar — MP2K2/MEK2 (log2FC −1.26, 95% CI [−1.64, −0.87], q=1.2e-4), USO1 (−0.78 [−1.02, −0.53]), VPS18 (+1.70 [1.10, 2.30]), PTN6/SHP1 (−0.89 [−1.21, −0.56]), TYPH (−1.32 [−1.85, −0.80]), EIF3L (+1.74 [1.05, 2.43]), IVD (+1.43 [0.85, 2.01]), KS6A3/RSK2 (−2.90 [−4.05, −1.76]), NDUV2/NDUFV2 (−1.84 [−2.60, −1.09]) — and the strongest single effect is the proteasome subunit **ADRM1 (+3.52 [2.00, 5.03], q=2.3e-3)**. These are kinases (MEK2, RSK2, SHP1), proteasome (ADRM1), mitochondrial (NDUFV2, IVD), and vesicle-trafficking (VPS18, USO1) proteins — the intracellular/organellar interior, not secreted or EV-cargo species.

The volcano shows both the size of the list and, crucially, where the known EV/AD markers land:

![Volcano, ADD vs pooled Other, moderated model — x = log2FC (ADD − Other), y = −log10(BH q), n=40; hits above q=0.05 colored by direction, 11 known plasma/EV AD markers overlaid.](../figures/analysis/differential-abundance/add-vs-other/volcano-add-vs-other-moderated-q05.png)

![Legend for Figure 1 — direction of significant hits (higher in ADD vs higher in Other) and the three known-marker highlight classes.](../figures/analysis/differential-abundance/add-vs-other/volcano-add-vs-other-moderated-q05.legend.png)

*Figure 1. Volcano of ADD vs pooled Other, moderated model. Produced by `scripts/figures/plot_volcano_add_vs_other.py` (uncommitted) from `results/differential_abundance/add-vs-other-moderated/contrast_table.csv`, data `sha256:90a0462eb03a02bd`.*

Each point is one of 2340 proteins, positioned by log2FC (x) against BH significance (y); the colored points above the q=0.05 line are the 363 hits, split by direction. The labeled extremes — MEK2, RSK2, ADRM1, NDUFV2 and the rest — are intracellular proteins, not EV cargo. The overlaid known plasma/EV AD markers all sit **at or below** the q=0.05 line: **none is significant** — CD9 (q=0.19), FLOT1 (0.29), TSG101 (0.38), ALIX (0.70), syntenin/SDCBP (0.052, the closest), GFAP (0.61), APP (0.95), L1CAM (0.82), NCAM1 (0.94), F13A1 (0.56), F13B (0.76). The reference panel is exactly where a null hypothesis would put it, which is what makes the intracellular-dominated hit list suspect as EV biology.

The calibration check is the main honesty diagnostic available without a null (no shuffle null is run for univariate DE):

![p-value histogram of the 2340-feature contrast — spike at p≈0 with a suppressed null tail (π0≈0.53).](../figures/analysis/differential-abundance/add-vs-other/pvalue-hist-add-vs-other.png)

![Legend for Figure 2 — observed p-value density vs the uniform-null reference line.](../figures/analysis/differential-abundance/add-vs-other/pvalue-hist-add-vs-other.legend.png)

*Figure 2. p-value calibration check. Produced by `scripts/figures/plot_pvalue_hist_add_vs_other.py` (uncommitted) from the same `contrast_table.csv`, data `sha256:90a0462eb03a02bd`.*

The first bin (p<0.05) stands at **6.5× uniform** — there is real small-p signal — but the null tail (p>0.5) is **depleted to ~0.53× uniform** (π0≈0.53). A well-calibrated mix of true nulls and true effects gives a flat tail at 1× with a spike; here the tail is *below* 1×, the anticonservative signature of variance deflation. Read this as: some genuine signal, but the BH hit **count is inflated** because the null is being under-counted — so effect sizes and CIs are the trustworthy currency, not the q-values.

## Methods / how to produce
Run `scripts/scratch/de_add_vs_other.py` (a call site for the project copy of the `differential-abundance` lib template, seeded from `lib/analysis/differential-abundance` @ 0.1) on the delivered analysis matrix `data/protein-data.csv` (`sha256:90a0462eb03a02bd`, 40 samples × 2340 tested features; controls already excluded and the yeast spike-in + two floored sentinels dropped upstream). The script builds the binary primary factor **ADD vs Other = {PDD, PDCN, HCN} pooled**, runs `differential_abundance(..., method="moderated")` with `reference="Other"` so `log2FC = ADD − Other`, and writes the full per-feature contrast table (effect + 95% CI + p + BH-q + mean abundance) to `results/differential_abundance/add-vs-other-moderated/contrast_table.csv` plus a JSON summary carrying the histogram calibration read and hit counts. **No cohort covariate** is fitted — ADD is 100% Kerchner, so Cohort is aliased with the contrast and the design would be rank-deficient. **No shuffle null** is run (that is the classifier/regressor discipline; the eBayes moderation is the variance control here). The environment is pinned by `uv.lock`; the analysis is deterministic (no RNG). Figures 1–2 are pure `load → plot` jobs over the cached `contrast_table.csv` (`scripts/figures/plot_volcano_add_vs_other.py`, `scripts/figures/plot_pvalue_hist_add_vs_other.py`); they recompute nothing.

**Not promotable as-is:** the analysis script is under `scripts/scratch/` (not code-reviewed, not promoted). This finding therefore **cannot be validated** until the script is reviewed and moved to `scripts/promoted/` (`conventions/findings.md` §2.2, §4).

## Discussion
The biology to be honest about is what the top of the list is *made of*. A plasma-EV proteomics contrast that is supposed to surface secreted or vesicle-cargo differences instead returns intracellular kinases (MEK2, RSK2, SHP1), a proteasome subunit (ADRM1, the single strongest effect), mitochondrial proteins (NDUFV2, IVD), and trafficking machinery (VPS18, USO1) — while the canonical EV/AD reference markers (CD9, FLOT1, TSG101, ALIX, syntenin, GFAP, APP, L1CAM, NCAM1) sit at the null. Read together with the anticonservative histogram and the single-cohort confound, the most parsimonious reading is that much of this list reflects site/processing structure and inflated significance rather than ADD EV biology. That does not make it worthless — there is real small-p signal — but it caps what may be claimed. (The marker-panel and EV-QC designations invoked here are background claims that require the research-reviewer's `references`; see Caveats and the report flag.)

## Caveats
- **Exploratory, small-n, imbalanced, confounded.** n=10 vs 30; the effect was derived on the same data it is observed in — the multiplicity context (this is the primary of several Stage-4 contrasts explored) is in `findings/exploration-log.md`. Never write this up with confirmatory confidence (`conventions/findings.md` §5).
- **Site confound (relates_to [finding 0001](0001-group-aliased-with-cohort.md)).** ADD is 100% Kerchner, so the ADD effect is inseparable from a collection-site effect and **no cohort covariate is possible** — adding Cohort makes the design singular (rank-deficient). Every ADD-specific claim must be qualified as inseparable from an ADD/Kerchner site effect.
- **Anticonservative p/q from batch-corrected input (relates_to [finding 0004](0004-batch-correction-extrapolated-single-cohort.md)).** The matrix is already `limma::removeBatchEffect`-corrected on Cohort; running inference on that output deflates within-feature variance, which is exactly the suppressed-null-tail signature in Figure 2. The ideal both-ways robustness check — uncorrected matrix with batch-as-covariate — is **impossible** because the pre-correction matrix is unavailable. Consequence: **the 363/594 hit count is inflated; trust effect sizes + CIs, not q-values.**
- **Entity ids** are the UniProt accessions carried in the contrast table's feature identifiers; the research-reviewer should confirm the symbol↔accession labels at promotion.

## Follow-ups
- Run a **label-shuffle-nulled classifier** (ADD vs Other) — the honest test that also serves as the leakage check the univariate DE cannot provide.
- Pursue **blinded analytic replication** and, for any `confirmatory` claim, **held-out / orthogonal data replication** (none available in-project yet).
- **Promote the analysis script** (`scripts/scratch/de_add_vs_other.py` → code review → `scripts/promoted/`) before any validation attempt.
- **Add and fact-check `references`** for the known plasma/EV AD marker panel and the MISEV2018 EV-QC / Nielsen-2021 candidate designations (research-reviewer).
- Read the ADD hits against the **four-arm structure** (is a top hit dementia-generic in PDD? does it track the HCN/site pattern?) to separate ADD-specific signal from site/processing structure.
- Normalize the full hit list into `entities` (only the named notable proteins and reference markers are captured here).

## Related findings
- **relates_to → [finding 0001](0001-group-aliased-with-cohort.md):** diagnostic group is aliased with source cohort. ADD is 100% Kerchner, so this contrast's effect is inseparable from a collection-site effect and no cohort covariate is admissible.
- **relates_to → [finding 0004](0004-batch-correction-extrapolated-single-cohort.md):** the matrix is already batch-corrected on cohort; DE on that output deflates variance and yields anticonservative p/q (the suppressed null tail of Figure 2), and the both-ways robustness check is impossible because the pre-correction matrix is unavailable.

## References
None yet — the marker-panel / EV-QC background claims require `references` added and fact-checked by the research-reviewer before promotion (see Caveats and Follow-ups).
</content>
</invoke>
