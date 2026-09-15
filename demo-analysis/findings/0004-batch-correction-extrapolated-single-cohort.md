---
id: 4
title: "Delivered matrix is limma batch-corrected on cohort — correction is extrapolated for the single-cohort groups (ADD, PDCN) and cannot resolve the design confound"
status: candidate
phase: exploratory
kind: caveat
created: 2026-09-15
updated: 2026-09-15

summary: "The analysis matrix (data/protein-data.csv, 2343 features × 40 samples, log2, median-normalized) was batch-corrected upstream with limma::removeBatchEffect on Cohort while protecting Condition. Within the only two groups that span cohorts the residual within-group Kerchner-vs-Poston mean difference is ~0 (HCN Δ=-0.002, PDD Δ=0.003), but ADD (100% Kerchner) and PDCN (100% Poston) are single-cohort, so their cohort adjustment is extrapolated from HCN/PDD — not estimated — and the design confound (finding 0001) is not resolved."
verdict: "The batch correction preserves group means by design but cannot separate cohort from condition for the single-cohort groups; its ADD/PDCN adjustment is unverifiable extrapolation. Any ADD-specific hit remains inseparable from an ADD/Kerchner site effect beyond what HCN/PDD estimate. Downstream inference on removeBatchEffect output tends to understate variance, so treat p-values on this matrix as anticonservative — prefer effect sizes + CIs, the label-shuffle null, and blinded replication; do not re-apply batch correction."

entities: []

relationships:
  - { type: relates_to, target: 1, note: "The batch correction attempts to mitigate the group×cohort/site confound recorded in 0001 but cannot resolve it, because ADD and PDCN are single-cohort so their cohort adjustment is extrapolated, not estimated." }
  - { type: relates_to, target: 3, note: "The anticonservative-variance risk from running inference on removeBatchEffect output reinforces 0003's requirement for a label-shuffle null and effect-size + CI reporting rather than reliance on p-values." }

provenance:
  data_version: "sha256:90a0462eb03a02bd"
  script: null
  params:
    sample_set: "experimental"
    n_experimental: 40
    verification: "within-group Kerchner-vs-Poston mean-delta after correction (HCN mean Δ=-0.002, median|Δ|=0.000; PDD mean Δ=0.003, median|Δ|=0.000)"
    single_cohort_groups: "ADD=10/10 Kerchner, PDCN=10/10 Poston — cohort effect not estimable within-group"
    cross_cohort_groups: "HCN, PDD — span both site clusters; the only within-condition cross-site estimates available"
    batch_correction_source: "limma::removeBatchEffect(nobatchProt, design = model.matrix(~Condition, data = nobatchMeta), batch = nobatchMeta$Cohort)"
    precorrection_matrix_available: false
    note: "the within-group cohort-delta verification was computed interactively (ad hoc Python) in Stage 2; no script artifact exists — provenance.script is null. This computation must be scripted (and the script promoted) with environment populated before any promotion."
  environment: null
  seeded_from: null
  seed: null
  result_id: null

evidence:
  - metric: "within-group Kerchner-vs-Poston mean Δ (log2, post-correction)"
    value: -0.002
    ci: null
    p_value: null
    p_adjusted: null
    correction: null
    test: "residual cohort mean difference within HCN (the cross-cohort group)"
    n: 10
    note: "HCN spans both cohorts; median |Δ| = 0.000 — consistent with removeBatchEffect having zeroed the within-condition cohort effect"
  - metric: "within-group Kerchner-vs-Poston mean Δ (log2, post-correction)"
    value: 0.003
    ci: null
    p_value: null
    p_adjusted: null
    correction: null
    test: "residual cohort mean difference within PDD (the cross-cohort group)"
    n: 10
    note: "PDD spans both cohorts; median |Δ| = 0.000 — consistent with the stated correction"
  - metric: "cohort composition (single-cohort)"
    value: null
    ci: null
    p_value: null
    p_adjusted: null
    correction: null
    test: "ADD cohort crosstab (see finding 0001)"
    n: 10
    note: "ADD is 10/10 Kerchner — cohort effect not estimable within group; correction is extrapolated from HCN/PDD"
  - metric: "cohort composition (single-cohort)"
    value: null
    ci: null
    p_value: null
    p_adjusted: null
    correction: null
    test: "PDCN cohort crosstab (see finding 0001)"
    n: 10
    note: "PDCN is 10/10 Poston — cohort effect not estimable within group; correction is extrapolated from HCN/PDD"

figures: []

references: []

validation:
  computational_reproduction: { status: not_attempted }
  analytic_replication:       { status: not_attempted }
  data_replication:           { status: not_attempted }

integrity_signoff: true
---

# Delivered matrix is limma batch-corrected on cohort — correction is extrapolated for the single-cohort groups (ADD, PDCN) and cannot resolve the design confound

## Summary
The analysis matrix delivered as `data/protein-data.csv` (2343 features × 40 experimental samples, log2, already median-normalized) was, per the scientist's provenance, batch-corrected upstream in R with `limma::removeBatchEffect(nobatchProt, design = model.matrix(~Condition, data = nobatchMeta), batch = nobatchMeta$Cohort)` — the Cohort effect removed while the Condition (diagnostic-group) effect was protected via `design`. In-data verification confirms the correction zeroed the within-condition cohort effect **where it was estimable**: within the only two groups that span cohorts, the residual Kerchner-vs-Poston mean difference is ~0 (HCN mean Δ = -0.002, median |Δ| = 0.000; PDD mean Δ = 0.003, median |Δ| = 0.000). But ADD is 100% Kerchner and PDCN is 100% Poston (single-cohort), so for those two groups the cohort effect is **not estimable within-group** and their correction is entirely **extrapolated** from HCN/PDD. The design confound recorded in [finding 0001](0001-group-aliased-with-cohort.md) is therefore not resolved by this correction.

## Verdict
`removeBatchEffect` preserves group means by design, but it cannot separate cohort from condition for the single-cohort groups: the ADD/PDCN cohort adjustment is an HCN/PDD-derived extrapolation imposed where it is unverifiable. Any ADD-specific hit remains inseparable from an ADD-specific (Kerchner-specific) site effect beyond what HCN/PDD estimate. Downstream inference run on `removeBatchEffect` output also tends to understate variance, so p-values on this matrix may be anticonservative — prefer effect sizes + CIs, the label-shuffle null (see [finding 0003](0003-contrasts-class-imbalanced-small-n.md)), and blinded replication. Do not re-apply batch correction downstream (double-correction).

## Evidence
The correction worked where the cohort effect was estimable. Within the only two diagnostic groups that span both site clusters, the residual within-group Kerchner-vs-Poston mean difference after correction is essentially zero: **HCN mean Δ = -0.002** (median |Δ| = 0.000) and **PDD mean Δ = 0.003** (median |Δ| = 0.000). This is consistent with the stated `removeBatchEffect` having removed the within-condition cohort effect for the groups where both cohorts are present.

The problem is the two groups where both cohorts are **not** present. ADD is entirely Kerchner (10/10) and PDCN is entirely Poston (10/10), so within either group the cohort effect cannot be estimated at all; whatever adjustment `removeBatchEffect` applied to them is extrapolated from the cross-cohort groups (HCN, PDD) rather than measured. This single-cohort structure is exactly the near-block-diagonal crosstab already recorded and reviewed in [finding 0001](0001-group-aliased-with-cohort.md): its group × cohort sample-count matrix (`figures/metadata/crosstabs/group-x-cohort-confound.png`) shows ADD = 10/10 Kerchner, PDCN = 10/10 Poston, and that only HCN and PDD span cohorts. No new figure is produced here — the claim rests on that already-reviewed crosstab, referenced rather than duplicated.

## Methods / how to produce
The batch-correction provenance is the scientist-reported upstream R step quoted in `state/DATA_DESCRIPTION.md`: `limma::removeBatchEffect(nobatchProt, design = model.matrix(~Condition, data = nobatchMeta), batch = nobatchMeta$Cohort)`. The residual within-group cohort effect was checked ad hoc in Stage 2 on the delivered matrix (`data/protein-data.csv`, `sha256:90a0462eb03a02bd`, mtime 2026-09-14) by computing, within each cross-cohort group, the mean per-feature Kerchner-vs-Poston difference (HCN, PDD). This was an ad-hoc interactive Python computation — **no script artifact exists**, so `provenance.script` is `null`; the computation must be scripted (and the script promoted, with `environment` populated) before any promotion. The single-cohort composition (ADD, PDCN) is read from the group × cohort crosstab in [finding 0001](0001-group-aliased-with-cohort.md). The pre-correction matrix (`nobatchProt`) is **not available** in the project — only the corrected matrix is delivered, so the correction cannot be re-derived or independently re-estimated here.

## Discussion
Because Cohort is aliased with Condition ([finding 0001](0001-group-aliased-with-cohort.md)), `removeBatchEffect` cannot separate the two for ADD and PDCN. Protecting Condition via `design` guarantees the group **means** are preserved, but the cohort adjustment applied to the single-cohort groups is borrowed from the geometry of the cross-cohort groups (HCN, PDD) and is unverifiable for ADD/PDCN. The design confound is therefore not dissolved by correction: an ADD-specific signal remains inseparable from an ADD/Kerchner site effect beyond what HCN/PDD can estimate.

There is a second, statistical hazard. `removeBatchEffect` is intended by the limma authors for visualization and clustering; the inferential route is to include batch as a term in the linear model, which is impossible here precisely because Cohort is aliased with the contrast. Running downstream inference on the `removeBatchEffect` **output** tends to understate residual variance, so p-values computed on this matrix can be anticonservative. The honest downstream posture is effect sizes with confidence intervals, an empirical label-shuffle null, and blinded replication — not naked p-values.

## Caveats
This is a data-provenance / processing caveat surfaced at Stage 2 (data understanding), recorded before the integrity gate; its `integrity_signoff` was certified at the Stage-3 integrity gate (signed off by Michael Riffle on 2026-09-15) for the certified `data_version`, and is now `true`. The evidence is descriptive (a residual within-group mean difference and cohort composition), not a hypothesis test, and makes no held-out claim — hence `phase: exploratory`. The pre-correction matrix is unavailable, so the correction itself cannot be independently reproduced; the within-group residual check verifies its effect only where cohort was estimable. **Do not re-apply batch correction downstream** — double-correction, and cohort is non-identifiable for ADD/PDCN in any case.

## Follow-ups
- Carry this caveat into every downstream group-contrast discovery (alongside [finding 0001](0001-group-aliased-with-cohort.md)): qualify ADD-specific claims as inseparable from an ADD/Kerchner site effect, and read them against the four-arm structure (is a feature dementia-generic in PDD? does it track the HCN cohort pattern?).
- Prefer effect sizes + CIs and the label-shuffle null over p-values on this matrix; treat inference on the corrected matrix as anticonservative.
- If a figure is wanted for the residual-cohort-effect verification (e.g. per-group within-group Kerchner-vs-Poston mean-delta), commission one when the ad hoc check is promoted to a script.

## Related findings
- **relates_to → [finding 0001](0001-group-aliased-with-cohort.md):** the batch correction attempts to mitigate the group × cohort/site confound recorded in 0001 but cannot resolve it, because ADD and PDCN are single-cohort so their cohort adjustment is extrapolated, not estimated. This finding's single-cohort evidence is read directly from 0001's crosstab.
- **relates_to → [finding 0003](0003-contrasts-class-imbalanced-small-n.md):** the anticonservative-variance risk from running inference on `removeBatchEffect` output reinforces 0003's requirement for a label-shuffle null and effect-size + CI reporting rather than reliance on p-values.

## References
None.
