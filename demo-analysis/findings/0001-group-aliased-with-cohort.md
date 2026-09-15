---
id: 1
title: "Diagnostic group is aliased with source cohort — group contrasts are confounded with collection site"
status: candidate
phase: exploratory
kind: caveat
created: 2026-09-15
updated: 2026-09-15

summary: "Diagnostic group is almost perfectly aliased with collection cohort/site (bias-corrected Cramér's V group×cohort = 0.515; ADD-vs-other = 0.620; PD-vs-other = 0.828). Any protein that separates a diagnostic group from the rest is inseparable from a Kerchner-vs-Poston collection-site/batch effect."
verdict: "The single most consequential design issue in this study: group contrasts cannot be cleanly attributed to disease biology because group and site are confounded. Every downstream group-contrast finding must carry this caveat."

entities: []

relationships: []

provenance:
  data_version: "sha256:62f10c7072fe613e"
  script: { path: "scripts/figures/plot_group_x_cohort_confound.py", commit: "uncommitted" }
  params: { sample_set: "experimental", n_experimental: 40, n_controls_excluded: 8, stats_source: "results/confounding_stats.json", crosstab_source: "results/crosstab_group_x_cohort.csv" }
  environment: null
  seeded_from: null
  seed: null
  result_id: null

evidence:
  - metric: "Cramér's V (bias-corrected)"
    value: 0.515
    ci: null
    p_value: null
    p_adjusted: null
    correction: null
    test: "association of diagnostic group (4 levels) with collection cohort/site (4 levels)"
    n: 40
    note: "group × cohort — overall aliasing of diagnosis with site"
  - metric: "Cramér's V (bias-corrected)"
    value: 0.620
    ci: null
    p_value: null
    p_adjusted: null
    correction: null
    test: "association of PRIMARY contrast (ADD vs other) with cohort/site"
    n: 40
    note: "primary contrast ADD-vs-other × cohort; ADD is 10/10 Kerchner"
  - metric: "Cramér's V (bias-corrected)"
    value: 0.828
    ci: null
    p_value: null
    p_adjusted: null
    correction: null
    test: "association of SECONDARY contrast (PD vs other) with cohort/site"
    n: 40
    note: "secondary contrast PD-vs-other × cohort; PD is 17/20 Poston, Other is 17/20 Kerchner"

figures:
  - png:        "figures/metadata/crosstabs/group-x-cohort-confound.png"
    svg:        "figures/metadata/crosstabs/group-x-cohort-confound.svg"
    legend_png: "figures/metadata/crosstabs/group-x-cohort-confound.legend.png"
    legend_svg: "figures/metadata/crosstabs/group-x-cohort-confound.legend.svg"
    caption:    "Sample-count matrix of diagnostic group (rows) by collection cohort/site (columns) for the experimental subset (n=40; 8 control pools excluded). Cell fill = sample count on a sequential viridis scale (0–10); integer count printed per cell; each row totals 10. Bias-corrected Cramér's V: group×cohort=0.515; ADD-vs-other=0.620; PD-vs-other=0.828."
    script:       { path: "scripts/figures/plot_group_x_cohort_confound.py", commit: "uncommitted" }
    data_version: "sha256:62f10c7072fe613e"
    result_id:    null
    params:       {}

references: []

validation:
  computational_reproduction: { status: not_attempted }
  analytic_replication:       { status: not_attempted }
  data_replication:           { status: not_attempted }

integrity_signoff: true
---

# Diagnostic group is aliased with source cohort — group contrasts are confounded with collection site

## Summary
In the 40-sample experimental subset (4 diagnostic groups × 10: ADD, PDD, PDCN, HCN), diagnostic group is almost perfectly aliased with the collection cohort/site. Bias-corrected Cramér's V for group × cohort is 0.515; it rises to 0.620 for the primary ADD-vs-other contrast and 0.828 for the secondary PD-vs-other contrast. Any protein that separates a diagnostic group from the others is therefore inseparable from a Kerchner-vs-Poston collection-site / batch effect.

## Verdict
This is the single most consequential design issue in the study. A feature that "distinguishes" ADD (or PD) could be a site/handling artifact rather than disease biology, and the two explanations cannot be pulled apart with the data as structured. Every downstream group-contrast finding must carry this caveat.

## Evidence
Diagnostic group and collection cohort/site are near-perfectly aliased, and the aliasing is worst for exactly the contrasts the study cares about: bias-corrected Cramér's V is 0.515 for group × cohort overall, 0.620 for the primary ADD-vs-other split (ADD is 10/10 Kerchner), and 0.828 for the secondary PD-vs-other split (PD is 17/20 Poston; Other is 17/20 Kerchner). The sample-count matrix shows why.

![Sample-count matrix of diagnostic group (rows) by collection cohort/site (columns) for the experimental subset (n=40; 8 control pools excluded). Cell fill = sample count on a sequential viridis scale (0–10); integer count printed per cell; each row totals 10. Bias-corrected Cramér's V: group×cohort=0.515; ADD-vs-other=0.620; PD-vs-other=0.828.](../figures/metadata/crosstabs/group-x-cohort-confound.png)

![Legend for the group-by-cohort confound matrix — cell fill encodes sample count on a sequential viridis scale (0–10).](../figures/metadata/crosstabs/group-x-cohort-confound.legend.png)

*Figure 1. Group × cohort/site sample-count matrix. Produced by `scripts/figures/plot_group_x_cohort_confound.py` (uncommitted) from data `sha256:62f10c7072fe613e`, with statistics from `results/confounding_stats.json` and counts from `results/crosstab_group_x_cohort.csv`.*

The matrix is near block-diagonal: ADD is entirely Kerchner (10/10), PDCN entirely Poston (10/10), PDD almost entirely Poston/PUC (7 Poston, 2 ADRC-PUC, 1 PUC-only), and only HCN spans both site clusters (7 Kerchner, 3 Poston). Diagnostic group is therefore almost perfectly aliased with collection site — any group contrast is inseparable from a site effect.

## Methods / how to produce
The confounding statistics (bias-corrected Cramér's V for group × cohort and for the primary/secondary pooled contrasts) are in `results/confounding_stats.json`; the group × cohort contingency counts are in `results/crosstab_group_x_cohort.csv`. The figure is rendered by `scripts/figures/plot_group_x_cohort_confound.py` (scratch) from `data/metadata.csv` (`sha256:62f10c7072fe613e`, mtime 2026-09-14), `results/confounding_stats.json`, and `state/color_registry.json`, restricted to the 40 experimental samples (8 pooled-reference controls excluded). See `provenance`.

## Discussion
When the covariate that defines the biological contrast is aliased with a technical/collection covariate, the two are not separately estimable: a marker that appears to distinguish ADD or PD from other groups is equally consistent with a Kerchner-vs-Poston site, batch, or handling artifact. The four-arm design offers partial leverage — a feature that also moves in PDD points to a dementia-generic or site effect rather than a disease-specific one, and a feature that tracks site across arms is suspect — but the aliasing cannot be dissolved by analysis alone.

## Caveats
This is a metadata (Stage 1) caveat recorded before the integrity gate; its `integrity_signoff` was certified at the Stage-3 integrity gate (signed off by Michael Riffle on 2026-09-15) for the certified `data_version`, and is now `true`. The evidence is descriptive (association statistics on the cohort as observed), not a hypothesis test, and makes no held-out claim — hence `phase: exploratory`.

**Mitigations to carry downstream:** (1) exploit the four-arm structure — check whether a candidate feature also separates PDD (dementia-generic) or tracks site rather than diagnosis; (2) run within-cohort checks where sample sizes allow — HCN is the only group that spans both site clusters, so it is the only within-group cross-site comparison available.

## Follow-ups
- Attach a `relates_to` edge from every downstream group-contrast discovery to this caveat so the confound travels into analysis and the final report's Limitations.
- Where feasible, add site as a modeled covariate or restrict to within-site comparisons, recognizing that full deconfounding is not possible with this design.

## Related findings
None yet. Downstream group-contrast discoveries should link here with a `relates_to` edge (§2.6). The collection-year skew caveat compounds this one (see that finding's `relates_to` edge into this record).

## References
None.
</content>
</invoke>
