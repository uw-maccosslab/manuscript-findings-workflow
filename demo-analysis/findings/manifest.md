---
schema_version: 1
generated: 2026-09-15
next_id: 6
engine_version: 0.5.1
---

# Findings manifest

Derived index of the findings graph — regenerable from the finding files. One row per finding.

| ID | Slug | Title | Status | Phase | Kind | Entities | Relationships | Updated | Data version |
|----|------|-------|--------|-------|------|----------|---------------|---------|--------------|
| 1 | group-aliased-with-cohort | Diagnostic group is aliased with source cohort — group contrasts are confounded with collection site | candidate | exploratory | caveat | — | — | 2026-09-15 | sha256:62f10c7072fe613e |
| 2 | collection-year-skews-across-groups | Plasma collection year skews across diagnostic groups — differential storage time confounds group contrasts | candidate | exploratory | caveat | — | relates_to:1 | 2026-09-15 | sha256:62f10c7072fe613e |
| 3 | contrasts-class-imbalanced-small-n | Primary/secondary contrasts are class-imbalanced and small-n; pooled references absent from the data matrix | candidate | exploratory | caveat | — | relates_to:1 | 2026-09-15 | sha256:62f10c7072fe613e |
| 4 | batch-correction-extrapolated-single-cohort | Delivered matrix is limma batch-corrected on cohort — correction is extrapolated for the single-cohort groups (ADD, PDCN) and cannot resolve the design confound | candidate | exploratory | caveat | — | relates_to:1, relates_to:3 | 2026-09-15 | sha256:90a0462eb03a02bd |
| 5 | tpad1012-sample-correlation-outlier | TPAD1012 is a global sample-correlation outlier — retained in the analysis set but flagged for sensitivity analysis | candidate | exploratory | caveat | — | relates_to:3 | 2026-09-15 | sha256:90a0462eb03a02bd |
