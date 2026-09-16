---
id: 8
title: "Boruta all-relevant selection confirms 40 proteins (18 up / 22 down in ADD) for ADD-vs-others — 33 shared with the elastic-net minimal-optimal set plus 7 redundant-but-real proteins the classifier zeroed — still confound-bounded"
status: candidate
phase: exploratory
kind: discovery
created: 2026-09-15
updated: 2026-09-15

summary: "Boruta all-relevant feature selection (random-forest importance vs shadow/permuted-copy features — a built-in null) on the full 40×2340 log2 experimental matrix confirms 40 proteins as all-relevant to the ADD-vs-Other {PDD,PDCN,HCN} partition (40 Confirmed / 15 Tentative / 2285 Rejected; shadow floor 0.0053; converged 99 iters, seed 20260915). By DE log2FC sign the signature is roughly balanced and bidirectional — 18 up / 22 down in ADD, 0 ambiguous. Of the 40 Confirmed, 33 lie in the elastic-net minimal-optimal set (257 features, finding 0007) and 7 are redundant-but-real — Boruta-confirmed but EN-zeroed (ANO10, ATG3, BAHD1, DOK3, NFX1, UBP4, VP13A) — correlated real signal the L1 penalty dropped; 38/40 are also DE hits (q<0.05, finding 0006), the 2 multivariate-only picks being NFX1 (q=0.054, up) and WDR13 (q=0.080, up). The Confirmed set reads predominantly intracellular/organellar by gene name, but that compartment impression was NOT enrichment-tested here."
verdict: "The all-relevant complement to the classifier's minimal-optimal view: Boruta's built-in shadow null confirms 40 proteins carry real (above-permuted-noise) signal for separating ADD from the pooled others, and it recovers 7 redundant-but-real proteins the elastic net zeroed for parsimony — the honest all-relevant vs minimal-optimal distinction. The signature is balanced and bidirectional (18 up / 22 down in ADD) and 38/40 corroborate the univariate DE hit list. But the ceiling is unchanged: the shadow null addresses relevance-vs-noise, NOT the confound — ADD is 100% Kerchner (finding 0001), so 'all-relevant to ADD' is all-relevant to the ADD-vs-collection-site partition, not necessarily to disease. Exploratory (n=10 positives; the shadow null does not address sampling stability at this n). Not promotable: the compute scripts are scratch (not code-reviewed / promoted). The intracellular/organellar compartment read is a gene-name impression only — a rigorous over-representation test was run separately and found NO selective enrichment (a follow-up finding will record that); it is not asserted here."

entities: []
  # FLAG (see report): the notable named proteins in this finding — the 7 redundant-but-real
  # EN-zeroed (ANO10, ATG3, BAHD1, DOK3, NFX1, UBP4, VP13A), the two multivariate-only picks
  # (NFX1, WDR13), and the secreted exception SERPINA10/ZPI — are not yet normalized to VERIFIED
  # UniProt accessions, and none appear in the sibling findings' entity lists. Left empty rather
  # than fabricating identifiers (conventions/findings.md §2.1: an unverifiable id is a defect).
  # The research-reviewer should normalize + verify the full Confirmed-40 set at promotion.

relationships:
  - { type: relates_to, target: 1, note: "The ceiling: ADD is 100% Kerchner (0001), so Boruta confirming features 'all-relevant to ADD' means all-relevant to the ADD-vs-collection-site partition — not necessarily to disease. The shadow null addresses relevance-vs-noise, NOT the confound." }
  - { type: relates_to, target: 7, note: "All-relevant complement to the classifier's minimal-optimal set (0007): 33/40 Confirmed are in the EN 257-feature set, and 7 are redundant-but-real — Boruta-confirmed but EN-zeroed (ANO10, ATG3, BAHD1, DOK3, NFX1, UBP4, VP13A) — correlated real signal the L1 penalty dropped." }
  - { type: relates_to, target: 6, note: "DE cross-check (0006): 38/40 Confirmed are BH-q<0.05 DE hits; the 2 multivariate-only picks are NFX1 (q=0.054, up) and WDR13 (q=0.080, up); the 18-up/22-down direction split is taken from the DE log2FC signs." }
  - { type: relates_to, target: 3, note: "Class imbalance / small n (0003): n=10 positives. The Boruta shadow null addresses relevance-vs-noise, NOT sampling stability at this n — hence exploratory." }

provenance:
  data_version: "protein-data.csv sha256:90a0462eb03a02bd"
  script: { path: "scripts/scratch/boruta_add_vs_other.py", commit: "uncommitted" }
  params:
    sample_set: "experimental"
    n_experimental: 40
    n_controls_excluded: 0
    n_add: 10
    n_other: 30
    outcome: "ADD-vs-Other"
    positive_label: "ADD"
    negative_levels: "PDD, PDCN, HCN (pooled)"
    n_features: 2340
    method: "boruta_all_relevant"
    estimator: "random_forest_classifier(class_weight=balanced)"
    max_iter: 100
    converged_iters: 99
    alpha: 0.05
    perc: 100
    two_step: true
    shadow_floor: 0.0053
    cv: "none — all-relevant selection runs once on the whole matrix (no CV); feeding the ComBat-corrected matrix is legitimate for selection"
    result_confirmed: 40
    result_tentative: 15
    result_rejected: 2285
    direction_split: "18 up / 22 down in ADD, 0 ambiguous (from DE log2FC sign, finding 0006)"
    selection_scripts:
      - "scripts/scratch/boruta_add_vs_other.py"
      - "scripts/scratch/boruta_selection.py"
    scripts_not_promoted: true
  environment: "uv.lock@uncommitted"
  seeded_from: { template: "boruta", version: "0.1" }
  seed: 20260915
  result_id: "762e16f294e3"

evidence:
  - metric: "Boruta all-relevant decision counts (ADD vs pooled Other)"
    value: 40
    ci: null
    p_value: null
    p_adjusted: null
    correction: "built-in shadow-feature (permuted-copy) null, alpha=0.05, perc=100"
    test: "Boruta all-relevant, random-forest importance vs shadow features; max_iter=100 (converged 99 iters), two-step correction"
    n: 40
    note: "40 Confirmed / 15 Tentative / 2285 Rejected; shadow-noise floor 0.0053; every Confirmed feature's RF importance exceeds the permuted-copy shadow maximum. Runs once on the whole 40×2340 log2 matrix (no CV); feeding the ComBat-corrected matrix is legitimate for selection. Result 762e16f294e3, seed 20260915. Selection, not a per-feature significance test — no per-feature q; minimal-optimal/all-relevant caveats apply."
  - metric: "direction split of the 40 Confirmed (from DE log2FC sign)"
    value: null
    ci: null
    p_value: null
    p_adjusted: null
    correction: null
    test: "sign of the finding-0006 DE log2FC for each Confirmed protein"
    n: 40
    note: "18 up-in-ADD / 22 down-in-ADD / 0 ambiguous — a roughly balanced bidirectional signature."
  - metric: "overlap with the elastic-net minimal-optimal set (finding 0007, 257 non-zero features)"
    value: 33
    ci: null
    p_value: null
    p_adjusted: null
    correction: null
    test: "set overlap, Boruta-Confirmed-40 vs EN 257 non-zero features"
    n: 40
    note: "33/40 Confirmed are in the EN set; the 7 remaining are 'redundant-but-real' — Boruta-confirmed but EN-ZEROED: ANO10, ATG3, BAHD1, DOK3, NFX1, UBP4, VP13A. Correlated real signal the L1 penalty dropped for parsimony — the all-relevant vs minimal-optimal complement."
  - metric: "overlap with the DE hit list (finding 0006, BH q<0.05)"
    value: 38
    ci: null
    p_value: null
    p_adjusted: null
    correction: null
    test: "set overlap, Boruta-Confirmed-40 vs DE BH-q<0.05 hits"
    n: 40
    note: "38/40 Confirmed are also DE hits (q<0.05). The 2 multivariate-only picks are NFX1 (q=0.054, up) and WDR13 (q=0.080, up) — just below the DE threshold but all-relevant to the multivariate partition."

figures:
  - png:        "figures/analysis/boruta/add-vs-other/boruta-importance-add-vs-other.png"
    svg:        "figures/analysis/boruta/add-vs-other/boruta-importance-add-vs-other.svg"
    legend_png: "figures/analysis/boruta/add-vs-other/boruta-importance-add-vs-other.legend.png"
    legend_svg: "figures/analysis/boruta/add-vs-other/boruta-importance-add-vs-other.legend.svg"
    caption:    "Per-iteration random-forest importance box plots for the ADD-vs-Other Boruta run, in Confirmed / Tentative / Rejected blocks with the Confirmed features colored by DE direction (18 up / 22 down in ADD) and the shadow-noise floor (max shadow importance, 0.0053) drawn across. Built-in shadow (permuted-copy) null, no CV; runs once on the whole 40×2340 log2 experimental matrix. EXPLORATORY (n=10 ADD); HARD CONFOUND — ADD=100% Kerchner (finding 0001), so all-relevant to ADD-vs-site not necessarily disease; 7 Confirmed are redundant-but-real (EN-zeroed)."
    script:       { path: "scripts/figures/plot_boruta_importance_add_vs_other.py", commit: "uncommitted" }
    data_version: "protein-data.csv sha256:90a0462eb03a02bd"
    result_id:    "762e16f294e3"
    params:       {}

references: []
  # FLAG (see report): the "predominantly intracellular/organellar" compartment reading in
  # ## Discussion is a gene-name impression, NOT an assertion of enrichment. A rigorous
  # over-representation test was run separately and found NO selective over-representation of
  # that compartment; a follow-up finding will record it. No compartment-enrichment claim is
  # made here, so no reference is asserted.

validation:
  computational_reproduction: { status: not_attempted }
  analytic_replication:       { status: not_attempted }
  data_replication:           { status: not_attempted }

integrity_signoff: true
---

# Boruta all-relevant selection confirms 40 proteins (18 up / 22 down in ADD) for ADD-vs-others — 33 shared with the elastic-net minimal-optimal set plus 7 redundant-but-real proteins the classifier zeroed — still confound-bounded

## Summary
Boruta **all-relevant** feature selection — random-forest importance compared against **shadow (permuted-copy) features**, a built-in null — was run once on the full **40×2340 log2 experimental matrix** for the **ADD vs Other = {PDD, PDCN, HCN}** partition (default params: max_iter=100, alpha=0.05, perc=100; converged at 99 iterations; shadow-noise floor 0.0053; seed 20260915; result `762e16f294e3`). It confirms **40 proteins** as all-relevant (**40 Confirmed / 15 Tentative / 2285 Rejected**). By DE log2FC sign ([finding 0006](0006-add-vs-other-differential-abundance.md)) the signature is roughly balanced and bidirectional: **18 up / 22 down in ADD, 0 ambiguous**. Against the elastic-net **minimal-optimal** set (257 features, [finding 0007](0007-add-vs-other-classification-confound-bounded.md)), **33/40** Confirmed overlap and **7** are **redundant-but-real** — Boruta-confirmed but EN-zeroed (ANO10, ATG3, BAHD1, DOK3, NFX1, UBP4, VP13A) — correlated real signal the L1 penalty dropped for parsimony. **38/40** are also BH-q<0.05 DE hits, the two multivariate-only picks being NFX1 (q=0.054, up) and WDR13 (q=0.080, up).

## Verdict
This is the **all-relevant complement** to the classifier's minimal-optimal picture. Boruta's shadow null confirms that 40 proteins carry real (above-permuted-noise) signal for the ADD-vs-Other split, and it recovers **7 redundant-but-real proteins the elastic net zeroed** — the honest all-relevant vs minimal-optimal distinction. The signature is **balanced and bidirectional** (18 up / 22 down in ADD) and **38/40** corroborate the univariate DE hit list. But the ceiling is unchanged: the shadow null addresses **relevance-vs-noise, not the confound**. ADD is 100% Kerchner ([finding 0001](0001-group-aliased-with-cohort.md)), so "all-relevant to ADD" is all-relevant to the **ADD-vs-collection-site** partition, not necessarily to disease. **Exploratory** — n=10 positives, and the shadow null says nothing about sampling stability at this n ([finding 0003](0003-contrasts-class-imbalanced-small-n.md)). **Not promotable as-is:** the compute scripts are under `scripts/scratch/`, so this finding cannot be validated until they are code-reviewed and moved to `scripts/promoted/`. The Confirmed set reads predominantly intracellular/organellar **by gene name only** — a separate rigorous over-representation test found NO selective enrichment of that compartment (a follow-up finding will record it), so no enrichment claim is made here.

## Evidence
Boruta confirms **40 proteins** as all-relevant to the ADD-vs-Other partition (**40 Confirmed / 15 Tentative / 2285 Rejected**), each Confirmed feature's random-forest importance exceeding the permuted-copy **shadow-noise floor of 0.0053**; the run converged at 99 of 100 iterations (result `762e16f294e3`, seed 20260915). This is a *selection*, not a per-feature significance test — there is no per-feature q. By DE log2FC sign the Confirmed set is **18 up / 22 down in ADD, 0 ambiguous** — a roughly balanced, bidirectional signature. Against the elastic-net minimal-optimal set (257 non-zero features, [finding 0007](0007-add-vs-other-classification-confound-bounded.md)), **33/40** overlap and **7** are **redundant-but-real** — Boruta-confirmed but EN-zeroed (ANO10, ATG3, BAHD1, DOK3, NFX1, UBP4, VP13A). Against the univariate contrast, **38/40** are also BH-q<0.05 DE hits ([finding 0006](0006-add-vs-other-differential-abundance.md)); the two multivariate-only picks are NFX1 (q=0.054, up) and WDR13 (q=0.080, up).

The all-relevant selection is shown as the per-iteration importance distribution against the shadow floor:

![Per-iteration RF importance box plots for ADD-vs-Other Boruta — Confirmed (colored by DE direction, 18 up / 22 down in ADD), Tentative, and Rejected blocks, with the shadow-noise floor (0.0053) drawn across; built-in shadow null, no CV; n=40 (10 ADD).](../figures/analysis/boruta/add-vs-other/boruta-importance-add-vs-other.png)

![Legend for Figure 1 — the Confirmed/Tentative/Rejected blocks, the DE-direction coloring (up vs down in ADD), and the shadow-noise floor line.](../figures/analysis/boruta/add-vs-other/boruta-importance-add-vs-other.legend.png)

*Figure 1. Boruta all-relevant importances for ADD vs pooled Other. Produced by `scripts/figures/plot_boruta_importance_add_vs_other.py` (uncommitted) from data `protein-data.csv sha256:90a0462eb03a02bd`, result `762e16f294e3`.*

Read the blocks left to right: every **Confirmed** box sits to the right of the shadow-noise floor — all 40 cleared the permuted-copy shadows, which is what "all-relevant" means here — while the **Rejected** mass sits at or below it. The Confirmed features split roughly evenly by color into up- and down-in-ADD (18 vs 22), a balanced bidirectional signature rather than a one-directional shift. Because the shadow null only distinguishes real signal from permuted noise, this figure establishes *relevance*, not disease-specificity: the same 40 proteins are all-relevant to the ADD-vs-Kerchner-site partition that ADD is aliased with, so the picture is **confound-bounded**.

## Methods / how to produce
Two scratch scripts drive the selection, both on the delivered analysis matrix `data/protein-data.csv` (`protein-data.csv sha256:90a0462eb03a02bd`, 40 samples × 2340 features), seed 20260915, adapted from the `boruta` template:

- `scripts/scratch/boruta_add_vs_other.py` — binarizes **ADD=1 vs {PDD,PDCN,HCN}=0** and runs Boruta all-relevant selection (random-forest estimator, `class_weight=balanced`; max_iter=100, alpha=0.05, perc=100, two-step correction), caching result `762e16f294e3` (40 Confirmed / 15 Tentative / 2285 Rejected; shadow floor 0.0053).
- `scripts/scratch/boruta_selection.py` — post-processes the cached decision into the Confirmed/Tentative/Rejected sets, cross-references the elastic-net non-zero set (result `902c1e98c018`, [finding 0007](0007-add-vs-other-classification-confound-bounded.md)) and the DE hits ([finding 0006](0006-add-vs-other-differential-abundance.md)) to produce the 33/40 overlap, the 7 redundant-but-real EN-zeroed proteins, the 38/40 DE overlap, and the 18-up/22-down direction split.

Boruta all-relevant selection runs **once on the whole matrix (no CV)** — feeding it the ComBat/limma batch-corrected matrix is legitimate for *selection* (there is no held-out generalization estimate to leak into). Figure 1 is a pure `load_cached_result → plot` job (`scripts/figures/plot_boruta_importance_add_vs_other.py`) over result `762e16f294e3` — it recomputes nothing. The environment is pinned by `uv.lock` (this run added `boruta>=0.4.3`); the stochastic RF/shadow steps are seeded at 20260915. The fingerprint, params, and `run_null?`/status flags are registered in `results/manifest.md`.

**Not promotable as-is:** both compute scripts are under `scripts/scratch/` (not code-reviewed, not promoted). This finding therefore **cannot be validated** until they are reviewed and moved to `scripts/promoted/` (`conventions/findings.md` §2.2, §4).

## Discussion
The value of the all-relevant view is that it names signal the minimal-optimal classifier deliberately discards. The elastic net's L1 penalty keeps a *parsimonious* subset and zeros correlated redundancies; Boruta, by testing every feature against permuted copies, recovers those redundancies as **redundant-but-real** — here the 7 EN-zeroed proteins (ANO10, ATG3, BAHD1, DOK3, NFX1, UBP4, VP13A). That a feature was zeroed by the classifier is therefore **not** evidence it is noise; it may be correlated-but-real, and Boruta is the test that tells the two apart. The strong convergence with the univariate contrast (38/40 also DE hits) and the balanced bidirectional split (18 up / 22 down) both reinforce that this is a coherent multivariate signature, not an artifact of one estimator's inductive bias.

The Confirmed set reads **predominantly intracellular/organellar by gene name** (proteasome, mitochondrial, ER/Golgi, autophagy/endo-lysosomal, nuclear), with the one notable secreted exception being SERPINA10/ZPI — the same qualitative pattern seen in the classifier and DE findings. This is a **gene-name impression, not an enrichment result**: a rigorous over-representation test of that compartment impression was run separately and found **no selective over-representation**, and it will be recorded as its own follow-up finding. No compartment-enrichment claim is asserted here.

## Caveats
- **The ceiling — collection-site confound (relates_to [finding 0001](0001-group-aliased-with-cohort.md)).** ADD is 100% Kerchner, so "all-relevant to ADD" is all-relevant to the ADD-vs-collection-site partition — not necessarily to disease. The Boruta shadow null addresses **relevance-vs-noise**, not the confound; it cannot separate disease biology from site.
- **Exploratory, small-n (relates_to [finding 0003](0003-contrasts-class-imbalanced-small-n.md)).** n=10 positives. The shadow null distinguishes real signal from permuted noise but says nothing about **sampling stability** at this n — a different Confirmed set could emerge on a resample. Never write this up with confirmatory confidence (`conventions/findings.md` §5). The multiplicity context (this is one of several Stage-4 analyses) lives in `findings/exploration-log.md`.
- **Selection, not per-feature significance.** There is no per-feature q. The **minimal-optimal** caveat is complemented, not removed: Boruta answers "all-relevant vs noise", the elastic net answers "minimal sufficient subset" — different questions with different membership.
- **Compartment read is not enrichment.** The predominantly intracellular/organellar impression is by gene name only; a rigorous over-representation test found no selective enrichment (a follow-up finding will record it), so no enrichment is claimed here.
- **Entities not yet normalized.** The notable named proteins (the 7 EN-zeroed, NFX1, WDR13, SERPINA10/ZPI) are not yet mapped to verified UniProt accessions; the research-reviewer should complete and verify the full Confirmed-40 `entities` set at promotion.

## Follow-ups
- **Promote the compute scripts** (`boruta_add_vs_other.py`, `boruta_selection.py` → code review → `scripts/promoted/`) before any validation attempt.
- **Record the compartment over-representation follow-up finding** (the rigorous test that found NO selective enrichment of the intracellular/organellar compartment), and link it here.
- **Normalize and verify the `entities` set** — map the full 40 Confirmed proteins (and the 7 redundant-but-real + NFX1/WDR13/SERPINA10) to verified UniProt accessions (research-reviewer).
- **Pursue blinded analytic replication** once scripts are promoted; and, for any `confirmatory` claim, held-out / orthogonal data replication (none available in-project yet).
- If an ADD cohort from a *different* site ever becomes available, re-run the design that breaks the ADD↔Kerchner alias — the only thing that can convert all-relevance-to-ADD into all-relevance-to-disease.

## Related findings
- **relates_to → [finding 0001](0001-group-aliased-with-cohort.md):** diagnostic group is aliased with source cohort. ADD is 100% Kerchner, so features "all-relevant to ADD" are all-relevant to the ADD-vs-collection-site partition — the shadow null cannot separate disease from site. The hard ceiling.
- **relates_to → [finding 0007](0007-add-vs-other-classification-confound-bounded.md):** the elastic-net / SVM classifier. This finding is its all-relevant complement — 33/40 Confirmed lie in the EN minimal-optimal 257-feature set, and 7 (ANO10, ATG3, BAHD1, DOK3, NFX1, UBP4, VP13A) are redundant-but-real correlated signal the L1 penalty zeroed.
- **relates_to → [finding 0006](0006-add-vs-other-differential-abundance.md):** the univariate ADD-vs-Other DE result. 38/40 Confirmed are BH-q<0.05 DE hits; the 2 multivariate-only picks are NFX1 (q=0.054, up) and WDR13 (q=0.080, up); the 18-up/22-down direction split is read from the DE log2FC signs.
- **relates_to → [finding 0003](0003-contrasts-class-imbalanced-small-n.md):** class imbalance / small n (10 ADD positives). The shadow null addresses relevance-vs-noise, not sampling stability at this n — which is why this finding stays exploratory.

## References
None. The "predominantly intracellular/organellar" reading is a gene-name impression, not an enrichment assertion — a rigorous over-representation test found no selective over-representation (to be recorded as a follow-up finding), so no background/interpretive claim requiring a reference is made here.
</content>
</invoke>
