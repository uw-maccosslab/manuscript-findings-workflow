---
id: 7
title: "The plasma-EV proteome predicts ADD-vs-others above chance (AUC ~0.96, beats label-shuffle null p=0.001) but the signal is confound-bounded — inseparable from collection site — and rides on intracellular proteins mostly reduced in ADD"
status: candidate
phase: exploratory
kind: discovery
created: 2026-09-15
updated: 2026-09-15

summary: "Two leakage-safe linear classifiers (elastic net, linear SVM) separate ADD (n=10) from pooled Other {PDD,PDCN,HCN} (n=30) on the delivered 40×2340 log2 plasma-EV matrix at nested-CV ROC-AUC 0.940 (EN) / 0.963 (SVM), and BOTH beat a 1000-permutation label-shuffle null (empirical p=0.001), ruling out leakage and pure overfitting. But ADD is 100% Kerchner cohort, so above-null discrimination cannot be attributed to disease vs collection site; a cohort-site probe among the non-ADD samples collapses to null (AUC 0.34, p=0.983), which bounds rather than clears the confound because ADD's batch correction was extrapolated, not estimated. The discriminative basis is predominantly intracellular/organellar proteins mostly reduced in ADD, not secreted EV markers, and 23/25 of the top EN features are also DE hits with 100% direction agreement."
verdict: "There is a genuine, reproducible, null-beating proteomic signal distinguishing ADD from the other groups — leakage and overfitting are ruled out. But this dataset cannot determine whether that signal is ADD biology or the Kerchner collection site: ADD and Kerchner are the same partition among these 40 samples, and the one test that could adjudicate it (the site probe) is blocked because ADD's on-cohort batch correction is extrapolated, not estimated. The site probe collapsing to null bounds residual site leakage in the non-ADD subset; it does not exonerate the ADD result. Coefficients are now licensed (null passed) and magnitudes trustworthy (converged), and they point at intracellular/organellar proteins mostly lower in ADD — not EV-surface cargo. This is hypothesis-generating only, and confound-bounded is the ceiling on any disease interpretation. Not promotable: the compute scripts are scratch (not reviewed/promoted), so this finding cannot be validated as-is."

entities:
  - { db: uniprot, id: "Q16822", label: "PCK2 (PEPCK-M, mitochondrial)" }
  - { db: uniprot, id: "P29350", label: "PTN6/SHP1" }
  - { db: uniprot, id: "O60763", label: "USO1" }
  - { db: uniprot, id: "Q9P253", label: "VPS18" }
  - { db: uniprot, id: "P26440", label: "IVD" }
  - { db: uniprot, id: "P19404", label: "NDUV2/NDUFV2" }
  - { db: uniprot, id: "P36507", label: "MP2K2/MEK2" }
  - { db: uniprot, id: "P40939", label: "HADHA" }
  - { db: uniprot, id: "P62333", label: "PSMC6" }
  - { db: uniprot, id: "O43242", label: "PSMD3" }
  - { db: uniprot, id: "P43686", label: "PSMC4" }
  - { db: uniprot, id: "O14653", label: "GOSR2" }
  - { db: uniprot, id: "Q9HC98", label: "NEK6" }
  - { db: uniprot, id: "Q8TDX7", label: "NEK7" }

relationships:
  - { type: relates_to, target: 1, note: "The hard ceiling: ADD is 100% Kerchner (0001), so above-null discrimination is inseparable from collection site; the label-shuffle null permutes that same partition and cannot resolve it." }
  - { type: relates_to, target: 4, note: "Why the site probe cannot clear the ADD result: ADD is single-cohort, so its limma on-cohort batch correction is extrapolated not estimated (0004) — a site-null in the non-ADD subset does not certify the ADD stratum." }
  - { type: relates_to, target: 3, note: "Class imbalance / small n: ADD 10 vs Other 30, and the site probe is n=27 (7 vs 20) — a null-level site AUC bounds rather than proves absence of residual site leakage (0003)." }
  - { type: relates_to, target: 6, note: "DE cross-check: 23/25 of the EN top-25 features are also BH-q<0.05 DE hits (0006), with 100% coefficient-sign vs DE-log2FC direction agreement — the classifier and the univariate contrast converge on the same intracellular, mostly-down-in-ADD basis." }

provenance:
  data_version: "protein-data.csv sha256:90a0462eb03a02bd"
  script: { path: "scripts/scratch/classify_add_vs_other.py", commit: "uncommitted" }
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
    generalization_target: "unseen individuals, same cohort structure"
    cv: "stratified 5-fold × 5-repeat nested CV, row-level (individual)"
    models:
      elastic_net:
        result_id: "902c1e98c018"
        script: "scripts/scratch/classify_add_vs_other.py"
        run_null: true
        max_iter: 20000
        converged: true
      linear_svm:
        result_id: "a0e23fb510a0"
        script: "scripts/scratch/classify_add_vs_other.py"
        run_null: true
        C: 0.001
        n_support_vectors: 22
        tuning_warning: "residual TuningNoiseWarning — 13/25 folds tune below the C plateau (AUC unchanged by grid fix); flag for stats-reviewer"
    null_model:
      type: "label-shuffle (unit-level)"
      script: "scripts/scratch/classify_add_vs_other_null.py"
      n_permutations: 1000
      null_repeats: 3
      result: "both models beat null, empirical permutation p=0.001"
    cohort_site_probe:
      script: "scripts/scratch/classify_cohort_site_probe.py"
      result_id: "b31f645d39ee"
      sample_subset: "non_ADD {PDD,PDCN,HCN}"
      outcome: "Kerchner (n=7) vs Poston (n=20)"
      run_null: true
      result: "AUC 0.34 (observed 0.475), null mean 0.500, permutation p=0.983 — collapses to null; 1 feature selected"
    batch_correct: { method: "limma::removeBatchEffect", batch_column: "Cohort", covariate_preserved: "Condition", applied: "upstream on delivered matrix", note: "ADD single-cohort → correction extrapolated not estimated; see finding 0004" }
    scripts_not_promoted: true
  environment: "uv.lock@uncommitted"
  seeded_from: { template: "classification", version: "0.1" }
  seed: 20260915
  result_id: "902c1e98c018"

evidence:
  - metric: "nested CV ROC-AUC (elastic net, ADD vs pooled Other)"
    value: 0.940
    ci: null
    p_value: 0.001
    p_adjusted: null
    correction: "empirical permutation (1000 label shuffles, unit-level)"
    test: "elastic-net logistic, stratified 5-fold×5-repeat nested CV, generalization target = unseen individuals (same cohort structure)"
    n: 40
    note: "Fold SD 0.09; observed fixed-hp AUC 0.956; balanced acc 0.83, PR-AUC 0.87, MCC 0.63; 257 non-zero features. Result 902c1e98c018 (run_null=True, converged max_iter=20000). Above-null p=0.001 ⇒ NOT leakage/overfitting."
  - metric: "nested CV ROC-AUC (linear SVM, ADD vs pooled Other)"
    value: 0.963
    ci: null
    p_value: 0.001
    p_adjusted: null
    correction: "empirical permutation (1000 label shuffles, unit-level)"
    test: "linear SVM, stratified 5-fold×5-repeat nested CV, generalization target = unseen individuals (same cohort structure)"
    n: 40
    note: "Observed AUC 0.972; balanced acc 0.74, PR-AUC 0.92, MCC 0.52; C=0.001, 22 support vectors. Result a0e23fb510a0 (run_null=True). Above-null p=0.001. NOTE residual TuningNoiseWarning (13/25 folds tune below C plateau) — AUC unchanged by grid fix; flag for stats-reviewer."
  - metric: "label-shuffle null mean ROC-AUC (EN)"
    value: 0.502
    ci: [0.244, 0.778]
    p_value: null
    p_adjusted: null
    correction: null
    test: "1000-permutation unit-level label-shuffle null, elastic net"
    n: 40
    note: "Null centered at chance; observed 0.940 far to the right → leakage and pure overfitting ruled out."
  - metric: "label-shuffle null mean ROC-AUC (SVM)"
    value: 0.500
    ci: [0.228, 0.756]
    p_value: null
    p_adjusted: null
    correction: null
    test: "1000-permutation unit-level label-shuffle null, linear SVM"
    n: 40
    note: "Null centered at chance; observed 0.963 far to the right."
  - metric: "EN vs SVM mean paired per-fold AUC difference"
    value: 0.023
    ci: null
    p_value: 0.068
    p_adjusted: null
    correction: null
    test: "Wilcoxon signed-rank over paired per-fold AUCs (25 folds)"
    n: 40
    note: "A TIE within seed noise: 21/25 folds identical; no evidence either model dominates."
  - metric: "cohort-site probe ROC-AUC (Kerchner vs Poston among non-ADD)"
    value: 0.34
    ci: null
    p_value: 0.983
    p_adjusted: null
    correction: "empirical permutation (1000 label shuffles, unit-level)"
    test: "elastic-net logistic, class_weight=balanced, stratified 5-fold×5-repeat nested CV; Kerchner(7) vs Poston(20)"
    n: 27
    note: "Observed 0.475, null mean 0.500; COLLAPSES to null (p=0.983); only 1 feature selected. Result b31f645d39ee. On-cohort batch correction removed detectable within-non-ADD site structure — but does NOT clear ADD (single-cohort → extrapolated correction, 0004); small/imbalanced (7v20) so this BOUNDS rather than proves absence of residual site leakage."
  - metric: "EN coefficient stability (top 25 selected proteins)"
    value: null
    ci: null
    p_value: null
    p_adjusted: null
    correction: null
    test: "selection frequency + sign consistency over stability resamples, elastic net (converged)"
    n: 40
    note: "Top 25 all have selection frequency ≥0.88 and sign consistency 1.00; 16 negative (down in ADD) + 9 positive (up in ADD). Top: PCK2(−), SHP1(−), CDC42SE2(+), UBN2(+), USO1(−), NEK6/7(−), MEK2(−), VPS18(+), IVD(+), NDUFV2(−). Selection, not a per-feature significance test — no per-feature q; minimal-optimal caveat applies. Coefficients now LICENSED (null passed) and magnitudes trustworthy (converged)."
  - metric: "classifier↔DE overlap (EN top-25 vs BH-q<0.05 DE hits, finding 0006)"
    value: 23
    ci: null
    p_value: null
    p_adjusted: null
    correction: null
    test: "set overlap + coefficient-sign vs DE-log2FC direction agreement"
    n: 25
    note: "23/25 of the EN top-25 are also DE hits (BH-q<0.05), 100% coefficient-sign vs DE-log2FC direction agreement. The classifier and the univariate contrast (0006) converge on the same basis."

figures:
  - png:        "figures/analysis/classification/add-vs-other/roc-add-vs-other-nullpass.png"
    svg:        "figures/analysis/classification/add-vs-other/roc-add-vs-other-nullpass.svg"
    caption:    "Nested-CV ROC (TPR vs FPR) for the primary ADD-vs-Other classifiers, from the current null-pass caches. Left: elastic-net logistic (scores = predicted probability); right: linear SVM (scores = uncalibrated decision-function margin). 25 per-fold ROC curves (thin) behind the mean ROC (bold) with ±1 SD band and chance diagonal. Nested ROC-AUC 0.940 ± 0.09 (EN), 0.963 ± 0.07 (SVM); per-fold AUC range 0.75–1.00; n=40 (10 ADD / 30 Other). The label-shuffle null passed (p=0.001) → not leakage/overfitting, but exploratory (small n) and confound-bounded: ADD=100% Kerchner (finding 0001), so above-null discrimination cannot be separated from collection site. ROC key on-axes."
    script:       { path: "scripts/figures/plot_roc_add_vs_other.py", commit: "uncommitted" }
    data_version: "protein-data.csv sha256:90a0462eb03a02bd"
    result_id:    "902c1e98c018"
    params:       { svm_result_id: "a0e23fb510a0", note: "EN result 902c1e98c018 (primary) + SVM result a0e23fb510a0; on-axes ROC key, no separate legend image (ROC exception)" }
  - png:        "figures/analysis/classification/add-vs-other/roc-vs-null-add-vs-other.png"
    svg:        "figures/analysis/classification/add-vs-other/roc-vs-null-add-vs-other.svg"
    legend_png: "figures/analysis/classification/add-vs-other/roc-vs-null-add-vs-other.legend.png"
    legend_svg: "figures/analysis/classification/add-vs-other/roc-vs-null-add-vs-other.legend.svg"
    caption:    "Label-shuffle null vs observed ROC-AUC, ADD vs pooled Other, both models: the 1000-permutation unit-level null distribution (centered ~0.5) against the observed nested-CV AUC (EN 0.940, SVM 0.963), empirical p=0.001 for each. Both beat the null → not leakage/overfitting; but ADD=100% Kerchner, so above-null ≠ disease vs collection site. Exploratory, confound-bounded."
    script:       { path: "scripts/figures/plot_roc_vs_null_add_vs_other.py", commit: "uncommitted" }
    data_version: "protein-data.csv sha256:90a0462eb03a02bd"
    result_id:    "902c1e98c018"
    params:       { second_result_id: "a0e23fb510a0", note: "EN result 902c1e98c018 + SVM result a0e23fb510a0 overlaid" }
  - png:        "figures/analysis/classification/add-vs-other/en-coefficients-add-vs-other-nullpass.png"
    svg:        "figures/analysis/classification/add-vs-other/en-coefficients-add-vs-other-nullpass.svg"
    legend_png: "figures/analysis/classification/add-vs-other/en-coefficients-add-vs-other-nullpass.legend.png"
    legend_svg: "figures/analysis/classification/add-vs-other/en-coefficients-add-vs-other-nullpass.legend.svg"
    caption:    "Top 25 of 257 non-zero elastic-net coefficients (signed, standardized), ADD vs pooled Other: 16 down / 9 up in ADD, all selection frequency ≥0.88 and sign consistency 1.00. Null passed → coefficients licensed; converged (max_iter=20000) → magnitudes trustworthy. Minimal-optimal (Boruta not run); confound-bounded; 23/25 also DE hits (finding 0006)."
    script:       { path: "scripts/figures/plot_en_coefficients_add_vs_other.py", commit: "uncommitted" }
    data_version: "protein-data.csv sha256:90a0462eb03a02bd"
    result_id:    "902c1e98c018"
    params:       {}
  - png:        "figures/analysis/classification/cohort-site-probe/site-probe-null.png"
    svg:        "figures/analysis/classification/cohort-site-probe/site-probe-null.svg"
    legend_png: "figures/analysis/classification/cohort-site-probe/site-probe-null.legend.png"
    legend_svg: "figures/analysis/classification/cohort-site-probe/site-probe-null.legend.svg"
    caption:    "Cohort-site probe vs the ADD signal: predicting Kerchner (n=7) vs Poston (n=20) among the 30 non-ADD samples collapses to null (AUC 0.34, observed 0.475, p=0.983) while ADD-vs-Other beats its null (p=0.001). The collapse does NOT clear ADD — ADD is single-cohort, so its batch correction is extrapolated (finding 0004) — and the probe is small/imbalanced (7 vs 20), so it BOUNDS rather than eliminates residual site leakage."
    script:       { path: "scripts/figures/plot_site_probe_null.py", commit: "uncommitted" }
    data_version: "protein-data.csv sha256:90a0462eb03a02bd"
    result_id:    "b31f645d39ee"
    params:       { reference_result_id: "902c1e98c018", note: "site probe b31f645d39ee shown against the ADD EN result 902c1e98c018" }

references: []
  # FLAG (see report): the intracellular/organellar vs secreted-EV-marker reading in
  # ## Discussion is a background/interpretive claim (conventions/findings.md §2.5) and a
  # rigorous GO-CC over-representation test has NOT been run. Add references + run the
  # enrichment test before promotion; left empty rather than fabricating identifiers.

validation:
  computational_reproduction: { status: not_attempted }
  analytic_replication:       { status: not_attempted }
  data_replication:           { status: not_attempted }

integrity_signoff: true
---

# The plasma-EV proteome predicts ADD-vs-others above chance (AUC ~0.96, beats label-shuffle null p=0.001) but the signal is confound-bounded — inseparable from collection site — and rides on intracellular proteins mostly reduced in ADD

## Summary
Two leakage-safe linear classifiers — an **elastic net** and a **linear SVM** — were trained to separate **ADD (n=10)** from a pooled **Other = {PDD, PDCN, HCN} (n=30)** group on the delivered 40×2340 log2 plasma-EV matrix, using stratified 5-fold × 5-repeat **nested** CV at the individual level (generalization target = unseen individuals, same cohort structure; seed 20260915). Both discriminate ADD well — nested ROC-AUC **0.940** (EN) and **0.963** (SVM) — and both **beat a 1000-permutation label-shuffle null** at empirical permutation **p=0.001**, which rules out leakage and pure overfitting. The catch is the design: **ADD is 100% Kerchner cohort** ([finding 0001](0001-group-aliased-with-cohort.md)), so among these 40 samples "ADD" and "Kerchner" are the same partition, and above-null discrimination **cannot** be attributed to disease vs collection site. The discriminative basis is predominantly **intracellular/organellar proteins mostly reduced in ADD** — not secreted EV-surface markers — and **23/25** of the top elastic-net features are also BH-q<0.05 DE hits ([finding 0006](0006-add-vs-other-differential-abundance.md)) with **100% direction agreement**.

## Verdict
There is a **genuine, reproducible, null-beating proteomic signal** distinguishing ADD from the other groups: leakage and overfitting are ruled out by the label-shuffle null. But this dataset **cannot determine whether the signal is ADD biology or the Kerchner collection site** — ADD and Kerchner are the same partition here, and the label-shuffle null does **not** resolve it because it permutes that same partition (it certifies *not-leakage / not-overfitting* only). The one test that could adjudicate — the cohort-site probe — is blocked: it collapses to null among the non-ADD samples, but because ADD is single-cohort its batch correction was **extrapolated, not estimated** ([finding 0004](0004-batch-correction-extrapolated-single-cohort.md)), so the probe **bounds** residual site leakage rather than clearing the ADD result. Coefficients are now **licensed** (null passed) and their **magnitudes trustworthy** (the elastic net converged), and they point at intracellular/organellar proteins mostly lower in ADD. This is **hypothesis-generating only**; confound-bounded is the ceiling on any disease interpretation. **Not promotable as-is:** the compute scripts are under `scripts/scratch/`, so this finding cannot be validated until they are code-reviewed and moved to `scripts/promoted/`.

## Evidence
Both classifiers separate ADD from pooled Other well and, critically, **beat a label-shuffle null**. The elastic net reaches nested ROC-AUC **0.940** (fold SD 0.09; observed fixed-hp 0.956; balanced acc 0.83, PR-AUC 0.87, MCC 0.63; 257 non-zero features; result `902c1e98c018`, converged at max_iter=20000). The linear SVM reaches **0.963** (observed 0.972; balanced acc 0.74, PR-AUC 0.92, MCC 0.52; C=0.001, 22 support vectors; result `a0e23fb510a0`). Under the 1000-permutation unit-level **label-shuffle null**, the EN null mean is 0.502 [0.244, 0.778] and the SVM null mean 0.500 [0.228, 0.756]; **both observed AUCs sit at empirical permutation p=0.001**. Leakage and pure overfitting are therefore ruled out — a real above-null signal separates ADD from Other. EN vs SVM is a **tie** (mean paired per-fold AUC difference 0.023, within seed noise; 21/25 folds identical; Wilcoxon p=0.068).

The standard performance view of this discrimination is the nested-CV ROC: both classifiers bow sharply toward the top-left, and the per-fold spread is the honest small-n uncertainty.

![Nested-CV ROC (TPR vs FPR), ADD vs pooled Other, both models — 25 per-fold curves (thin) behind the mean ROC (bold) with ±1 SD band and chance diagonal; nested AUC 0.940±0.09 (EN, left) / 0.963±0.07 (SVM, right); n=40.](../figures/analysis/classification/add-vs-other/roc-add-vs-other-nullpass.png)

*Figure 1. Nested-CV ROC for the primary ADD-vs-Other classifiers (EN + SVM), from the current null-pass caches; ROC key on-axes (no separate legend image — the documented ROC exception). Produced by `scripts/figures/plot_roc_add_vs_other.py` (uncommitted) from data `protein-data.csv sha256:90a0462eb03a02bd`, results `902c1e98c018` (EN) + `a0e23fb510a0` (SVM).*

Both classifiers bow sharply toward the top-left — EN nested AUC 0.940±0.09, SVM 0.963±0.07 — while the per-fold curves fan out to AUC 0.75–1.00, the honest small-n uncertainty (10 ADD / 30 Other). This is the standard performance view of the discrimination whose reality is established by the permutation null (next figure) and whose interpretation is bounded by the cohort confound.

![Label-shuffle null vs observed ROC-AUC, ADD vs pooled Other — the permutation null distributions (centered ~0.5) against the observed nested-CV AUC (EN 0.940, SVM 0.963), empirical p=0.001; n=40.](../figures/analysis/classification/add-vs-other/roc-vs-null-add-vs-other.png)

![Legend for Figure 2 — the null distribution, the observed-AUC markers for EN and SVM, and the empirical-p annotation.](../figures/analysis/classification/add-vs-other/roc-vs-null-add-vs-other.legend.png)

*Figure 2. Label-shuffle null vs observed AUC, both models. Produced by `scripts/figures/plot_roc_vs_null_add_vs_other.py` (uncommitted) from data `protein-data.csv sha256:90a0462eb03a02bd`, results `902c1e98c018` (EN) + `a0e23fb510a0` (SVM).*

The null mass piles up at 0.5 while both observed AUCs land far to the right at ~0.96 with p=0.001 — the picture of a signal that is real (not leakage, not overfitting), but which, because ADD=100% Kerchner, is above-null with respect to *some* ADD/Kerchner-aligned axis the null cannot tease apart from disease.

With the null passed, the elastic-net coefficients are **licensed** and, because the model converged, their **magnitudes are trustworthy**. The top 25 of 257 non-zero features are highly stable (selection frequency ≥0.88, sign consistency 1.00) and lean negative: **16 down / 9 up in ADD**. The leading features are PCK2(−), SHP1(−), CDC42SE2(+), UBN2(+), USO1(−), NEK6/7(−), MEK2(−), VPS18(+), IVD(+), NDUFV2(−).

![Top 25 of 257 signed elastic-net coefficients, ADD vs pooled Other — 16 down / 9 up in ADD, selection frequency ≥0.88, sign consistency 1.00; null passed (licensed), converged (magnitudes trustworthy).](../figures/analysis/classification/add-vs-other/en-coefficients-add-vs-other-nullpass.png)

![Legend for Figure 3 — coefficient sign (down vs up in ADD) and the stability annotation (selection frequency, sign consistency).](../figures/analysis/classification/add-vs-other/en-coefficients-add-vs-other-nullpass.legend.png)

*Figure 3. Top-25 elastic-net coefficients (null-pass). Produced by `scripts/figures/plot_en_coefficients_add_vs_other.py` (uncommitted) from data `protein-data.csv sha256:90a0462eb03a02bd`, result `902c1e98c018`.*

Read the bars: the discriminative basis is **predominantly proteins lower in ADD**, and by identity they are intracellular/organellar (mitochondrial, proteasome, kinases, vesicle-trafficking) rather than EV-surface cargo. This is a *selection*, not a per-feature significance test — there is no per-feature q, and the minimal-optimal caveat applies (a lower-ranked protein is not thereby unimportant; correlated features share credit). The convergence with the univariate contrast is strong: **23/25** of these top features are also BH-q<0.05 DE hits ([finding 0006](0006-add-vs-other-differential-abundance.md)), with **100%** coefficient-sign vs DE-log2FC direction agreement.

The cohort-site probe is the test that would adjudicate disease vs site — and its result must be read in **both directions, without overclaim**. Predicting **Kerchner (n=7) vs Poston (n=20)** among the 30 non-ADD samples **collapses to null**: ROC-AUC 0.34 (observed 0.475), null mean 0.500, permutation **p=0.983**, only 1 feature selected (result `b31f645d39ee`).

![Cohort-site probe vs the ADD signal — the site probe (Kerchner vs Poston among non-ADD) collapses to null (AUC 0.34, p=0.983) while ADD-vs-Other beats its null (p=0.001); n=27 for the probe (7 vs 20), n=40 for ADD.](../figures/analysis/classification/cohort-site-probe/site-probe-null.png)

![Legend for Figure 4 — the site-probe null and observed AUC against the ADD-vs-Other observed AUC.](../figures/analysis/classification/cohort-site-probe/site-probe-null.legend.png)

*Figure 4. Site probe vs the ADD signal. Produced by `scripts/figures/plot_site_probe_null.py` (uncommitted) from data `protein-data.csv sha256:90a0462eb03a02bd`, result `b31f645d39ee` (+ `902c1e98c018` as the ADD reference).*

The probe sits squarely in its null while the ADD signal sits far right of its own null. Read one way, the on-cohort limma `removeBatchEffect` **removed detectable within-non-ADD site structure** — site is not predictable among the non-ADD samples. Read honestly, this does **not** clear the ADD result: ADD is single-cohort, so its correction was **extrapolated not estimated** ([finding 0004](0004-batch-correction-extrapolated-single-cohort.md)), and the probe is small and imbalanced (n=27, 7 vs 20), so a null-level AUC **bounds** residual site leakage rather than proving none. The site probe cannot exonerate a confound in a stratum it does not contain.

## Methods / how to produce
Three scratch compute scripts drive this finding, all on the delivered analysis matrix `data/protein-data.csv` (`protein-data.csv sha256:90a0462eb03a02bd`, 40 samples × 2340 features), seed 20260915, from project copies of the `classification` / `classification-svm` templates:

- `scripts/scratch/classify_add_vs_other.py` — builds the binary factor **ADD vs Other = {PDD,PDCN,HCN} pooled** and runs the elastic-net and linear-SVM classifiers under stratified 5-fold×5-repeat nested CV (generalization target = unseen individuals), caching results `902c1e98c018` (EN) and `a0e23fb510a0` (SVM).
- `scripts/scratch/classify_add_vs_other_null.py` — runs the 1000-permutation unit-level **label-shuffle null** (null_repeats=3) for both models; both observed AUCs land at empirical p=0.001. This is the leakage check the univariate DE cannot provide.
- `scripts/scratch/classify_cohort_site_probe.py` — the **cohort-site probe**: restricts to the 30 non-ADD samples and predicts **Kerchner vs Poston** (class_weight=balanced), run against its own label-shuffle null; caches result `b31f645d39ee`.

Figures 1–3 are pure `load_cached_result → plot` jobs (`scripts/figures/plot_roc_vs_null_add_vs_other.py`, `plot_en_coefficients_add_vs_other.py`, `plot_site_probe_null.py`) over the named cached results — they recompute nothing. The environment is pinned by `uv.lock`; stochastic steps are seeded at 20260915. Result fingerprints, params, and the `run_null?`/status flags are registered in `results/manifest.md`.

**Not promotable as-is:** all three compute scripts are under `scripts/scratch/` (not code-reviewed, not promoted). This finding therefore **cannot be validated** until they are reviewed and moved to `scripts/promoted/` (`conventions/findings.md` §2.2, §4). The SVM result carries a residual **TuningNoiseWarning** (13/25 folds tune below the C plateau; AUC unchanged by the grid fix) that the stats-reviewer should assess.

## Discussion
The honest reading is a two-part one. First, the signal is **real in the not-leakage / not-overfitting sense** — the label-shuffle null settles that decisively, and it is what licenses the coefficients and makes their converged magnitudes trustworthy. Second, the signal is **confound-bounded**: ADD and the Kerchner collection site are one and the same partition among these 40 samples, so nothing in this dataset can tell disease biology apart from site. The cohort-site probe is the natural adjudicator and it collapses to null, but the collapse is *bounded evidence* — ADD's single-cohort correction is extrapolated, and the probe is small — so it lowers, without eliminating, the plausibility of residual site leakage. The discriminative proteins compound the caution: they are predominantly **intracellular/organellar** (mitochondrial PCK2/NDUFV2/IVD/HADHA, proteasome PSMC6/PSMD3/PSMC4, vesicle-trafficking USO1/VPS18/GOSR2, signalling MEK2/SHP1) and **mostly down in ADD** — not the secreted plasma/EV-surface markers a clean EV-biology story would predict. This is a **gene-name identity read**: a rigorous GO-CC over-representation test has **not** been run, and doing so is the honest next step (the intracellular/organellar interpretation here is a background claim that requires `references` and the enrichment test before promotion — see Caveats and the report flag).

## Caveats
- **Exploratory, small-n, imbalanced (relates_to [finding 0003](0003-contrasts-class-imbalanced-small-n.md)).** ADD 10 vs Other 30; the models were derived on the same data they are read on. The multiplicity context — this is the primary of several Stage-4 analyses explored — lives in `findings/exploration-log.md`. Never write this up with confirmatory confidence (`conventions/findings.md` §5).
- **The central limitation — collection-site confound (relates_to [finding 0001](0001-group-aliased-with-cohort.md)).** ADD is 100% Kerchner, so above-null discrimination **cannot** be attributed to disease vs site. The label-shuffle null does **not** resolve this — it permutes that same partition, so it establishes *not-leakage / not-overfitting only*. This is the ceiling on any disease interpretation.
- **The site probe bounds, it does not clear (relates_to [finding 0004](0004-batch-correction-extrapolated-single-cohort.md)).** The site-null in the non-ADD subset shows the on-cohort correction removed detectable within-non-ADD site structure, but ADD is single-cohort so its correction is **extrapolated not estimated**, and the probe is small/imbalanced (n=27, 7 vs 20). A null-level site AUC therefore bounds residual site leakage rather than proving none.
- **Selection, not per-feature significance.** The coefficient ranking is a stability read (selection frequency, sign consistency), not a per-feature test — there is no per-feature q. The **minimal-optimal** caveat applies: a low-ranked or unselected protein is not thereby unimportant (correlated features), and Boruta (all-relevant) was not run.
- **SVM tuning warning.** Result `a0e23fb510a0` carries a residual TuningNoiseWarning (13/25 folds tune below the C plateau); AUC is unchanged by the grid fix, but the stats-reviewer should confirm.
- **Intracellular-vs-EV reading is not yet enrichment-tested.** The organelle/compartment identity read rests on gene names; a GO-CC over-representation test has not been run, and the interpretive claim needs `references`.
- **Entity ids** cover only the named notable proteins; several top features (e.g. CDC42SE2, UBN2, JMJD6) are not yet normalized to verified accessions — the research-reviewer should complete and confirm the `entities` set at promotion.

## Follow-ups
- **Promote the compute scripts** (`classify_add_vs_other.py`, `classify_add_vs_other_null.py`, `classify_cohort_site_probe.py` → code review → `scripts/promoted/`) before any validation attempt.
- **Have the stats-reviewer assess the SVM TuningNoiseWarning** on result `a0e23fb510a0`.
- **Run a GO-CC (cellular-component) over-representation test** on the discriminative features to test the intracellular/organellar reading rigorously, and add fact-checked `references`.
- **Pursue blinded analytic replication** and, for any `confirmatory` claim, held-out / orthogonal data replication (none available in-project yet).
- **Complete and verify the `entities` set** (normalize the full top-25 to UniProt accessions; research-reviewer).
- If an ADD cohort from a *different* site ever becomes available, re-run the design that would actually break the ADD↔Kerchner alias — the only thing that can convert this from confound-bounded to disease-attributable.

## Related findings
- **relates_to → [finding 0001](0001-group-aliased-with-cohort.md):** diagnostic group is aliased with source cohort. ADD is 100% Kerchner, so above-null discrimination is inseparable from a collection-site effect — the hard ceiling on this finding.
- **relates_to → [finding 0004](0004-batch-correction-extrapolated-single-cohort.md):** the matrix is batch-corrected on cohort, and ADD's single-cohort correction is extrapolated not estimated — which is exactly why the site probe collapsing to null cannot clear the ADD result.
- **relates_to → [finding 0003](0003-contrasts-class-imbalanced-small-n.md):** class imbalance / small n (ADD 10 vs 30; site probe 7 vs 20), which is why the site-null bounds rather than proves absence of residual site leakage.
- **relates_to → [finding 0006](0006-add-vs-other-differential-abundance.md):** the univariate ADD-vs-Other DE result. 23/25 of the EN top-25 are also BH-q<0.05 DE hits with 100% direction agreement — the classifier and the contrast converge on the same intracellular, mostly-down-in-ADD basis.

## References
None yet — the intracellular/organellar-vs-secreted-EV-marker reading in the Discussion is a background/interpretive claim that requires `references` added and fact-checked by the research-reviewer, and a GO-CC enrichment test run, before promotion (see Caveats and Follow-ups).
</content>
