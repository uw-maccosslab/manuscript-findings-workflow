# QC processing states

This project's proteomics data has **only one delivered processing state**. There is
no raw -> normalized -> batch-corrected series to materialize here, because the
pre-normalization and pre-batch-correction matrices were never delivered to this
project (see `state/DATA_DESCRIPTION.md` and finding
[0004](../../findings/0004-batch-correction-extrapolated-single-cohort.md)).

## What "delivered" means

`data/protein-data.csv` arrives:

- **log2**-scale (not linear);
- **already median-normalized** upstream (per-sample medians span 0.15 log2 units);
- **already batch-corrected on `Cohort`** upstream, via
  `limma::removeBatchEffect(..., design = ~Condition, batch = Cohort)` (protecting
  `Condition`) -- see finding 0004. The pre-correction matrix (`nobatchProt` in the
  scientist's R pipeline) is **not available** in this project.

Consequently `scripts/scratch/qc_prep.py` does **not** call
`common.normalize.normalize()` / `log2_transform()` (the data is already log2;
those functions require a **linear**-scale input and would refuse to run -- this is a
correctness guard, not an oversight) and does **not** run ComBat /
`removeBatchEffect` again (the matrix is already corrected; re-running it would
double-correct).

## The two states materialized here

- **`delivered_log/`** (`scale="log2"`) -- the Stage-4 analysis-subset `Dataset` (40
  samples x 2340 proteins) exactly as loaded, no transform applied. Consumed by: PCA,
  sample-correlation, abundance-boxplot.
- **`delivered_linear/`** (`scale="linear"`) --
  `common.normalize.to_linear(delivered_log)`, the mathematical de-log (`2**x`) of
  the delivered matrix. Consumed by: CV, dynamic-range.

**`delivered_linear` is NOT raw / unnormalized intensity.** It is the de-logged image
of the already-normalized, already-batch-corrected delivered matrix -- a linear-scale
view materialized only because coefficient-of-variation and dynamic-range are
undefined on a log scale (a per-feature mean near 0 makes `std/mean` explode, and
"detected" is a linear-scale concept). Do not read or label it as raw intensity.

## Controls are absent

The 8 pooled-reference control samples (`CP01`-`CP04`, `KP01`-`KP04`) exist only in
`data/metadata.csv`; they have no columns in `data/protein-data.csv` (finding
[0003](../../findings/0003-contrasts-class-imbalanced-small-n.md)). Both states above
therefore contain **only the 40 experimental samples** -- every QC figure built from
them is experimental-only; there is no pooled-reference / run-drift comparison
available from this matrix.

## Batch confounding re-confirmation

`batch_confounding_report.json` (written by the same script) re-runs
`common.batch_correct.assess_batch_confounding(delivered_log, batch_column="Cohort",
covariate_columns=["Condition"])` on the delivered matrix. This **RE-CONFIRMS** the
Cohort <-> Condition confound already recorded in findings
[0001](../../findings/0001-group-aliased-with-cohort.md) and
[0004](../../findings/0004-batch-correction-extrapolated-single-cohort.md) -- it is a
corroborating check computed for provenance, not a new finding.

## Regenerating

```
./.venv/bin/python scripts/scratch/qc_prep.py
```

See `--help` for parameters. Figure-generator jobs should
`common.dataset_io.load_dataset(...)` the state they need from this directory rather
than re-deriving it.
