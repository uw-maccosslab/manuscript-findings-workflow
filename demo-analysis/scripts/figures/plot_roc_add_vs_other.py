"""Nested-CV ROC curves for the PRIMARY ADD-vs-Other classifiers — NULL PASS.

Renders the two CACHED, NULL-PASSED classifier results (it recomputes NOTHING —
conventions/results-cache.md; a *load -> plot* job). These SUPERSEDE the first-pass
caches (``cbd809ac6168`` / ``cba84421b85e``); this figure replaces the earlier
first-pass ROC render:

  * elastic-net logistic  — ``results/classification/902c1e98c018/``
    (nested ROC-AUC 0.940 ± 0.090 over 25 outer folds; per-fold range 0.75-1.00;
    run_null=True, converged at max_iter=20000);
  * linear SVM            — ``results/classification_svm/a0e23fb510a0/``
    (nested ROC-AUC 0.963 ± 0.067; per-fold range 0.75-1.00; run_null=True).

Both were produced by ``scripts/scratch/classify_add_vs_other_null.py`` on the delivered
log2 40x2340 analysis subset, ADD (n=10, positive) vs pooled Other {PDD,PDCN,HCN}
(n=30, negative), stratified 5-fold x 5 repeats at the individual (row) level, SAME
seed / SAME outer folds. Seeded from ``lib/figures/classification`` (``plot_roc``) —
this reuses that template's per-fold ROC construction (roc_curve -> interpolate onto a
common FPR grid -> mean ± 1 SD), extended to two side-by-side panels so the two models
compare directly and each panel shows the honest per-fold spread.

Two panels (elastic net | linear SVM); each draws the 25 per-fold ROC curves (thin,
translucent), their mean ROC (bold) with a ±1 SD band, and a chance diagonal. The ROC
key sits ON-AXES in the empty lower-right corner (the documented separate-legend
exception for ROC, conventions/visualization.md), so this figure emits NO separate
legend image. The SVM panel is drawn from the decision-function MARGIN (uncalibrated;
AUC is rank-based, so no probability is needed); the elastic net from predicted
probability.

LOAD-BEARING HONESTY (title + caption; findings 0001/0004):
  * The label-shuffle null HAS now been run and PASSED (permutation p=0.001 for each
    model; observed AUC far above null ~0.50) -> the discrimination is NOT leakage or
    overfitting. Still EXPLORATORY (small n; confound below).
  * HARD CONFOUND (UNRESOLVED) — ADD is 100% Kerchner cohort (finding 0001). Among these
    40 samples "ADD" and "Kerchner" are the SAME partition, so above-null discrimination
    canNOT be attributed to disease vs collection-site; the unit-level shuffle permutes
    that shared partition, so the passing null does not resolve the confound.
  * Small n (10 ADD / 30 Other) -> wide per-fold AUC (0.75-1.00).

Run:  ./.venv/bin/python scripts/figures/plot_roc_add_vs_other.py
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from sklearn.metrics import auc, roc_curve

# Reuse the project shared figure-IO + result-cache modules, and the scratch analysis
# modules that define the cached result dataclasses (type-directed reload).
_SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_SCRIPTS))
sys.path.insert(0, str(_SCRIPTS / "scratch"))
from classification import ClassificationResult  # noqa: E402
from classification_svm import SVMClassificationResult  # noqa: E402
from common.figure_io import publication_style, save_figure  # noqa: E402
from common.result_io import load_cached_result  # noqa: E402

LOGGER = logging.getLogger("plot_roc_add_vs_other")

__script_meta__: dict[str, object] = {
    "task": "plot-roc-add-vs-other-nullpass",
    "kind": "figure",
    "provides": ["build_roc_figure", "main"],
    "uses": [
        "common.figure_io",
        "common.result_io",
        "classification",
        "classification_svm",
    ],
    "seeded_from": {"template": "classification-figures", "version": "0.2"},
    "description": (
        "Two-panel nested-CV ROC (elastic net | linear SVM) for the NULL-PASSED "
        "ADD-vs-Other classifiers, rendered from the current cached results "
        "(classification/902c1e98c018, classification_svm/a0e23fb510a0; these "
        "supersede the first-pass cbd809ac6168 / cba84421b85e). Per-fold curves + "
        "mean ± 1 SD + chance diagonal; per-fold AUC range annotated. Null run and "
        "PASSED (p=0.001) -> not leakage/overfitting, still EXPLORATORY; ADD == "
        "Kerchner cohort confound UNRESOLVED. ROC key on-axes (no separate legend "
        "image). Load -> plot; recomputes nothing."
    ),
}

# Cached-result identity (results/manifest.md; provenance for the finding). These are
# the CURRENT null-pass results; they supersede the first-pass cbd809ac6168 /
# cba84421b85e that the earlier ROC render used.
DEFAULT_CACHE_ROOT = Path("results")
EN_ANALYSIS = "classification"
EN_FINGERPRINT = "902c1e98c018"
SVM_ANALYSIS = "classification_svm"
SVM_FINGERPRINT = "a0e23fb510a0"
DATA_VERSION = "protein-data.csv sha256:90a0462eb03a02bd"
# Empirical label-shuffle permutation p (identical for both models; results/manifest).
PERM_P = 0.001

DEFAULT_OUTPUT_DIR = Path("figures") / "analysis" / "classification" / "add-vs-other"
DEFAULT_BASE_NAME = "roc-add-vs-other-nullpass"

# House-style (non-categorical) ROC colors, matching lib/figures/classification. These
# encode no metadata category, so they are house-style constants, not registry slots.
_ROC_COLOR = "#0072B2"  # Okabe-Ito blue — mean ROC + band + per-fold curves
_CHANCE_COLOR = "#999999"  # Okabe-Ito neutral gray — chance diagonal

_ROC_GRID = np.linspace(0.0, 1.0, 100)


def _fold_roc(
    folds: object, score_attr: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-fold interpolated TPRs + per-fold AUCs, mirroring lib plot_roc.

    Each outer fold contributes one ROC curve (its held-out ``score_attr`` — predicted
    probability for the elastic net, decision-function margin for the SVM); curves are
    interpolated onto a common FPR grid. Returns the per-fold TPR stack, the mean TPR
    (with mean_tpr[-1] pinned to 1.0), and the per-fold AUCs.
    """
    tprs: list[np.ndarray] = []
    aucs: list[float] = []
    for fold in folds:  # type: ignore[attr-defined]
        score = getattr(fold, score_attr)
        fpr, tpr, _ = roc_curve(fold.y_true, score)
        interp = np.interp(_ROC_GRID, fpr, tpr)
        interp[0] = 0.0
        tprs.append(interp)
        aucs.append(float(auc(fpr, tpr)))
    tpr_stack = np.asarray(tprs, dtype=float)
    mean_tpr = tpr_stack.mean(axis=0)
    mean_tpr[-1] = 1.0
    return tpr_stack, mean_tpr, np.asarray(aucs, dtype=float)


def _draw_panel(
    ax: Axes,
    *,
    panel_title: str,
    tpr_stack: np.ndarray,
    mean_tpr: np.ndarray,
    auc_arr: np.ndarray,
    nested_auc: float,
    nested_auc_sd: float,
    n_positive: int,
    n_negative: int,
    positive_label: str,
    negative_label: str,
    score_note: str,
) -> None:
    """One model's ROC panel: per-fold curves + mean ± 1 SD + chance, key on-axes."""
    n_folds = tpr_stack.shape[0]
    sd_tpr = tpr_stack.std(axis=0)
    fold_mean = float(auc_arr.mean())
    fold_sd = float(auc_arr.std())

    ax.plot([0, 1], [0, 1], "--", color=_CHANCE_COLOR, lw=1.2, label="chance", zorder=1)
    # Per-fold ROC curves — the honest spread behind the mean.
    for i in range(n_folds):
        ax.plot(
            _ROC_GRID,
            tpr_stack[i],
            color=_ROC_COLOR,
            lw=0.6,
            alpha=0.18,
            zorder=2,
            label=f"per-fold ROC ({n_folds} folds)" if i == 0 else None,
        )
    ax.fill_between(
        _ROC_GRID,
        np.clip(mean_tpr - sd_tpr, 0.0, 1.0),
        np.clip(mean_tpr + sd_tpr, 0.0, 1.0),
        color=_ROC_COLOR,
        alpha=0.20,
        label="± 1 SD (folds)",
        zorder=3,
    )
    ax.plot(
        _ROC_GRID,
        mean_tpr,
        color=_ROC_COLOR,
        lw=2.4,
        label=f"mean ROC (AUC = {fold_mean:.3f} ± {fold_sd:.3f})",
        zorder=4,
    )

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title(panel_title, fontsize=12, fontweight="bold")
    # Key (chance / per-fold / ±SD / mean ROC) pinned to the bottom of the empty
    # lower-right corner (documented on-axes ROC exception).
    ax.legend(
        loc="lower right", bbox_to_anchor=(1.0, 0.0), fontsize=8, borderaxespad=0.6
    )

    # Terse load-bearing box STACKED ABOVE the key so the two never overlap (nested AUC
    # + per-fold RANGE = the honest uncertainty + per-class N + the EXPLORATORY caveat
    # marker). A short marker, not a paragraph.
    text = (
        f"nested ROC-AUC = {nested_auc:.3f} ± {nested_auc_sd:.2f}\n"
        f"per-fold AUC: {auc_arr.min():.2f}–{auc_arr.max():.2f}  ({n_folds} folds)\n"  # noqa: RUF001
        f"{positive_label}: N={n_positive}  |  {negative_label}: N={n_negative}\n"
        f"{score_note}\n"
        f"EXPLORATORY · above null (p={PERM_P:.3f})"
    )
    ax.text(
        0.97,
        0.30,
        text,
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.6,
        bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "lightgray"},
    )


def build_roc_figure(
    en: ClassificationResult, svm: SVMClassificationResult
) -> Figure:
    """Two-panel nested-CV ROC (elastic net | linear SVM) from the cached results."""
    if not en.fold_predictions or not svm.fold_predictions:
        raise ValueError("a cached result has no fold predictions for a ROC.")
    en_stack, en_mean, en_aucs = _fold_roc(en.fold_predictions, "y_prob")
    svm_stack, svm_mean, svm_aucs = _fold_roc(svm.fold_predictions, "y_score")

    with publication_style():
        fig, (ax_en, ax_svm) = plt.subplots(1, 2, figsize=(12.6, 6.5))
        try:
            _draw_panel(
                ax_en,
                panel_title="Elastic-net logistic",
                tpr_stack=en_stack,
                mean_tpr=en_mean,
                auc_arr=en_aucs,
                nested_auc=en.cv_auc,
                nested_auc_sd=en.cv_auc_sd,
                n_positive=en.n_positive,
                n_negative=en.n_negative,
                positive_label=en.positive_label,
                negative_label=en.negative_label,
                score_note="scores = predicted probability",
            )
            _draw_panel(
                ax_svm,
                panel_title="Linear SVM",
                tpr_stack=svm_stack,
                mean_tpr=svm_mean,
                auc_arr=svm_aucs,
                nested_auc=svm.cv_auc,
                nested_auc_sd=svm.cv_auc_sd,
                n_positive=svm.n_positive,
                n_negative=svm.n_negative,
                positive_label=svm.positive_label,
                negative_label=svm.negative_label,
                score_note="scores = SVM margin (uncalibrated)",
            )
            fig.suptitle(
                f"{en.positive_label} vs {en.negative_label} — nested-CV ROC "
                f"(elastic net vs linear SVM)",
                fontsize=14,
                fontweight="bold",
            )
            # One terse caveat line (a mandatory marker, not a paragraph): the full
            # confound explanation lives in the finding caption.
            fig.text(
                0.5,
                0.005,
                f"EXPLORATORY · beats label-shuffle null (p={PERM_P:.3f}) → not "
                "leakage/overfitting  ·  BUT hard confound UNRESOLVED: ADD = 100% "
                "Kerchner cohort (finding 0001) — above-null ≠ disease vs collection "
                "site",
                ha="center",
                va="bottom",
                fontsize=9,
                color="#8a3d00",
            )
            fig.tight_layout(rect=(0, 0.035, 1, 0.97))
        except BaseException:
            plt.close(fig)
            raise
    return fig


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(name)s: %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-root", type=Path, default=DEFAULT_CACHE_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--base-name", type=str, default=DEFAULT_BASE_NAME)
    args = parser.parse_args(argv)

    en = load_cached_result(
        ClassificationResult,
        cache_root=args.cache_root,
        analysis=EN_ANALYSIS,
        fingerprint=EN_FINGERPRINT,
    )
    svm = load_cached_result(
        SVMClassificationResult,
        cache_root=args.cache_root,
        analysis=SVM_ANALYSIS,
        fingerprint=SVM_FINGERPRINT,
    )
    LOGGER.info(
        "loaded cached results: EN %s (AUC %.3f), SVM %s (AUC %.3f); data_version=%s",
        EN_FINGERPRINT,
        en.cv_auc,
        SVM_FINGERPRINT,
        svm.cv_auc,
        DATA_VERSION,
    )

    fig = build_roc_figure(en, svm)
    artifacts = save_figure(fig, args.output_dir, args.base_name, dpi=300)
    for p in (artifacts.svg, artifacts.png):
        LOGGER.info("wrote %s", p)
    LOGGER.info("ROC key on-axes (documented exception) — no separate legend image.")


if __name__ == "__main__":
    main()
