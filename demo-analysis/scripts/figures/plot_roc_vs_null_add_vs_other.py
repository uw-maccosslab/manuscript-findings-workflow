"""Permutation-null vs observed AUC for the ADD-vs-Other classifiers (NULL PASS).

Renders the two CACHED null-pass classifier results (it recomputes NOTHING —
conventions/results-cache.md; a pure *load -> plot* job):

  * elastic-net logistic  — ``results/classification/902c1e98c018/``
    (observed nested ROC-AUC 0.940 ± 0.090; fixed-hp observed 0.956; label-shuffle
    NULL AUC mean 0.502, 95% [0.244, 0.778], max 0.883; permutation p = 0.001);
  * linear SVM            — ``results/classification_svm/a0e23fb510a0/``
    (observed nested 0.963 ± 0.067; fixed-hp observed 0.972; NULL AUC mean 0.500,
    95% [0.228, 0.756]; permutation p = 0.001).

Both were produced on the delivered log2 40x2340 analysis subset, ADD (n=10, positive)
vs pooled Other {PDD,PDCN,HCN} (n=30, negative), with a label-shuffle null of 1000
permutations performed at the UNIT (individual) level and evaluated at fixed
hyperparameters. Seeded from ``lib/figures/classification`` (``plot_null``) — this
reuses that template's null-histogram + observed-AUC-marker + empirical-p construction,
extended to two side-by-side panels so the two models compare directly.

Two panels (elastic net | linear SVM); each draws the 1000 permuted-label AUCs as a
histogram (neutral gray), a chance reference at 0.5, and the observed AUC as a bold
vermillion line carrying the empirical p. The observed AUC sits FAR to the right of the
null (null centered at 0.5, observed ~0.96) — the visual message that above-null signal
exists. The key (null / chance / observed) is a SEPARATE legend image; the load-bearing
numbers (nested + fixed-hp AUC, null mean/95%, p, N, caveat) are a terse on-panel box.

LOAD-BEARING HONESTY (title + caption + finding prose; findings 0001/0004):
  * Observed AUC beats the label-shuffle null (p = 0.001) -> the performance is NOT an
    artifact of leakage and NOT pure overfitting: a real above-null signal separates ADD
    from Other in these data.
  * This does NOT establish disease biology. ADD is 100% Kerchner cohort (finding 0001),
    so the label the null permutes is BOTH "ADD" and "Kerchner" (the SAME partition of
    these 40 samples). An above-null result cannot separate disease from collection
    site.
  * Small n (10 ADD / 30 Other) -> EXPLORATORY.

Run:  ./.venv/bin/python scripts/figures/plot_roc_vs_null_add_vs_other.py
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
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# Reuse the project shared figure-IO + result-cache modules, and the scratch analysis
# modules that define the cached result dataclasses (type-directed reload).
_SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_SCRIPTS))
sys.path.insert(0, str(_SCRIPTS / "scratch"))
from classification import ClassificationResult  # noqa: E402
from classification_svm import SVMClassificationResult  # noqa: E402
from common.figure_io import publication_style, save_figure  # noqa: E402
from common.result_io import load_cached_result  # noqa: E402

LOGGER = logging.getLogger("plot_roc_vs_null_add_vs_other")

__script_meta__: dict[str, object] = {
    "task": "plot-roc-vs-null-add-vs-other",
    "kind": "figure",
    "provides": ["build_null_figure", "build_legend_figure", "main"],
    "uses": [
        "common.figure_io",
        "common.result_io",
        "classification",
        "classification_svm",
    ],
    "seeded_from": {"template": "classification-figures", "version": "0.2"},
    "description": (
        "Two-panel permutation-null vs observed AUC (elastic net | linear SVM) for the "
        "NULL-PASS ADD-vs-Other classifiers, rendered from the cached results "
        "(classification/902c1e98c018, classification_svm/a0e23fb510a0). Histogram of "
        "1000 unit-level label-shuffle AUCs + chance @0.5 + observed-AUC marker with "
        "empirical p. Observed AUC (~0.96) sits far above the null (mean ~0.50) -> not "
        "leakage / not pure overfitting; BUT ADD == Kerchner cohort confound (finding "
        "0001) means above-null != disease-vs-site. EXPLORATORY. Load -> plot; "
        "recomputes nothing. Separate legend image."
    ),
}

# Cached-result identity (results/manifest.md; provenance for the finding).
DEFAULT_CACHE_ROOT = Path("results")
EN_ANALYSIS = "classification"
EN_FINGERPRINT = "902c1e98c018"
SVM_ANALYSIS = "classification_svm"
SVM_FINGERPRINT = "a0e23fb510a0"
DATA_VERSION = "protein-data.csv sha256:90a0462eb03a02bd"

DEFAULT_OUTPUT_DIR = Path("figures") / "analysis" / "classification" / "add-vs-other"
DEFAULT_BASE_NAME = "roc-vs-null-add-vs-other"

# Colors. The observed-AUC marker is drawn in Okabe-Ito vermillion #D55E00 — the SAME
# hue the color registry assigns to the ADD group / primary contrast (state/
# color_registry.json: group.ADD, primary_contrast.ADD), so the marker for the ADD
# classifier's observed performance stays consistent with ADD's color elsewhere. The
# null histogram + chance line use Okabe-Ito neutral gray #999999 (registry
# primary_contrast."Other"), which reads as the no-signal reference. These follow the
# lib plot_null house style (_OBSERVED_COLOR / _NULL_FILL / _CHANCE_COLOR).
_OBSERVED_COLOR = "#D55E00"  # observed-AUC marker (== registry ADD)
_NULL_FILL = "#999999"  # null-AUC histogram (== registry "Other" / no-signal)
_CHANCE_COLOR = "#999999"  # chance reference at AUC = 0.5

# Shared x-limits so the two panels compare on one scale; the gap between the null mass
# and the observed marker is the message.
_XLIM = (0.1, 1.02)
_N_BINS = 25


def _draw_panel(
    ax: Axes,
    *,
    panel_title: str,
    null_aucs: np.ndarray,
    observed_auc: float,
    nested_auc: float,
    nested_auc_sd: float,
    null_p: float,
    n_positive: int,
    n_negative: int,
    positive_label: str,
    negative_label: str,
) -> None:
    """One model's null panel: null histogram + chance @0.5 + observed marker + box."""
    nulls = np.asarray(null_aucs, dtype=float)
    n_perm = int(nulls.size)
    null_mean = float(nulls.mean())
    lo, hi = (float(x) for x in np.percentile(nulls, [2.5, 97.5]))

    ax.hist(
        nulls,
        bins=_N_BINS,
        range=_XLIM,
        color=_NULL_FILL,
        edgecolor="white",
        zorder=2,
    )
    ax.axvline(0.5, color=_CHANCE_COLOR, ls="--", lw=1.2, zorder=3)
    ax.axvline(observed_auc, color=_OBSERVED_COLOR, lw=2.4, zorder=4)

    ax.set_xlim(*_XLIM)
    ax.set_xlabel("ROC AUC under permuted labels")
    ax.set_ylabel("permutations (count)")
    ax.set_title(panel_title, fontsize=12, fontweight="bold")

    # Terse load-bearing box (a marker set, not a paragraph): the honest nested-CV AUC,
    # the fixed-hp observed AUC the p is computed against, the null center/spread, the
    # empirical p, per-class N, and the EXPLORATORY caveat marker. Interpretation +
    # the confound live in the finding caption/prose, not on the canvas.
    text = (
        f"observed AUC = {observed_auc:.3f}  (p = {null_p:.3f})\n"
        f"nested-CV AUC = {nested_auc:.3f} ± {nested_auc_sd:.2f}\n"
        f"null: mean {null_mean:.3f}, 95% [{lo:.2f}, {hi:.2f}]  ({n_perm} perms)\n"
        f"{positive_label}: N={n_positive}  |  {negative_label}: N={n_negative}\n"
        f"EXPLORATORY"
    )
    ax.text(
        0.03,
        0.97,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.8,
        bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "lightgray"},
    )


def build_null_figure(
    en: ClassificationResult, svm: SVMClassificationResult
) -> Figure:
    """Two-panel permutation-null vs observed AUC (elastic net | linear SVM)."""
    for name, res in (("elastic net", en), ("linear SVM", svm)):
        if res.null_aucs is None or res.observed_auc is None or res.null_p is None:
            raise ValueError(
                f"the {name} cached result has no label-shuffle null "
                f"(run_null=False?); nothing to plot."
            )

    with publication_style():
        fig, (ax_en, ax_svm) = plt.subplots(1, 2, figsize=(12.6, 5.6))
        try:
            _draw_panel(
                ax_en,
                panel_title="Elastic-net logistic",
                null_aucs=en.null_aucs,
                observed_auc=en.observed_auc,
                nested_auc=en.cv_auc,
                nested_auc_sd=en.cv_auc_sd,
                null_p=en.null_p,
                n_positive=en.n_positive,
                n_negative=en.n_negative,
                positive_label=en.positive_label,
                negative_label=en.negative_label,
            )
            _draw_panel(
                ax_svm,
                panel_title="Linear SVM",
                null_aucs=svm.null_aucs,
                observed_auc=svm.observed_auc,
                nested_auc=svm.cv_auc,
                nested_auc_sd=svm.cv_auc_sd,
                null_p=svm.null_p,
                n_positive=svm.n_positive,
                n_negative=svm.n_negative,
                positive_label=svm.positive_label,
                negative_label=svm.negative_label,
            )
            fig.suptitle(
                f"{en.positive_label} vs {en.negative_label} — label-shuffle null vs "
                f"observed AUC (elastic net vs linear SVM)",
                fontsize=14,
                fontweight="bold",
            )
            # One terse caveat line (a mandatory marker, not a paragraph): observed
            # beats null (not leakage / not pure overfitting) BUT the confound means
            # above-null cannot separate disease from site. Full reasoning -> finding.
            fig.text(
                0.5,
                0.005,
                "Label-shuffle null (1000 perms, unit-level) · observed ≫ null "
                "(p=0.001) → NOT leakage / NOT pure overfitting  ·  BUT ADD = 100% "
                "Kerchner cohort (finding 0001): above-null ≠ disease vs site  ·  "
                "n=40 · EXPLORATORY",
                ha="center",
                va="bottom",
                fontsize=8.5,
                color="#8a3d00",
            )
            fig.tight_layout(rect=(0, 0.045, 1, 0.96))
        except BaseException:
            plt.close(fig)
            raise
    return fig


def build_legend_figure() -> Figure:
    """Standalone swatch/line key (null histogram / chance / observed marker).

    Rendered as its own figure (conventions/visualization.md) so it never overlaps the
    panels; embedded beneath the figure in the finding.
    """
    handles = [
        Patch(
            facecolor=_NULL_FILL,
            edgecolor="white",
            label="null AUC (permuted labels)",
        ),
        Line2D(
            [0], [0], color=_CHANCE_COLOR, ls="--", lw=1.2, label="chance (AUC = 0.5)"
        ),
        Line2D([0], [0], color=_OBSERVED_COLOR, lw=2.4, label="observed AUC (with p)"),
    ]
    with publication_style():
        fig = plt.figure(figsize=(4.6, 1.1))
        try:
            fig.legend(
                handles=handles, loc="center", ncol=1, frameon=False, fontsize=10
            )
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
        "loaded cached NULL-pass results: EN %s (obs AUC %.3f, p %.3f), "
        "SVM %s (obs AUC %.3f, p %.3f); data_version=%s",
        EN_FINGERPRINT,
        en.observed_auc,
        en.null_p,
        SVM_FINGERPRINT,
        svm.observed_auc,
        svm.null_p,
        DATA_VERSION,
    )

    fig = build_null_figure(en, svm)
    legend = build_legend_figure()
    artifacts = save_figure(
        fig, args.output_dir, args.base_name, legend_fig=legend, dpi=300
    )
    for p in (artifacts.svg, artifacts.png, artifacts.legend_svg, artifacts.legend_png):
        LOGGER.info("wrote %s", p)


if __name__ == "__main__":
    main()
