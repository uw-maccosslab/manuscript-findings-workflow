"""COHORT/SITE confound probe — observed site-AUC vs its label-shuffle null.

Renders the CACHED cohort-site confound-probe result (it recomputes NOTHING —
conventions/results-cache.md; a *load -> plot* job):

  * elastic-net logistic, SITE among the NON-ADD samples
    (``results/classification_cohort_site/b31f645d39ee/``): Cohort LevelMap
    **Kerchner (n=7, positive) vs Poston (n=20, negative)** among the 30 non-ADD
    samples (3 other-site samples dropped -> n=27). Nested ROC-AUC 0.34; observed
    (fixed-hp) AUC 0.475 vs a 1000-permutation label-shuffle null (mean 0.500,
    95% [0.475, 0.517], max 0.767); permutation p = 0.983; only 1 feature selected.

The honest question this figure answers: can the DELIVERED (already batch-corrected)
proteome predict collection SITE where site is NOT aliased with disease? The answer is
the OPPOSITE of the ADD classifier: the observed site AUC sits INSIDE its own null
(p = 0.98) -> the proteome does NOT predict collection site among non-ADD samples. To
make that contrast concrete the figure carries a second reference panel: the PRIMARY,
null-PASSED ADD-vs-Other elastic net (``results/classification/902c1e98c018/``;
observed AUC 0.956, null mean 0.502, p = 0.001) whose observed AUC lands far to the
RIGHT of its null. Same x-axis (ROC-AUC) on both panels, so the eye reads directly:
site collapses to null, ADD beats it.

Two stacked panels sharing the AUC x-axis:
  * TOP  — SITE among non-ADD (Kerchner vs Poston): null AUC histogram + 95% band +
    chance (0.5) + null mean; the observed AUC marked as a bold line SITTING INSIDE the
    null. p = 0.983.
  * BOTTOM — ADD vs Other (reference, null-passed): its null AUC histogram + the
    observed AUC (vermillion — the ADD `group` registry color) marked far to the right.
    p = 0.001.

The key (null / 95% band / chance / null mean / the two observed markers) is a separate
legend image (``<base>.legend.{svg,png}``), so it never overlaps the histograms.

LOAD-BEARING HONESTY (title + caption; findings 0001/0004/0006), stated BOTH directions
so nothing is overclaimed:
  * The delivered matrix is ALREADY limma removeBatchEffect-corrected on Cohort. A
    null-level site AUC means that correction removed detectable WITHIN-non-ADD site
    structure -- it does NOT clear the ADD result.
  * ADD is single-cohort (100% Kerchner, finding 0001), so ITS correction was
    EXTRAPOLATED, not estimated (finding 0004). Absence of residual site signal among
    cohort-spanning groups CANNOT prove the ADD separation is disease rather than site.
  * Small / imbalanced (n=27, 7 vs 20) -> a bounded probe, wide CIs. EXPLORATORY.

Run:  ./.venv/bin/python scripts/figures/plot_site_probe_null.py
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure

# Reuse the project shared figure-IO / color / result-cache modules, and the scratch
# analysis module that defines the cached result dataclass (type-directed reload).
_SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_SCRIPTS))
sys.path.insert(0, str(_SCRIPTS / "scratch"))
from classification import ClassificationResult  # noqa: E402
from common.colors import get_color  # noqa: E402
from common.figure_io import publication_style, save_figure  # noqa: E402
from common.result_io import load_cached_result  # noqa: E402

LOGGER = logging.getLogger("plot_site_probe_null")

__script_meta__: dict[str, object] = {
    "task": "plot-site-probe-null",
    "kind": "figure",
    "provides": ["build_site_probe_figure", "main"],
    "uses": [
        "common.figure_io",
        "common.colors",
        "common.result_io",
        "classification",
    ],
    "seeded_from": {"template": "classification-figures", "version": "0.2"},
    "description": (
        "Cohort/site confound probe: observed site ROC-AUC vs its label-shuffle null, "
        "rendered from the cached result (classification_cohort_site/b31f645d39ee). "
        "Two panels sharing the AUC x-axis -- SITE among non-ADD (Kerchner vs Poston; "
        "observed 0.475 INSIDE the null, p=0.983 -> collapses to null) over the "
        "null-passed ADD-vs-Other reference (classification/902c1e98c018; observed "
        "0.956, p=0.001 -> beats null). Null histograms + 95% band + chance(0.5); "
        "observed markers as bold lines (ADD in the group registry vermillion). "
        "Separate legend image. Load -> plot; recomputes nothing. Both-directions "
        "honesty: null-level site AUC does NOT clear the extrapolated ADD correction "
        "(findings 0001/0004). EXPLORATORY."
    ),
}

# Cached-result identity (results/manifest.md; provenance for the finding).
DEFAULT_CACHE_ROOT = Path("results")
SITE_ANALYSIS = "classification_cohort_site"
SITE_FINGERPRINT = "b31f645d39ee"
ADD_ANALYSIS = "classification"
ADD_FINGERPRINT = "902c1e98c018"  # null-PASSED ADD-vs-Other EN (the reference)
DATA_VERSION = "protein-data.csv sha256:90a0462eb03a02bd"

DEFAULT_OUTPUT_DIR = (
    Path("figures") / "analysis" / "classification" / "cohort-site-probe"
)
DEFAULT_BASE_NAME = "site-probe-null"
DEFAULT_REGISTRY = Path("state") / "color_registry.json"

# House-style (non-categorical) constants — these encode no metadata category, so they
# are house-style constants, not registry slots (mirrors lib/figures/classification).
_NULL_FILL = "#BBBBBB"  # null-distribution histogram fill (neutral gray)
_NULL_EDGE = "#7A7A7A"  # null-distribution histogram edge
_BAND_COLOR = "#999999"  # central 95% of the null (shaded span)
_CHANCE_COLOR = "#555555"  # chance line at AUC = 0.5
_NULLMEAN_COLOR = "#333333"  # null mean marker
_SITE_OBS_COLOR = "#0072B2"  # observed SITE AUC — Okabe-Ito blue house accent

_N_BINS = 34
_XLIM = (0.10, 1.0)


def _draw_panel(
    ax: Axes,
    *,
    panel_title: str,
    null_aucs: np.ndarray,
    observed: float,
    obs_color: str,
    verdict: str,
    verdict_xy: tuple[float, float],
    verdict_ha: str,
    info_lines: list[str],
) -> None:
    """One classifier's null-vs-observed panel: null histogram + observed marker."""
    lo, hi = float(np.percentile(null_aucs, 2.5)), float(np.percentile(null_aucs, 97.5))
    null_mean = float(null_aucs.mean())

    # Central 95% of the null (shaded), then the null AUC histogram on top.
    ax.axvspan(lo, hi, color=_BAND_COLOR, alpha=0.18, zorder=1)
    ax.hist(
        null_aucs,
        bins=_N_BINS,
        range=_XLIM,
        density=True,
        color=_NULL_FILL,
        edgecolor=_NULL_EDGE,
        linewidth=0.4,
        zorder=2,
    )
    # Chance (0.5) and the null mean — thin reference lines.
    ax.axvline(0.5, color=_CHANCE_COLOR, ls="--", lw=1.2, zorder=3)
    ax.axvline(null_mean, color=_NULLMEAN_COLOR, ls=":", lw=1.4, zorder=3)
    # The observed AUC — the load-bearing marker.
    ax.axvline(observed, color=obs_color, lw=3.0, zorder=5)

    ax.set_xlim(*_XLIM)
    ax.set_ylabel("null density")
    ax.set_title(panel_title, fontsize=12, fontweight="bold", loc="left")

    # The verdict pinned by the observed line (terse: which side of the null it lands).
    ax.text(
        verdict_xy[0],
        verdict_xy[1],
        verdict,
        transform=ax.transAxes,
        ha=verdict_ha,
        va="top",
        fontsize=8.5,
        fontweight="bold",
        color=obs_color,
    )
    # Terse load-bearing box (observed / null summary / p / n) — a marker, not prose.
    ax.text(
        0.015,
        0.96,
        "\n".join(info_lines),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.8,
        bbox={"boxstyle": "round,pad=0.35", "fc": "white", "ec": "lightgray"},
    )


def _build_legend(add_obs_color: str) -> Figure:
    """Standalone key for the null-vs-observed panels, as its own figure image."""
    with publication_style():
        lf, la = plt.subplots(figsize=(4.8, 2.4))
        la.axis("off")
        handles = [
            mpatches.Patch(
                fc=_NULL_FILL, ec=_NULL_EDGE, label="label-shuffle null AUCs"
            ),
            mpatches.Patch(fc=_BAND_COLOR, alpha=0.18, label="central 95% of null"),
            mlines.Line2D(
                [], [], color=_CHANCE_COLOR, ls="--", lw=1.2, label="chance (AUC = 0.5)"
            ),
            mlines.Line2D(
                [], [], color=_NULLMEAN_COLOR, ls=":", lw=1.4, label="null mean"
            ),
            mlines.Line2D(
                [], [], color=_SITE_OBS_COLOR, lw=3.0,
                label="observed site AUC (non-ADD)",
            ),
            mlines.Line2D(
                [], [], color=add_obs_color, lw=3.0,
                label="observed ADD-vs-Other AUC",
            ),
        ]
        la.legend(
            handles=handles,
            loc="center",
            frameon=False,
            fontsize=10,
            handletextpad=0.8,
            labelspacing=0.8,
            title="site-probe null key",
            title_fontsize=11,
        )
    return lf


def build_site_probe_figure(
    site: ClassificationResult, add: ClassificationResult, *, registry: Path
) -> tuple[Figure, Figure, str]:
    """Two-panel site-probe null figure (site over ADD reference) + its legend."""
    for name, r in (("site probe", site), ("ADD reference", add)):
        if r.null_aucs is None or r.observed_auc is None or r.null_p is None:
            raise ValueError(f"cached {name} result has no null distribution to plot.")

    # ADD observed marker uses the diagnostic-group registry color for ADD (vermillion),
    # so "ADD" reads the same here as in every other figure. Read-only (persist=False).
    add_color = get_color("group", "ADD", registry_path=registry, persist=False)

    site_null = np.asarray(site.null_aucs, dtype=float)
    add_null = np.asarray(add.null_aucs, dtype=float)
    site_lo, site_hi = np.percentile(site_null, [2.5, 97.5])
    add_lo, add_hi = np.percentile(add_null, [2.5, 97.5])

    site_info = [
        f"observed AUC = {site.observed_auc:.3f}",
        f"null mean = {site_null.mean():.3f}   95% [{site_lo:.3f}, {site_hi:.3f}]",
        f"permutation p = {site.null_p:.3f}   (1000 shuffles)",
        f"n = {site.n_samples}  ({site.n_positive} Kerchner / "
        f"{site.n_negative} Poston)",
        "1 feature selected (essentially nothing)",
        "small/imbalanced probe: null-level AUC BOUNDS",
        "residual site leakage — does not prove none",
    ]
    add_info = [
        f"observed AUC = {add.observed_auc:.3f}",
        f"null mean = {add_null.mean():.3f}   95% [{add_lo:.3f}, {add_hi:.3f}]",
        f"permutation p = {add.null_p:.3f}   (1000 shuffles)",
        f"n = {add.n_samples}  ({add.n_positive} ADD / {add.n_negative} Other)",
    ]

    with publication_style():
        fig, (ax_site, ax_add) = plt.subplots(
            2, 1, figsize=(9.2, 7.4), sharex=True
        )
        try:
            _draw_panel(
                ax_site,
                panel_title="Site among non-ADD  ·  Kerchner vs Poston",
                null_aucs=site_null,
                observed=float(site.observed_auc),
                obs_color=_SITE_OBS_COLOR,
                verdict="observed INSIDE null →\ncollapses to null (p = 0.98)",
                verdict_xy=(0.995, 0.60),
                verdict_ha="right",
                info_lines=site_info,
            )
            _draw_panel(
                ax_add,
                panel_title="ADD vs Other  ·  reference (null-passed)",
                null_aucs=add_null,
                observed=float(add.observed_auc),
                obs_color=add_color,
                verdict="observed BEYOND null\nbeats null (p = 0.001)",
                verdict_xy=(0.905, 0.96),
                verdict_ha="right",
                info_lines=add_info,
            )
            ax_add.set_xlabel(
                "ROC-AUC  (elastic-net logistic; label-shuffle null vs observed)"
            )

            fig.suptitle(
                "Cohort/site confound probe — does the corrected proteome predict "
                "collection site?",
                fontsize=13.5,
                fontweight="bold",
            )
            # One terse caveat marker (both directions); full explanation -> caption.
            fig.text(
                0.5,
                0.005,
                "EXPLORATORY · null-level site AUC does NOT clear the ADD result — ADD "
                "is single-cohort, so its correction was EXTRAPOLATED, not estimated "
                "(findings 0001/0004)",
                ha="center",
                va="bottom",
                fontsize=8.5,
                color="#8a3d00",
            )
            fig.tight_layout(rect=(0, 0.03, 1, 0.96))
        except BaseException:
            plt.close(fig)
            raise

    legend_fig = _build_legend(add_color)
    return fig, legend_fig, add_color


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(name)s: %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-root", type=Path, default=DEFAULT_CACHE_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--base-name", type=str, default=DEFAULT_BASE_NAME)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    args = parser.parse_args(argv)

    site = load_cached_result(
        ClassificationResult,
        cache_root=args.cache_root,
        analysis=SITE_ANALYSIS,
        fingerprint=SITE_FINGERPRINT,
    )
    add = load_cached_result(
        ClassificationResult,
        cache_root=args.cache_root,
        analysis=ADD_ANALYSIS,
        fingerprint=ADD_FINGERPRINT,
    )
    LOGGER.info(
        "loaded cached results: SITE %s (observed AUC %.3f, p=%.3f), "
        "ADD %s (observed AUC %.3f, p=%.3f); data_version=%s",
        SITE_FINGERPRINT,
        site.observed_auc,
        site.null_p,
        ADD_FINGERPRINT,
        add.observed_auc,
        add.null_p,
        DATA_VERSION,
    )

    fig, legend_fig, add_color = build_site_probe_figure(
        site, add, registry=args.registry
    )
    artifacts = save_figure(
        fig, args.output_dir, args.base_name, legend_fig=legend_fig, dpi=300
    )
    LOGGER.info("ADD observed marker color (group=ADD, registry): %s", add_color)
    for p in (artifacts.svg, artifacts.png, artifacts.legend_svg, artifacts.legend_png):
        LOGGER.info("wrote %s", p)


if __name__ == "__main__":
    main()
