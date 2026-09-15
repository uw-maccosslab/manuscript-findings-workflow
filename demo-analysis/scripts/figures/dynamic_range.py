"""Dynamic-range / rank-abundance QC figures for a :class:`Dataset`.

PROJECT COPY of the ``dynamic-range`` lib/ template (``lib/figures/dynamic_range.py``
@ 0.1). The compute + plot body is the vetted template, adapted only where this project's
layout differs from the template's assumptions:

  * imports the project's shared figure helpers from ``common.*`` (``common.colors``,
    ``common.figure_io``) rather than the template's ``figures.*`` package — this project
    keeps its seeded figure modules under ``scripts/common/`` (see scripts/manifest.md);
  * appends a ``main()`` CLI runner that renders the **Stage-3 QC whole-cohort**
    dynamic-range figure for the delivered-linear processing state.

Held to the correctness charter: assume nothing, verify everything, fail loud. This
module makes NO study decisions in its compute/plot layer — the detection threshold, the
class split, and any highlights are the caller's; the QC render choices live in ``main``.

What it draws (one figure): features ranked by abundance (most -> least) on the x-axis vs
their log2 abundance on the y-axis -- the classic rank-abundance curve, showing the
**dynamic range** quantified (orders of magnitude) and whether a handful of hyper-abundant
features dominate. Whole-cohort mode (default): the per-feature median of detected values
with a shaded per-sample IQR band; optional ``highlight_features`` mark named proteins of
interest (empty at QC, populated downstream at Stage-4).

Scale is a HARD refuse: the dynamic range is a raw-linear property, so
:func:`compute_dynamic_range` raises :class:`DynamicRangeScaleError` on a non-linear
``Dataset.scale``. NOTE for this project: the linear matrix consumed here is the
**de-logged image of the delivered normalized + batch-corrected matrix** (produced by
``scripts/scratch/qc_prep.py`` via ``to_linear``), NOT raw instrument intensity. The
figure title/caption say so; read the curve as relative dynamic range on processed data.
"""

from __future__ import annotations

import argparse
import logging
import sys
import warnings
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

# Import the project's shared modules (seeded from the lib/ templates).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.colors import DEFAULT_REGISTRY_PATH, assign_colors
from common.data_loading import Dataset
from common.dataset_io import load_dataset
from common.figure_io import FigureArtifacts, publication_style, save_figure

__script_meta__: dict[str, object] = {
    "template": {"name": "dynamic-range", "version": "0.1"},
    "kind": "figure",
    "provides": [
        "DynamicRangeScaleError",
        "DynamicRangeResult",
        "DynamicRangePlot",
        "compute_dynamic_range",
        "plot_dynamic_range",
        "save_dynamic_range",
    ],
    "uses": ["common.data_loading", "common.dataset_io", "common.colors",
             "common.figure_io"],
    "seeded_from": {"template": "dynamic-range", "version": "0.1"},
    "description": (
        "Dynamic-range / rank-abundance QC figure from a linear Dataset: features "
        "ranked by abundance vs log2 abundance, showing the quantified dynamic range "
        "and any hyper-abundant dominance. Whole-cohort median + per-sample IQR band, "
        "optional leader-labelled highlights (empty at QC). Hard-refuses a non-linear "
        "scale. main() renders the Stage-3 QC delivered-linear state (de-logged image "
        "of the delivered normalized+batch-corrected matrix -- NOT raw intensity)."
    ),
}

DEFAULT_HIGHLIGHT_COLOR = "#0072B2"
_BAND_COLOR = "#bdbdbd"
_CURVE_COLOR = "#333333"

LOGGER = logging.getLogger("dynamic_range")


class DynamicRangeScaleError(ValueError):
    """Raised when the dynamic range is requested on a non-linear abundance scale.

    The rank-abundance view is read on a **linear matrix**: "detected"
    (``> min_intensity``) is a linear concept and the y-axis is the log2 of the linear
    abundance. On a log or centered scale those no longer hold.
    """


@dataclass(frozen=True)
class DynamicRangeResult:
    """The ranked per-feature abundance summary underlying a dynamic-range figure.

    Only features detected in at least one sample appear; arrays are ordered
    most -> least abundant.
    """

    feature_names_ranked: np.ndarray
    log2_median: np.ndarray
    log2_q25: np.ndarray
    log2_q75: np.ndarray
    dynamic_range_orders: float
    n_features_total: int
    n_features_detected: int

    @property
    def ranks(self) -> np.ndarray:
        """``(k,)`` ranks ``1..k`` aligned to the ranked arrays."""
        return np.arange(1, self.n_features_detected + 1)


@dataclass
class DynamicRangePlot:
    """A rendered dynamic-range figure plus its companion legend figure (or ``None``)."""

    figure: Figure
    legend_figure: Figure | None
    result: DynamicRangeResult
    color_map: dict[str, str]


# --------------------------------------------------------------------------- #
# Compute
# --------------------------------------------------------------------------- #


def compute_dynamic_range(
    dataset: Dataset, *, min_intensity: float = 0.0
) -> DynamicRangeResult:
    """Rank features by median detected abundance; return the ranked log2 summary.

    Raises :class:`DynamicRangeScaleError` unless ``dataset.scale == "linear"``.
    Features detected in no sample are excluded.
    """
    if dataset.scale != "linear":
        raise DynamicRangeScaleError(
            f"Dynamic range requires linear-scale abundances but the Dataset is on "
            f"scale {dataset.scale!r}. 'detected' (> {min_intensity}) and the log2 "
            f"y-axis are raw-linear properties. Run it on the linear matrix."
        )
    abundances = np.asarray(dataset.abundances, dtype=float)
    if abundances.ndim != 2:
        raise ValueError(
            f"abundances must be 2D (n_samples, n_features); got {abundances.shape}."
        )
    names = np.asarray(dataset.feature_names)
    n_total = abundances.shape[1]
    if names.shape[0] != n_total:
        raise ValueError(
            f"feature_names has {names.shape[0]} entries but abundances has {n_total} "
            f"features; they are not aligned."
        )
    if abundances.shape[0] < 1 or n_total < 1:
        raise ValueError(
            f"Dataset is empty ({abundances.shape[0]} samples x {n_total} features)."
        )

    median, q25, q75, detected_count = _detected_quartiles(abundances, min_intensity)
    keep = detected_count >= 1
    if not bool(keep.any()):
        raise ValueError(
            "no feature is detected in any sample (all <= min_intensity); none to rank."
        )

    med_keep = median[keep]
    order = np.argsort(med_keep)[::-1]
    med_ranked = med_keep[order]
    orders = (
        float(np.log10(med_ranked[0] / med_ranked[-1]))
        if med_ranked[-1] > 0
        else float("nan")
    )

    return DynamicRangeResult(
        feature_names_ranked=names[keep][order],
        log2_median=np.log2(med_ranked),
        log2_q25=np.log2(q25[keep][order]),
        log2_q75=np.log2(q75[keep][order]),
        dynamic_range_orders=orders,
        n_features_total=int(n_total),
        n_features_detected=int(keep.sum()),
    )


def _detected_quartiles(
    abundances: np.ndarray, min_intensity: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-feature median / q25 / q75 of detected (> min_intensity) values + count."""
    mask = np.isfinite(abundances) & (abundances > min_intensity)
    positive = np.where(mask, abundances, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        median = np.nanmedian(positive, axis=0)
        q25 = np.nanpercentile(positive, 25, axis=0)
        q75 = np.nanpercentile(positive, 75, axis=0)
    return median, q25, q75, mask.sum(axis=0)


# --------------------------------------------------------------------------- #
# Plot
# --------------------------------------------------------------------------- #


def plot_dynamic_range(
    dataset: Dataset,
    *,
    class_by: str | None = None,
    highlight_features: Mapping[str, str] | None = None,
    highlight_groups: Mapping[str, str] | None = None,
    highlight_markers: Mapping[str, str] | None = None,
    highlight_color: str = DEFAULT_HIGHLIGHT_COLOR,
    min_intensity: float = 0.0,
    log_rank: bool = False,
    title: str | None = None,
    legend_title: str | None = None,
    registry_path: str | Path = DEFAULT_REGISTRY_PATH,
    persist_colors: bool = True,
) -> DynamicRangePlot:
    """Draw the rank-abundance curve for one Dataset (whole-cohort or per-class).

    ``highlight_markers`` optionally maps each ``highlight_groups`` group to a matplotlib
    marker so the groups are distinguishable by shape as well as color (a second,
    colorblind-safe channel); it is echoed in the companion legend. Ignored when there
    are no grouped highlights.
    """
    if class_by is not None and highlight_features:
        raise ValueError(
            "class_by and highlight_features are mutually exclusive: per-class curves "
            "are ranked independently, so a highlight has no single rank. Pass one."
        )
    result = compute_dynamic_range(dataset, min_intensity=min_intensity)

    color_map: dict[str, str] = {}
    legend_figure: Figure | None = None
    with publication_style():
        fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
        try:
            if class_by is not None:
                color_map = _draw_per_class(
                    ax,
                    dataset,
                    class_by,
                    min_intensity=min_intensity,
                    registry_path=registry_path,
                    persist_colors=persist_colors,
                )
                legend_figure = _legend_figure(
                    color_map, legend_title if legend_title is not None else class_by
                )
            else:
                _draw_curve(ax, result)
                if highlight_features:
                    color_map = _draw_highlights(
                        ax,
                        dataset,
                        result,
                        highlight_features,
                        highlight_groups,
                        highlight_markers,
                        highlight_color,
                        registry_path=registry_path,
                        persist_colors=persist_colors,
                    )
                    if highlight_groups:
                        legend_figure = _legend_figure(
                            color_map,
                            legend_title if legend_title is not None else "highlights",
                            marker_map=dict(highlight_markers)
                            if highlight_markers
                            else None,
                        )

            if log_rank:
                ax.set_xscale("log")
            ax.set_xlabel("abundance rank (1 = most abundant)")
            ax.set_ylabel("log2 abundance (median of detected)")
            ax.grid(True, linestyle=":", alpha=0.4)
            if title is not None:
                fig.suptitle(title, fontsize=14, weight="bold")
        except BaseException:
            plt.close(fig)
            raise

    return DynamicRangePlot(
        figure=fig, legend_figure=legend_figure, result=result, color_map=color_map
    )


def save_dynamic_range(
    dataset: Dataset,
    output_dir: str | Path,
    base_name: str,
    *,
    class_by: str | None = None,
    highlight_features: Mapping[str, str] | None = None,
    highlight_groups: Mapping[str, str] | None = None,
    highlight_markers: Mapping[str, str] | None = None,
    highlight_color: str = DEFAULT_HIGHLIGHT_COLOR,
    min_intensity: float = 0.0,
    log_rank: bool = False,
    title: str | None = None,
    legend_title: str | None = None,
    registry_path: str | Path = DEFAULT_REGISTRY_PATH,
    persist_colors: bool = True,
    dpi: int = 300,
) -> tuple[FigureArtifacts, DynamicRangeResult]:
    """Render a dynamic-range figure and save it; return the artifacts + the result."""
    plot = plot_dynamic_range(
        dataset,
        class_by=class_by,
        highlight_features=highlight_features,
        highlight_groups=highlight_groups,
        highlight_markers=highlight_markers,
        highlight_color=highlight_color,
        min_intensity=min_intensity,
        log_rank=log_rank,
        title=title,
        legend_title=legend_title,
        registry_path=registry_path,
        persist_colors=persist_colors,
    )
    artifacts = save_figure(
        plot.figure, output_dir, base_name, legend_fig=plot.legend_figure, dpi=dpi
    )
    return artifacts, plot.result


# --------------------------------------------------------------------------- #
# Drawing
# --------------------------------------------------------------------------- #


def _draw_curve(ax: Axes, result: DynamicRangeResult) -> None:
    """The whole-cohort median curve + IQR band (2-entry key kept on-axes)."""
    ranks = result.ranks
    ax.fill_between(
        ranks,
        result.log2_q25,
        result.log2_q75,
        color=_BAND_COLOR,
        alpha=0.55,
        label="per-sample IQR",
    )
    ax.plot(
        ranks, result.log2_median, color=_CURVE_COLOR, linewidth=1.7, label="median"
    )
    ax.set_title(
        f"Dynamic range  ({result.dynamic_range_orders:.1f} orders, "
        f"n={result.n_features_detected} features)",
        loc="left",
        fontsize=13,
        weight="bold",
    )
    ax.legend(loc="upper right", frameon=True, fontsize=10)


def _draw_per_class(
    ax: Axes,
    dataset: Dataset,
    class_by: str,
    *,
    min_intensity: float,
    registry_path: str | Path,
    persist_colors: bool,
) -> dict[str, str]:
    """One independently-ranked median curve per class; returns the color map."""
    if class_by not in dataset.metadata.columns:
        raise ValueError(
            f"class_by {class_by!r} is not a metadata column "
            f"{list(dataset.metadata.columns)}."
        )
    series = dataset.metadata[class_by]
    if bool(series.isna().any()):
        raise ValueError(
            f"class_by {class_by!r} has {int(series.isna().sum())} missing value(s); "
            f"resolve or relabel them first."
        )
    classes = series.to_numpy().astype(str)
    abundances = np.asarray(dataset.abundances, dtype=float)
    ordered = _ordered_unique(classes)
    color_map = assign_colors(
        class_by, ordered, registry_path=registry_path, persist=persist_colors
    )
    for cls in ordered:
        members = classes == cls
        median, _, _, count = _detected_quartiles(abundances[members], min_intensity)
        med = median[count >= 1]
        med.sort()
        ax.plot(
            np.arange(1, med.size + 1),
            np.log2(med[::-1]),
            color=color_map[cls],
            linewidth=2.0,
            label=f"{cls} (n={int(members.sum())})",
        )
    ax.set_title(
        "Dynamic range by sample class", loc="left", fontsize=13, weight="bold"
    )
    return color_map


@dataclass(frozen=True)
class _Highlight:
    """One resolved highlight: its rank position, log2 abundance, label, and color."""

    x: int
    y: float
    text: str
    color: str


def _draw_highlights(
    ax: Axes,
    dataset: Dataset,
    result: DynamicRangeResult,
    highlight_features: Mapping[str, str],
    highlight_groups: Mapping[str, str] | None,
    highlight_markers: Mapping[str, str] | None,
    highlight_color: str,
    *,
    registry_path: str | Path,
    persist_colors: bool,
) -> dict[str, str]:
    """Mark + label the highlighted features; return the group color map if grouped."""
    rank_of = {name: i + 1 for i, name in enumerate(result.feature_names_ranked)}
    known = set(np.asarray(dataset.feature_names).astype(str))
    color_map: dict[str, str] = {}
    if highlight_groups is not None:
        missing = [f for f in highlight_features if f not in highlight_groups]
        if missing:
            raise ValueError(
                f"highlight_groups is missing group(s) for {missing[:5]}; every "
                f"highlighted feature needs a group when highlight_groups is given."
            )
        color_map = assign_colors(
            "highlight_group",
            _ordered_unique(np.asarray(list(highlight_groups.values()))),
            registry_path=registry_path,
            persist=persist_colors,
        )

    items: list[_Highlight] = []
    for fid, label in highlight_features.items():
        if fid not in known:
            raise ValueError(f"highlighted feature {fid!r} is not in the Dataset.")
        if fid not in rank_of:
            raise ValueError(
                f"highlighted feature {fid!r} is never detected, so it has no "
                f"abundance to mark; drop it from highlight_features."
            )
        idx = rank_of[fid] - 1
        color = (
            color_map[highlight_groups[fid]] if highlight_groups else highlight_color
        )
        marker = (
            highlight_markers.get(highlight_groups[fid], "o")
            if (highlight_markers and highlight_groups)
            else "o"
        )
        ax.scatter(
            [rank_of[fid]],
            [result.log2_median[idx]],
            s=70,
            color=color,
            marker=marker,
            edgecolor="black",
            linewidth=0.6,
            zorder=3,
        )
        items.append(
            _Highlight(
                x=rank_of[fid],
                y=float(result.log2_median[idx]),
                text=label,
                color=color,
            )
        )

    finite = result.log2_median[np.isfinite(result.log2_median)]
    _place_labels(
        ax, items, result.n_features_detected, float(finite.max()), float(finite.min())
    )
    return color_map


def _place_labels(
    ax: Axes, items: list[_Highlight], n: int, y_top: float, y_bot: float
) -> None:
    """Leader-line labels, de-collided and clamped inside the axes."""
    ax.set_ylim(y_bot - 1.5, y_top + 1.0)
    gap = (y_top - y_bot) * 0.052

    head = sorted((h for h in items if h.x < n * 0.5), key=lambda h: h.y, reverse=True)
    col_x = n * 0.14
    cur = y_top
    for h in head:
        ly = min(h.y, cur)
        cur = ly - gap
        _annotate(ax, h, col_x, ly, "left")

    tail = sorted((h for h in items if h.x >= n * 0.5), key=lambda h: h.y, reverse=True)
    cur = min((h.y for h in tail), default=y_bot) - gap * 1.6
    for h in tail:
        ly = max(cur, y_bot + gap)
        cur = ly - gap
        _annotate(ax, h, float(h.x - n * 0.01), ly, "right")


def _annotate(
    ax: Axes, item: _Highlight, label_x: float, label_y: float, ha: str
) -> None:
    """One leader-line annotation from its label position back to the point."""
    ax.annotate(
        item.text,
        xy=(item.x, item.y),
        xytext=(label_x, label_y),
        fontsize=9,
        color=item.color,
        va="center",
        ha=ha,
        arrowprops={"arrowstyle": "-", "color": item.color, "lw": 0.7},
    )


def _ordered_unique(values: np.ndarray) -> list[str]:
    """Distinct values in first-seen order (stable, study-agnostic)."""
    seen: list[str] = []
    for value in values.astype(str):
        if value not in seen:
            seen.append(value)
    return seen


# --------------------------------------------------------------------------- #
# Legend figure (rendered separately so it never overlaps the plot)
# --------------------------------------------------------------------------- #


def _legend_figure(
    color_map: dict[str, str],
    legend_title: str,
    marker_map: dict[str, str] | None = None,
) -> Figure:
    """Standalone swatch legend: one entry per class / highlight group.

    When ``marker_map`` is given, each entry is a marker glyph (color + shape) matching
    the plotted points; otherwise it is a colored line (the per-class curve key).
    """
    labels = list(color_map.keys())
    height = max(1.4, 0.35 * len(labels) + 0.8)
    width = 3.4 if marker_map is None else max(3.4, 0.11 * max(len(k) for k in labels) + 2.0)
    fig, ax = plt.subplots(figsize=(width, height))
    ax.axis("off")
    if marker_map is None:
        handles = [Line2D([0], [0], color=color_map[k], linewidth=2.5) for k in labels]
    else:
        handles = [
            Line2D(
                [0], [0],
                color=color_map[k],
                marker=marker_map.get(k, "o"),
                markeredgecolor="black",
                markeredgewidth=0.6,
                markersize=10,
                linestyle="None",
            )
            for k in labels
        ]
    ax.legend(
        handles,
        labels,
        title=legend_title or None,
        loc="center",
        frameon=True,
        fontsize=11,
        title_fontsize=12,
    )
    return fig


# --------------------------------------------------------------------------- #
# CLI runner: Stage-3 QC whole-cohort render for the delivered-linear state
# --------------------------------------------------------------------------- #

DEFAULT_STATE_DIR = Path("results") / "qc_states" / "delivered_linear"
DEFAULT_OUTPUT_DIR = Path("figures") / "qc" / "dynamic-range"
DEFAULT_BASE_NAME = "dynamic-range-delivered-linear"
DEFAULT_REGISTRY = Path("state") / "color_registry.json"

# The delivered matrix is normalized + batch-corrected (on Cohort), then de-logged to
# linear (scripts/scratch/qc_prep.py, finding 0004) -- it is NOT raw instrument
# intensity. The title says so honestly.
_STATE_LABEL = "delivered (normalized, batch-corrected) matrix, de-logged to linear"


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(name)s: %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state-dir", type=Path, default=DEFAULT_STATE_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--base-name", type=str, default=DEFAULT_BASE_NAME)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    args = parser.parse_args()

    dataset = load_dataset(args.state_dir)
    n_samples = int(dataset.abundances.shape[0])
    LOGGER.info(
        "loaded %s: %d samples x %d features, scale=%s",
        args.state_dir, n_samples, int(dataset.abundances.shape[1]), dataset.scale,
    )

    title = (
        f"Dynamic range (rank-abundance) -- {_STATE_LABEL}\n"
        f"{n_samples} experimental samples (median + per-sample IQR); "
        f"processed data, not raw intensity"
    )

    # QC whole-cohort render: no highlights (named proteins-of-interest annotation is a
    # downstream Stage-4 concern, per the dispatch); persist_colors=False since the
    # whole-cohort mode uses fixed non-registry curve/band colors.
    artifacts, result = save_dynamic_range(
        dataset,
        args.output_dir,
        args.base_name,
        highlight_features=None,
        min_intensity=0.0,
        title=title,
        registry_path=args.registry,
        persist_colors=False,
    )

    top = result.feature_names_ranked[:10]
    LOGGER.info("dynamic range: %.2f orders of magnitude (log10)",
                result.dynamic_range_orders)
    LOGGER.info("features detected: %d / %d",
                result.n_features_detected, result.n_features_total)
    LOGGER.info("top 10 most abundant: %s", ", ".join(str(t) for t in top))
    LOGGER.info("wrote %s", artifacts.svg)
    LOGGER.info("wrote %s", artifacts.png)
    LOGGER.info("legend image: %s (on-axes key: median / per-sample IQR)",
                artifacts.legend_png)


if __name__ == "__main__":
    main()
