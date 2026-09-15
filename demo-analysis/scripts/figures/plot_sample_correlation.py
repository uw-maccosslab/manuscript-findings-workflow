"""Clustered sample x sample correlation heatmap for the delivered QC state.

Seeded from the ``sample-correlation`` lib/ template
(``lib/figures/correlation.py`` @ 0.1) and adapted to this project: the two shared
helpers are the project's ``common.colors`` / ``common.figure_io`` (not the template's
``figures.*``), and a study ``main()`` loads the materialized Stage-3 QC state and wires
the three design factors to their registry color namespaces. The plotting body is the
vetted template body; it is held to the correctness charter (assume nothing, verify
everything, fail loud).

What it draws (one figure): each cell is the Pearson correlation between two samples'
whole 2340-protein log2 profiles, hierarchically clustered (scipy average linkage), with
three annotation stripes above it -- ``group`` (diagnostic Condition), ``sex``, and
``cohort`` -- so you can read whether the clustering tracks the design. This is the
Stage-3 design sanity check for the delivered matrix (log2, median-normalized,
batch-corrected on cohort); controls are absent (finding 0003), so all 40 samples are
experimental.

Registry wiring (conventions/visualization.md): the metadata columns ``Condition`` /
``Sex`` / ``Cohort`` are mapped onto the registry namespaces ``group`` / ``sex`` /
``cohort`` (with ``Sex`` lower-cased to the universal ``male``/``female`` keys) so every
value keeps the color it has in every other figure. Colors are read, never invented
(``persist=False`` -- a render must not mutate the shared registry). The correlation
colorbar stays on the heatmap; the three annotation keys are a separate legend image.

Scale note: Pearson on a log2 matrix is the correct choice (the template WARNS only when
Pearson runs on a non-log scale, where abundant proteins dominate); ``delivered_log`` is
``scale="log2"``, so no warning fires.

Run:  ./.venv/bin/python scripts/figures/plot_sample_correlation.py
"""

from __future__ import annotations

import argparse
import logging
import sys
import warnings
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, cast

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, to_rgba
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import dendrogram, linkage

# Import the project's shared helpers (seeded from the lib/ templates).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.colors import DEFAULT_REGISTRY_PATH, assign_colors
from common.data_loading import LOG_SCALES, Dataset
from common.dataset_io import load_dataset
from common.figure_io import FigureArtifacts, publication_style, save_figure

LOGGER = logging.getLogger("plot_sample_correlation")

__script_meta__: dict[str, object] = {
    "task": "figure",
    "kind": "figure",
    "provides": [
        "compute_correlation",
        "plot_sample_correlation",
        "save_sample_correlation",
        "correlation_outliers",
    ],
    "uses": [
        "common.data_loading",
        "common.dataset_io",
        "common.colors",
        "common.figure_io",
    ],
    "seeded_from": {"template": "sample-correlation", "version": "0.1"},
    "description": (
        "Clustered sample x sample Pearson correlation heatmap for the delivered QC "
        "state (results/qc_states/delivered_log, log2), with group/sex/cohort registry-"
        "colored annotation stripes; flags globally-low-correlation outlier samples by "
        "robust-Z of mean off-diagonal r; separate legend image; dual-export."
    ),
}

# --------------------------------------------------------------------------- #
# Template body (correlation compute + plot). Unmodified except imports.
# --------------------------------------------------------------------------- #

CorrelationMethod = str
_METHODS: frozenset[str] = frozenset({"pearson", "spearman"})

_HEATMAP_CMAP = "RdYlBu_r"
_VMAX = 1.0
_LINKAGE_METRIC = "euclidean"
_LINKAGE_METHOD = "average"
_METHOD_SYMBOL: dict[str, str] = {"pearson": "Pearson r", "spearman": "Spearman rho"}
_CONTINUOUS_CMAP = "viridis"
_BAD_COLOR = "#bdbdbd"


class CorrelationScaleWarning(UserWarning):
    """Warning that Pearson correlation is running on a non-log (e.g. linear) scale."""


@dataclass(frozen=True)
class _Annotation:
    """One resolved annotation dimension to draw as a stripe above the heatmap."""

    name: str
    continuous: bool
    per_sample: np.ndarray


@dataclass(frozen=True)
class CorrelationResult:
    """The sample-vs-sample correlation underlying the figure (input sample order)."""

    matrix: np.ndarray
    sample_ids: np.ndarray
    method: str
    order: np.ndarray


@dataclass
class CorrelationPlot:
    """A rendered correlation figure plus its companion legend figure."""

    figure: Figure
    legend_figure: Figure | None
    result: CorrelationResult
    color_maps: dict[str, dict[str, str]]


def compute_correlation(
    dataset: Dataset,
    *,
    method: CorrelationMethod = "pearson",
    sample_label_column: str | None = None,
    linkage_metric: str = _LINKAGE_METRIC,
    linkage_method: str = _LINKAGE_METHOD,
    cluster: bool = True,
) -> CorrelationResult:
    """Correlate every pair of samples over their full feature profiles (fail-loud)."""
    if method not in _METHODS:
        raise ValueError(
            f"method {method!r} is not supported; use one of {sorted(_METHODS)}."
        )
    abundances = np.asarray(dataset.abundances, dtype=float)
    if abundances.ndim != 2:
        raise ValueError(
            f"abundances must be 2D (n_samples, n_features); got {abundances.shape}."
        )
    n_samples = abundances.shape[0]
    if n_samples < 2:
        raise ValueError(f"correlation needs >= 2 samples to compare; got {n_samples}.")
    if not np.isfinite(abundances).all():
        n_bad = int((~np.isfinite(abundances)).sum())
        raise ValueError(
            f"compute_correlation requires finite abundances but found {n_bad} "
            f"non-finite value(s) (NaN/Inf). Resolve missingness first."
        )

    sample_ids = dataset.metadata.index.to_numpy().astype(str)
    if len(sample_ids) != n_samples:
        raise ValueError(
            f"metadata has {len(sample_ids)} rows but abundances has {n_samples} "
            f"samples; the Dataset is not row-aligned."
        )
    duplicate_ids = sample_ids[pd.Index(sample_ids).duplicated(keep=False)]
    if duplicate_ids.size:
        raise ValueError(
            f"sample identifiers (metadata index) must be unique; duplicates: "
            f"{sorted(set(duplicate_ids.tolist()))}."
        )
    if sample_label_column is not None and sample_label_column not in (
        dataset.metadata.columns
    ):
        raise ValueError(
            f"sample_label_column {sample_label_column!r} is not a metadata column "
            f"{list(dataset.metadata.columns)}."
        )

    constant = np.nanstd(abundances, axis=1) == 0.0
    if bool(constant.any()):
        bad = sample_ids[constant].tolist()
        raise ValueError(
            f"{len(bad)} sample(s) are constant across all features, so their "
            f"correlation is undefined: {bad[:5]}. Drop or fix them before correlating."
        )

    frame = pd.DataFrame(abundances.T, columns=pd.Index(sample_ids))
    corr_method = cast(Literal["pearson", "kendall", "spearman"], method)
    matrix = frame.corr(method=corr_method).to_numpy(dtype=float)

    if cluster:
        order = _linkage_order(matrix, linkage_metric, linkage_method)
    else:
        order = np.arange(n_samples)
    return CorrelationResult(
        matrix=matrix, sample_ids=sample_ids, method=method, order=order
    )


def _linkage_order(matrix: np.ndarray, metric: str, method: str) -> np.ndarray:
    linkage_z = linkage(matrix, metric=metric, method=method)
    leaves = dendrogram(linkage_z, no_plot=True)["leaves"]
    return np.asarray(leaves, dtype=int)


def _linkage_matrix(matrix: np.ndarray, metric: str, method: str) -> np.ndarray:
    return np.asarray(linkage(matrix, metric=metric, method=method), dtype=float)


def plot_sample_correlation(
    dataset: Dataset,
    *,
    categorical_annotations: Sequence[str] = (),
    continuous_annotations: Sequence[str] = (),
    method: CorrelationMethod = "pearson",
    cluster: bool = True,
    sample_label_column: str | None = None,
    feature_type: str = "feature",
    continuous_colormap: str = _CONTINUOUS_CMAP,
    linkage_metric: str = _LINKAGE_METRIC,
    linkage_method: str = _LINKAGE_METHOD,
    show_tick_labels: bool = True,
    title: str | None = None,
    registry_path: str | Path = DEFAULT_REGISTRY_PATH,
    persist_colors: bool = True,
) -> CorrelationPlot:
    """Render a clustered sample-correlation heatmap with annotation stripes."""
    annotations = _resolve_annotations(
        dataset, categorical_annotations, continuous_annotations
    )
    result = compute_correlation(
        dataset,
        method=method,
        sample_label_column=sample_label_column,
        linkage_metric=linkage_metric,
        linkage_method=linkage_method,
        cluster=cluster,
    )
    if method == "pearson" and dataset.scale not in LOG_SCALES:
        warnings.warn(
            f"Pearson correlation is running on scale {dataset.scale!r}, which is not "
            f"log-ish {sorted(LOG_SCALES)}. On a linear scale a few very abundant "
            f"{feature_type}s dominate, so r can read misleadingly low for "
            f"biologically-similar samples; normalize + log2_transform (or VSN) first, "
            f"or use method='spearman' (rank-based, scale-invariant).",
            CorrelationScaleWarning,
            stacklevel=2,
        )

    order = result.order
    matrix_ordered = result.matrix[np.ix_(order, order)]
    display_labels = _display_labels(dataset, sample_label_column)[order]

    linkage_z = (
        _linkage_matrix(result.matrix, linkage_metric, linkage_method)
        if cluster
        else None
    )

    color_maps: dict[str, dict[str, str]] = {}
    with publication_style():
        fig = plt.figure(figsize=_figsize(len(annotations)))
        try:
            axes = _build_layout(fig, n_annotations=len(annotations), cluster=cluster)
            if title is not None:
                fig.suptitle(title, fontsize=15, weight="bold", y=0.995)

            if linkage_z is not None and axes.dendrogram is not None:
                _draw_dendrogram(axes.dendrogram, linkage_z, n_samples=len(order))

            color_maps = _draw_annotation_stripes(
                axes.stripes,
                annotations,
                order=order,
                continuous_colormap=continuous_colormap,
                registry_path=registry_path,
                persist_colors=persist_colors,
            )

            _draw_heatmap(
                axes.heatmap,
                axes.colorbar,
                matrix_ordered,
                display_labels=display_labels,
                method=result.method,
                feature_type=feature_type,
                show_tick_labels=show_tick_labels,
            )

            legend_figure = (
                _legend_figure(annotations, color_maps, continuous_colormap)
                if annotations
                else None
            )
        except BaseException:
            plt.close(fig)
            raise

    return CorrelationPlot(
        figure=fig,
        legend_figure=legend_figure,
        result=result,
        color_maps=color_maps,
    )


def save_sample_correlation(
    dataset: Dataset,
    output_dir: str | Path,
    base_name: str,
    *,
    categorical_annotations: Sequence[str] = (),
    continuous_annotations: Sequence[str] = (),
    method: CorrelationMethod = "pearson",
    cluster: bool = True,
    sample_label_column: str | None = None,
    feature_type: str = "feature",
    continuous_colormap: str = _CONTINUOUS_CMAP,
    linkage_metric: str = _LINKAGE_METRIC,
    linkage_method: str = _LINKAGE_METHOD,
    show_tick_labels: bool = True,
    title: str | None = None,
    registry_path: str | Path = DEFAULT_REGISTRY_PATH,
    persist_colors: bool = True,
    dpi: int = 300,
) -> FigureArtifacts:
    """Render a correlation figure and dual-export it (+ its separate legend image)."""
    plot = plot_sample_correlation(
        dataset,
        categorical_annotations=categorical_annotations,
        continuous_annotations=continuous_annotations,
        method=method,
        cluster=cluster,
        sample_label_column=sample_label_column,
        feature_type=feature_type,
        continuous_colormap=continuous_colormap,
        linkage_metric=linkage_metric,
        linkage_method=linkage_method,
        show_tick_labels=show_tick_labels,
        title=title,
        registry_path=registry_path,
        persist_colors=persist_colors,
    )
    return save_figure(
        plot.figure,
        output_dir,
        base_name,
        legend_fig=plot.legend_figure,
        dpi=dpi,
    )


def _resolve_annotations(
    dataset: Dataset,
    categorical: Sequence[str],
    continuous: Sequence[str],
) -> list[_Annotation]:
    """Validate + extract each annotation column; categorical first, then continuous."""
    overlap = set(categorical) & set(continuous)
    if overlap:
        raise ValueError(
            f"annotation column(s) {sorted(overlap)} appear in both "
            f"categorical_annotations and continuous_annotations; choose one."
        )
    resolved: list[_Annotation] = []
    for name in categorical:
        if name not in dataset.metadata.columns:
            raise ValueError(
                f"categorical annotation {name!r} is not a metadata column "
                f"{list(dataset.metadata.columns)}."
            )
        series = dataset.metadata[name]
        if series.isna().any():
            n_bad = int(series.isna().sum())
            raise ValueError(
                f"categorical annotation {name!r} has {n_bad} missing value(s); a "
                f"missing label cannot be colored. Resolve or relabel them first."
            )
        resolved.append(
            _Annotation(
                name=name, continuous=False, per_sample=series.to_numpy().astype(str)
            )
        )
    for name in continuous:
        if name not in dataset.metadata.columns:
            raise ValueError(
                f"continuous annotation {name!r} is not a metadata column "
                f"{list(dataset.metadata.columns)}."
            )
        numeric = pd.to_numeric(dataset.metadata[name], errors="coerce")
        if numeric.isna().any():
            n_bad = int(numeric.isna().sum())
            raise ValueError(
                f"continuous annotation {name!r} has {n_bad} value(s) that are missing "
                f"or non-numeric; cannot map them to a colormap."
            )
        resolved.append(
            _Annotation(
                name=name, continuous=True, per_sample=numeric.to_numpy(dtype=float)
            )
        )
    return resolved


def _display_labels(dataset: Dataset, sample_label_column: str | None) -> np.ndarray:
    if sample_label_column is None:
        labels: np.ndarray = dataset.metadata.index.to_numpy().astype(str)
    else:
        labels = dataset.metadata[sample_label_column].to_numpy().astype(str)
    return labels


@dataclass(frozen=True)
class _Axes:
    """The axes of the correlation layout."""

    heatmap: Axes
    colorbar: Axes
    stripes: list[Axes]
    dendrogram: Axes | None


def _figsize(n_annotations: int) -> tuple[float, float]:
    return (16.0, 14.0 + 0.35 * n_annotations)


def _build_layout(fig: Figure, *, n_annotations: int, cluster: bool) -> _Axes:
    """Stack: [dendrogram?] / [annotation stripes...] / [heatmap], + a side colorbar."""
    has_dendro = cluster
    n_rows = (1 if has_dendro else 0) + n_annotations + 1

    dendro_h = 2.2
    stripe_h = 0.45
    heatmap_h = 16.0
    height_ratios: list[float] = []
    if has_dendro:
        height_ratios.append(dendro_h)
    height_ratios.extend([stripe_h] * n_annotations)
    height_ratios.append(heatmap_h)

    gs = GridSpec(
        n_rows,
        2,
        figure=fig,
        width_ratios=[40, 1],
        height_ratios=height_ratios,
        hspace=0.06,
        wspace=0.02,
        left=0.18,
        right=0.90,
        top=0.95,
        bottom=0.12,
    )

    row = 0
    dendro_ax: Axes | None = None
    if has_dendro:
        dendro_ax = fig.add_subplot(gs[row, 0])
        dendro_ax.axis("off")
        row += 1

    stripe_axes: list[Axes] = []
    for _ in range(n_annotations):
        ax = fig.add_subplot(gs[row, 0])
        stripe_axes.append(ax)
        row += 1

    heatmap_ax = fig.add_subplot(gs[row, 0])
    colorbar_ax = fig.add_subplot(gs[row, 1])
    return _Axes(
        heatmap=heatmap_ax,
        colorbar=colorbar_ax,
        stripes=stripe_axes,
        dendrogram=dendro_ax,
    )


def _draw_dendrogram(ax: Axes, linkage_z: np.ndarray, *, n_samples: int) -> None:
    dendro = dendrogram(linkage_z, no_plot=True)
    for xs, ys in zip(dendro["icoord"], dendro["dcoord"], strict=True):
        x_mapped = [(x - 5.0) / 10.0 for x in xs]
        ax.plot(x_mapped, ys, color="#444444", linewidth=1.0)
    ax.set_xlim(-0.5, n_samples - 0.5)
    ax.set_ylim(bottom=0.0)
    ax.axis("off")


def _draw_annotation_stripes(
    axes: list[Axes],
    annotations: list[_Annotation],
    *,
    order: np.ndarray,
    continuous_colormap: str,
    registry_path: str | Path,
    persist_colors: bool,
) -> dict[str, dict[str, str]]:
    """Draw one colored stripe per annotation; return the categorical color maps."""
    color_maps: dict[str, dict[str, str]] = {}
    for ax, annotation in zip(axes, annotations, strict=True):
        ordered = annotation.per_sample[order]
        if annotation.continuous:
            rgba = _continuous_rgba(ordered, continuous_colormap)
        else:
            color_map = assign_colors(
                annotation.name,
                [str(v) for v in annotation.per_sample],
                registry_path=registry_path,
                persist=persist_colors,
            )
            color_maps[annotation.name] = color_map
            rgba = np.array([to_rgba(color_map[str(v)]) for v in ordered])

        ax.imshow(
            rgba.reshape(1, -1, 4),
            aspect="auto",
            interpolation="none",
            extent=(-0.5, len(ordered) - 0.5, -0.5, 0.5),
        )
        ax.set_xlim(-0.5, len(ordered) - 0.5)
        ax.set_yticks([])
        ax.set_xticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_ylabel(
            annotation.name,
            rotation=0,
            ha="right",
            va="center",
            fontsize=12,
            labelpad=10,
        )
    return color_maps


def _continuous_rgba(values: np.ndarray, colormap: str) -> np.ndarray:
    vmin = float(values.min())
    vmax = float(values.max())
    norm = Normalize(vmin=vmin, vmax=vmax if vmax > vmin else vmin + 1.0)
    cmap = plt.get_cmap(colormap)
    return np.asarray(cmap(norm(values)), dtype=float)


def _draw_heatmap(
    ax: Axes,
    colorbar_ax: Axes,
    matrix: np.ndarray,
    *,
    display_labels: np.ndarray,
    method: str,
    feature_type: str,
    show_tick_labels: bool,
) -> None:
    """Draw the correlation heatmap (off-diagonal-min vmin) + its colorbar."""
    off_diagonal = matrix[~np.eye(len(matrix), dtype=bool)]
    finite_off = off_diagonal[np.isfinite(off_diagonal)]
    vmin = float(finite_off.min()) if finite_off.size else 0.0

    cmap = plt.get_cmap(_HEATMAP_CMAP).with_extremes(bad=_BAD_COLOR)
    image = ax.imshow(
        matrix,
        cmap=cmap,
        vmin=vmin,
        vmax=_VMAX,
        aspect="auto",
        interpolation="none",
    )

    n = len(matrix)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    if show_tick_labels:
        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels(display_labels, rotation=90, fontsize=5)
        ax.set_yticklabels(display_labels, fontsize=5)
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    ax.set_xlabel(f"sample ({feature_type} profile)", fontsize=13)

    cbar = ax.figure.colorbar(image, cax=colorbar_ax)
    cbar.set_label(_METHOD_SYMBOL.get(method, method), fontsize=13)


def _legend_figure(
    annotations: list[_Annotation],
    color_maps: dict[str, dict[str, str]],
    continuous_colormap: str,
) -> Figure:
    """A standalone legend: a swatch block per categorical (colorbar per continuous)."""
    n = len(annotations)
    fig, raw_axes = plt.subplots(n, 1, figsize=(4.0, max(1.6, 1.3 * n)))
    axes_list = list(np.atleast_1d(raw_axes))

    for ax, annotation in zip(axes_list, annotations, strict=True):
        if annotation.continuous:
            values = annotation.per_sample.astype(float)
            vmin = float(values.min())
            vmax = float(values.max())
            norm = Normalize(vmin=vmin, vmax=vmax if vmax > vmin else vmin + 1.0)
            mappable = ScalarMappable(cmap=plt.get_cmap(continuous_colormap), norm=norm)
            mappable.set_array(np.asarray([], dtype=float))
            cbar = ax.figure.colorbar(mappable, cax=ax, orientation="horizontal")
            cbar.set_label(annotation.name, fontsize=12)
        else:
            ax.axis("off")
            color_map = color_maps[annotation.name]
            seen: list[str] = []
            for value in annotation.per_sample.astype(str):
                if value not in seen:
                    seen.append(value)
            handles = [
                Line2D(
                    [0],
                    [0],
                    marker="s",
                    linestyle="",
                    markersize=11,
                    markerfacecolor=color_map[value],
                    markeredgecolor="none",
                )
                for value in seen
            ]
            ax.legend(
                handles,
                seen,
                title=annotation.name,
                loc="center",
                frameon=True,
                fontsize=11,
                title_fontsize=12,
                ncol=1 if len(seen) <= 4 else 2,
            )
    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------- #
# Study wiring
# --------------------------------------------------------------------------- #

# Metadata column -> registry color namespace. The registry dims are group / sex /
# cohort (state/color_registry.json); the delivered metadata carries them under the
# raw column names Condition / Sex / Cohort.
_ANNOTATION_NAMESPACES: tuple[tuple[str, str], ...] = (
    ("Condition", "group"),
    ("Sex", "sex"),
    ("Cohort", "cohort"),
)

DEFAULT_STATE_DIR = Path("results") / "qc_states" / "delivered_log"
DEFAULT_OUTPUT_DIR = Path("figures") / "qc" / "sample-correlation"
DEFAULT_BASE_NAME = "sample-correlation-delivered"
DEFAULT_REGISTRY = Path("state") / "color_registry.json"
# Processing-state provenance, surfaced in the on-canvas title (finding 0004).
STATE_DESCRIPTION = "delivered: log2, median-normalized, batch-corrected on cohort"
# Robust-Z threshold on per-sample mean off-diagonal r below which a sample is a clear
# global-correlation outlier (a candidate QC problem worth chasing).
OUTLIER_ROBUST_Z = -3.5


def _prepare_annotated_dataset(dataset: Dataset) -> Dataset:
    """Attach the three registry-namespaced annotation columns to a Dataset copy.

    Maps Condition->group / Sex->sex / Cohort->cohort so the template's per-column
    ``assign_colors`` call resolves each stripe against the correct registry namespace,
    lower-casing Sex to the universal male/female keys. Fails loud on a missing column.
    """
    md = dataset.metadata.copy()
    for source, namespace in _ANNOTATION_NAMESPACES:
        if source not in md.columns:
            raise ValueError(
                f"expected metadata column {source!r} for the {namespace!r} stripe; "
                f"have {list(md.columns)}."
            )
        values = md[source].astype(str)
        if namespace == "sex":
            values = values.str.lower()
        md[namespace] = values
    return Dataset(
        abundances=dataset.abundances,
        feature_names=dataset.feature_names,
        feature_metadata=dataset.feature_metadata,
        metadata=md,
        scale=dataset.scale,
    )


def correlation_outliers(
    result: CorrelationResult, *, robust_z: float = OUTLIER_ROBUST_Z
) -> list[tuple[str, float, float]]:
    """Samples whose mean off-diagonal correlation is a low outlier (robust-Z/MAD).

    Returns ``[(sample_id, mean_off_diagonal_r, robust_z), ...]`` sorted ascending by r,
    for every sample at or below ``robust_z``. An empty list means no sample stands out
    as globally poorly-correlated.
    """
    matrix = result.matrix
    off = matrix.copy()
    np.fill_diagonal(off, np.nan)
    mean_off = np.nanmean(off, axis=1)
    median = float(np.median(mean_off))
    mad = float(np.median(np.abs(mean_off - median)))
    scale = 1.4826 * mad if mad > 0 else 0.0
    flagged: list[tuple[str, float, float]] = []
    for sid, value in zip(result.sample_ids, mean_off, strict=True):
        z = (value - median) / scale if scale > 0 else 0.0
        if z <= robust_z:
            flagged.append((str(sid), float(value), float(z)))
    flagged.sort(key=lambda item: item[1])
    return flagged


def build_title(outliers: Sequence[tuple[str, float, float]], *, n: int) -> str:
    """Two-line suptitle: the comparison + state, then the terse load-bearing n/r."""
    line1 = f"Sample-sample correlation ({STATE_DESCRIPTION})"
    detail = f"Pearson r on log2 protein profiles - n={n} experimental - 2340 proteins"
    if outliers:
        marker = ", ".join(f"{sid} (mean r={r:.2f})" for sid, r, _ in outliers)
        detail += f" - correlation outlier: {marker}"
    return f"{line1}\n{detail}"


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(name)s: %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state-dir", type=Path, default=DEFAULT_STATE_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--base-name", type=str, default=DEFAULT_BASE_NAME)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument(
        "--method", choices=sorted(_METHODS), default="pearson"
    )
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()

    dataset = load_dataset(args.state_dir)
    LOGGER.info(
        "loaded %s: %d samples x %d proteins, scale=%s",
        args.state_dir,
        dataset.abundances.shape[0],
        dataset.abundances.shape[1],
        dataset.scale,
    )
    annotated = _prepare_annotated_dataset(dataset)

    # Detect outliers up front so the SampleID can be surfaced in the title.
    result = compute_correlation(annotated, method=args.method)
    outliers = correlation_outliers(result)
    if outliers:
        for sid, r, z in outliers:
            LOGGER.info(
                "CORRELATION OUTLIER: %s mean off-diagonal r=%.4f (robust-Z=%.2f)",
                sid,
                r,
                z,
            )
    else:
        LOGGER.info("no global-correlation outlier samples (all within robust-Z floor)")

    n_samples = dataset.abundances.shape[0]
    title = build_title(outliers, n=n_samples)

    artifacts = save_sample_correlation(
        annotated,
        args.output_dir,
        args.base_name,
        categorical_annotations=[ns for _, ns in _ANNOTATION_NAMESPACES],
        method=args.method,
        feature_type="protein",
        title=title,
        registry_path=args.registry,
        persist_colors=False,  # colors pre-exist; a render must not mutate the registry
        dpi=args.dpi,
    )
    LOGGER.info("wrote %s", artifacts.svg)
    LOGGER.info("wrote %s", artifacts.png)
    LOGGER.info("wrote %s", artifacts.legend_svg)
    LOGGER.info("wrote %s", artifacts.legend_png)


if __name__ == "__main__":
    main()
