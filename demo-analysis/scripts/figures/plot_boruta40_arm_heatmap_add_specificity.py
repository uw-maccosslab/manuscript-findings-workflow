"""Per-protein arm heatmap of the Boruta-40 ADD signature, with a cohort color bar.

Renders the 40-protein ADD-discriminative signature (Boruta-40 Confirmed,
ADD-vs-Other) across all 40 samples, grouped/blocked by arm (ADD, PDD, PDCN, HCN),
with two annotation bars above the columns (arm + collection cohort). The visual
point (findings 0001/0004): the ADD column-block is coherent (up-block high,
down-block low) AND the cohort bar shows ADD is 100% Kerchner while PDD/PDCN are
Poston -- the ADD signature is ALIASED with the Kerchner site.

Data (load -> subset -> plot; recomputes nothing):
  * Matrix + labels: scripts/scratch/project_loader.py -- load_paired_dataset() then
    load_analysis_subset() -> Dataset (40 samples x 2340 proteins, log2). Arm =
    metadata['Condition']; cohort/site = metadata['Cohort'].
  * The 40 signature proteins + direction: results/four_arm/8d418bb8509f/
    signature_direction.parquet (Boruta source of truth: results/boruta/762e16f294e3/).

Encoding:
  * Heatmap body = ROW-Z-SCORED abundance (per protein, z across the 40 samples), a
    diverging colormap centered at 0 (the one legitimate diverging use -- signed z).
  * Rows = the 40 Boruta-confirmed proteins in two blocks: up-in-ADD (18) then
    down-in-ADD (22), each ordered by |effect|. Row labels = UniProt entry-name
    mnemonic (e.g. CSN3), trimmed.
  * Columns = the 40 samples, blocked by arm in order ADD, PDD, PDCN, HCN; within an
    arm, ordered by cohort then by composite signature score.
  * Annotation bars (colors READ from state/color_registry.json, registry not
    mutated): arm via the 'group' category; cohort via the 'cohort' category.

Run:  ./.venv/bin/python scripts/figures/plot_boruta40_arm_heatmap_add_specificity.py
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

# Reuse the project shared figure-IO + color-registry modules and the Stage-3 loader.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scratch"))
import project_loader as pl
from common.colors import get_color
from common.figure_io import publication_style, save_figure

LOGGER = logging.getLogger("plot_boruta40_arm_heatmap_add_specificity")

DEFAULT_SIGNATURE = (
    Path("results") / "four_arm" / "8d418bb8509f" / "signature_direction.parquet"
)
DEFAULT_OUTPUT_DIR = Path("figures") / "analysis" / "four_arm" / "add-specificity"
DEFAULT_BASE_NAME = "boruta40-arm-heatmap-add-specificity"
DEFAULT_REGISTRY = Path("state") / "color_registry.json"

# Provenance (recorded on the figure for regenerability / staleness tracking).
DATA_VERSION = "protein-data.csv sha256:90a0462eb03a02bd"
FINGERPRINT = "8d418bb8509f"

# Column (sample) block order and row (direction) block order.
ARM_ORDER = ("ADD", "PDD", "PDCN", "HCN")
DIRECTION_ORDER = ("up_in_ADD", "down_in_ADD")
DIRECTION_LABEL = {"up_in_ADD": "Up in ADD", "down_in_ADD": "Down in ADD"}

# Diverging colormap for signed row-z (centered at 0). Colorbrewer RdBu reversed:
# high z -> red, low z -> blue. NOT the qualitative registry palette.
HEATMAP_CMAP = "RdBu_r"


def mnemonic(feature_id: str) -> str:
    """UniProt entry-name mnemonic (e.g. 'sp|Q9UNS2|CSN3_HUMAN' -> 'CSN3').

    Protein groups joined by ' / ' are collapsed to 'FIRST/SECOND(/+n)'.
    """
    parts = [p.strip() for p in str(feature_id).split("/") if p.strip()]
    names: list[str] = []
    for part in parts:
        fields = part.split("|")
        name = fields[2] if len(fields) >= 3 else part
        names.append(name.replace("_HUMAN", ""))
    if not names:
        return str(feature_id)
    if len(names) == 1:
        return names[0]
    label = "/".join(names[:2])
    if len(names) > 2:
        label += f"/+{len(names) - 2}"
    return label


def _load_signature(path: Path) -> pd.DataFrame:
    """Read + validate the 40-protein signature table (fail loud)."""
    sig = pd.read_parquet(path)
    required = {"feature", "effect_add_minus_other", "direction"}
    missing = required - set(sig.columns)
    if missing:
        raise ValueError(f"Signature table {path} missing columns: {sorted(missing)}.")
    if len(sig) != 40:
        raise ValueError(f"Expected 40 signature proteins in {path}; got {len(sig)}.")
    counts = sig["direction"].value_counts().to_dict()
    if counts.get("up_in_ADD") != 18 or counts.get("down_in_ADD") != 22:
        raise ValueError(
            f"Expected 18 up / 22 down in ADD; got {counts} in {path}."
        )
    return sig


def _order_rows(sig: pd.DataFrame) -> pd.DataFrame:
    """Row order: up-in-ADD block (by effect desc) then down-in-ADD (by effect asc)."""
    blocks = []
    for direction in DIRECTION_ORDER:
        block = sig[sig["direction"] == direction].copy()
        ascending = direction == "down_in_ADD"
        block = block.sort_values("effect_add_minus_other", ascending=ascending)
        blocks.append(block)
    return pd.concat(blocks, ignore_index=True)


def _order_samples(
    meta: pd.DataFrame, composite: np.ndarray
) -> tuple[np.ndarray, list[str], list[str]]:
    """Column order: arms in ARM_ORDER; within arm, by cohort then composite desc.

    Returns positional sample indices, and the ordered arm / cohort label lists.
    """
    frame = meta.reset_index(drop=True).copy()
    frame["_pos"] = np.arange(len(frame))
    frame["_composite"] = composite
    arm_rank = {arm: i for i, arm in enumerate(ARM_ORDER)}
    unknown = set(frame["Condition"]) - set(arm_rank)
    if unknown:
        raise ValueError(f"Unexpected arm(s) not in ARM_ORDER: {sorted(unknown)}.")
    frame["_arm_rank"] = frame["Condition"].map(arm_rank)
    frame = frame.sort_values(
        by=["_arm_rank", "Cohort", "_composite"],
        ascending=[True, True, False],
    )
    order = frame["_pos"].to_numpy()
    arms = frame["Condition"].tolist()
    cohorts = frame["Cohort"].tolist()
    return order, arms, cohorts


def _annotation_row_rgba(
    labels: list[str], category: str, registry_path: Path
) -> np.ndarray:
    """(1, N, 4) RGBA strip of registry colors for a categorical annotation bar."""
    colors = [
        mpl.colors.to_rgba(
            get_color(category, lab, registry_path=registry_path, persist=False)
        )
        for lab in labels
    ]
    return np.array(colors).reshape(1, len(labels), 4)


def _arm_boundaries(arms: list[str]) -> list[int]:
    """Column indices where the arm changes (for group separators)."""
    return [i for i in range(1, len(arms)) if arms[i] != arms[i - 1]]


def make_heatmap(
    Z: np.ndarray,
    row_labels: list[str],
    n_up: int,
    arms: list[str],
    cohorts: list[str],
    vlim: float,
    registry_path: Path,
) -> Figure:
    """Build the arm heatmap with arm + cohort annotation bars above the columns."""
    n_rows, n_cols = Z.shape
    fig = Figure(figsize=(11.0, 12.5))
    gs = GridSpec(
        3,
        1,
        figure=fig,
        height_ratios=[0.5, 0.5, 22.0],
        hspace=0.06,
        left=0.16,
        right=0.98,
        top=0.95,
        bottom=0.11,
    )
    ax_arm = fig.add_subplot(gs[0, 0])
    ax_cohort = fig.add_subplot(gs[1, 0])
    ax_hm = fig.add_subplot(gs[2, 0])

    arm_boundaries = _arm_boundaries(arms)

    # --- annotation bars (imshow of RGBA strips) ---
    for ax, labels, category, title in (
        (ax_arm, arms, "group", "Arm"),
        (ax_cohort, cohorts, "cohort", "Cohort"),
    ):
        strip = _annotation_row_rgba(labels, category, registry_path)
        ax.imshow(strip, aspect="auto", extent=(0, n_cols, 0, 1), zorder=1)
        ax.set_xlim(0, n_cols)
        ax.set_ylim(0, 1)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_ylabel(
            title, rotation=0, ha="right", va="center", fontsize=10, labelpad=8
        )
        for b in arm_boundaries:
            ax.axvline(b, color="white", lw=1.5, zorder=2)

    # --- heatmap body ---
    im = ax_hm.imshow(
        Z,
        aspect="auto",
        cmap=HEATMAP_CMAP,
        vmin=-vlim,
        vmax=vlim,
        extent=(0, n_cols, n_rows, 0),
        interpolation="nearest",
    )
    im.set_rasterized(True)
    ax_hm.set_xlim(0, n_cols)
    ax_hm.set_ylim(n_rows, 0)

    # row (protein) labels
    ax_hm.set_yticks(np.arange(n_rows) + 0.5)
    ax_hm.set_yticklabels(row_labels, fontsize=6)
    ax_hm.tick_params(axis="y", length=0)

    # column ticks: none (samples anonymous); label arm blocks below instead
    ax_hm.set_xticks([])

    # arm-block separators (vertical) + up/down block separator (horizontal)
    for b in arm_boundaries:
        ax_hm.axvline(b, color="white", lw=1.5)
    ax_hm.axhline(n_up, color="white", lw=2.5)

    # arm block labels centered under each arm group
    bounds = [0, *arm_boundaries, n_cols]
    for i in range(len(bounds) - 1):
        lo, hi = bounds[i], bounds[i + 1]
        arm = arms[lo]
        ax_hm.text(
            (lo + hi) / 2.0,
            n_rows + 0.9,
            f"{arm}\n(n={hi - lo})",
            ha="center",
            va="top",
            fontsize=10,
            fontweight="bold",
            color=get_color("group", arm, registry_path=registry_path, persist=False),
        )
    ax_hm.set_xlabel("Samples, grouped by arm", labelpad=58)

    # up / down block brackets on the left
    for label, y0, y1 in (
        (DIRECTION_LABEL["up_in_ADD"], 0, n_up),
        (DIRECTION_LABEL["down_in_ADD"], n_up, n_rows),
    ):
        ax_hm.annotate(
            f"{label}\n({y1 - y0})",
            xy=(-0.145, 1 - (y0 + y1) / 2.0 / n_rows),
            xycoords="axes fraction",
            ha="center",
            va="center",
            fontsize=9,
            fontweight="bold",
            rotation=90,
        )

    for spine in ax_hm.spines.values():
        spine.set_visible(False)

    ax_arm.set_title(
        "Boruta-40 ADD signature across arms  ·  row-z abundance  ·  "
        "EXPLORATORY / confound-bounded (ADD aliased with Kerchner site)",
        fontsize=12,
        fontweight="bold",
        pad=10,
    )
    return fig


def make_legend(vlim: float, registry_path: Path) -> Figure:
    """Standalone legend: diverging z colorbar + arm key + cohort key."""
    fig = Figure(figsize=(8.5, 2.2))
    gs = GridSpec(
        1, 3, figure=fig, width_ratios=[1.1, 1.0, 1.2], wspace=0.5,
        left=0.04, right=0.97, top=0.80, bottom=0.28,
    )

    # colorbar
    ax_cb = fig.add_subplot(gs[0, 0])
    norm = mpl.colors.Normalize(vmin=-vlim, vmax=vlim)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=HEATMAP_CMAP)
    cb = fig.colorbar(sm, cax=ax_cb, orientation="horizontal")
    cb.set_label("Row z-score (per-protein, across 40 samples)", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    # arm key
    ax_arm = fig.add_subplot(gs[0, 1])
    _swatch_key(
        ax_arm, "Arm", list(ARM_ORDER), "group", registry_path
    )

    # cohort key (registry cohort values, in a stable order)
    ax_cohort = fig.add_subplot(gs[0, 2])
    _swatch_key(
        ax_cohort,
        "Cohort",
        ["Kerchner", "Poston", "ADRC/PUC", "PUC only"],
        "cohort",
        registry_path,
    )
    return fig


def _swatch_key(
    ax: mpl.axes.Axes,
    title: str,
    labels: list[str],
    category: str,
    registry_path: Path,
) -> None:
    """Draw a titled vertical swatch key on ``ax``."""
    ax.set_title(title, fontsize=10, fontweight="bold", loc="left")
    n = len(labels)
    for i, lab in enumerate(labels):
        y = n - 1 - i
        color = get_color(category, lab, registry_path=registry_path, persist=False)
        ax.add_patch(
            Rectangle((0.0, y + 0.15), 0.6, 0.7, facecolor=color, edgecolor="none")
        )
        ax.text(0.75, y + 0.5, lab, va="center", ha="left", fontsize=9)
    ax.set_xlim(0, 3.0)
    ax.set_ylim(-0.2, n)
    ax.axis("off")


def build(
    signature_path: Path, registry_path: Path
) -> tuple[Figure, Figure]:
    """Load data, assemble the row-z matrix, and render heatmap + legend."""
    dataset = pl.load_analysis_subset(pl.load_paired_dataset())
    if dataset.scale != "log2":
        raise ValueError(f"Expected log2 Dataset; got scale={dataset.scale!r}.")

    sig = _order_rows(_load_signature(signature_path))
    feature_index = {name: i for i, name in enumerate(dataset.feature_names.tolist())}
    missing = [f for f in sig["feature"] if f not in feature_index]
    if missing:
        raise ValueError(
            f"{len(missing)} signature feature(s) absent from the matrix: {missing[:3]}"
        )
    col_idx = np.array([feature_index[f] for f in sig["feature"]], dtype=int)

    # proteins x samples, then row-z across the 40 samples
    M = dataset.abundances[:, col_idx].T.astype(float)  # (40 proteins, 40 samples)
    mean = M.mean(axis=1, keepdims=True)
    std = M.std(axis=1, ddof=0, keepdims=True)
    if np.any(std == 0):
        raise ValueError("A signature protein has zero variance; cannot row-z-score.")
    Z = (M - mean) / std

    n_up = int((sig["direction"] == "up_in_ADD").sum())
    # composite signature score per sample = mean up-block z - mean down-block z
    composite = Z[:n_up].mean(axis=0) - Z[n_up:].mean(axis=0)

    order, arms, cohorts = _order_samples(dataset.metadata, composite)
    Z = Z[:, order]

    row_labels = [mnemonic(f) for f in sig["feature"]]
    vlim = float(np.clip(np.percentile(np.abs(Z), 99), 1.5, 3.0))
    vlim = round(vlim, 1)

    with publication_style():
        fig = make_heatmap(
            Z, row_labels, n_up, arms, cohorts, vlim, registry_path
        )
        legend = make_legend(vlim, registry_path)
    return fig, legend


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--signature", type=Path, default=DEFAULT_SIGNATURE)
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    p.add_argument("--base-name", type=str, default=DEFAULT_BASE_NAME)
    p.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = parse_args(argv)
    fig, legend = build(args.signature, args.registry)
    artifacts = save_figure(
        fig, args.output_dir, args.base_name, legend_fig=legend
    )
    LOGGER.info("wrote %s", artifacts.svg)
    LOGGER.info("wrote %s", artifacts.png)
    LOGGER.info("wrote %s", artifacts.legend_svg)
    LOGGER.info("wrote %s", artifacts.legend_png)
    LOGGER.info(
        "provenance: data_version=%s fingerprint=%s", DATA_VERSION, FINGERPRINT
    )


if __name__ == "__main__":
    main()
