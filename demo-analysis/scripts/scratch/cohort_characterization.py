"""Stage-1 cohort characterization: a visual "Table 1" for the plasma-EV study.

Renders a four-panel figure summarizing the 40 experimental samples (four
diagnostic groups) of the plasma-EV proteomics cohort:

  (a) samples per diagnostic group (bar; documents the n=10 balance),
  (b) sex composition per group (stacked Female/Male bar),
  (c) age (years) by group (box + jittered strip),
  (d) education (years) by group (box + jittered strip; one NA dropped).

Group-colored panels use the registry ``group`` dimension; the sex panel uses
the universal ``sex`` colors. Dual-exports SVG + 300-DPI PNG plus a separate
legend image via the project's shared figure-IO module.

Held to the correctness charter: control references are excluded explicitly,
NA handling is logged (never silent), and expectations are asserted.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.patches import Patch

from scripts.common.colors import get_color
from scripts.common.figure_io import publication_style, save_figure

__script_meta__: dict[str, object] = {
    "template": {"name": "cohort-characterization", "version": "0.1"},
    "kind": "figure",
    "provides": ["cohort-characterization"],
    "uses": ["figure-io", "okabe-ito-colors"],
    "seeded_from": None,
    "description": (
        "Stage-1 cohort 'Table 1' multi-panel figure (counts, sex composition, "
        "age, education) for the four experimental diagnostic groups."
    ),
}

logger = logging.getLogger(__name__)

# The four experimental diagnostic groups, in the requested display order.
GROUP_ORDER: tuple[str, ...] = ("ADD", "PDD", "PDCN", "HCN")
# Control references to exclude from every panel (Stage-1 experimental/control split).
CONTROL_CONDITIONS: frozenset[str] = frozenset(
    {"Plasma Reference", "TPAD plasma pooled reference"}
)
SEX_ORDER: tuple[str, ...] = ("Female", "Male")

DEFAULT_DATA = Path("data/metadata.csv")
DEFAULT_OUTPUT_DIR = Path("figures/metadata/distributions")
DEFAULT_STEM = "cohort-characterization"
JITTER_SEED = 20260915


def load_experimental(data_path: Path) -> pd.DataFrame:
    """Load metadata and return the experimental subset, fail-loud on surprises."""
    raw = pd.read_csv(data_path)
    for col in ("Condition", "Age", "Sex", "Education"):
        if col not in raw.columns:
            raise ValueError(f"metadata is missing required column {col!r}")

    exp = raw[raw["Condition"].isin(GROUP_ORDER)].copy()
    # Guard: nothing that is a control reference may survive the group filter.
    leaked = set(exp["Condition"]).intersection(CONTROL_CONDITIONS)
    if leaked:
        raise ValueError(f"control conditions leaked into experimental data: {leaked}")

    unexpected_sex = set(exp["Sex"]).difference(SEX_ORDER)
    if unexpected_sex:
        raise ValueError(f"unexpected Sex values: {sorted(unexpected_sex)}")

    exp["Age"] = pd.to_numeric(exp["Age"], errors="raise")
    logger.info("experimental samples: %d across groups %s", len(exp), GROUP_ORDER)
    return exp


def _group_colors() -> dict[str, str]:
    return {g: get_color("group", g) for g in GROUP_ORDER}


def _sex_colors() -> dict[str, str]:
    # Registry stores sex keys lowercase; map the capitalized data labels onto them.
    return {label: get_color("sex", label.lower()) for label in SEX_ORDER}


def _panel_counts(ax: plt.Axes, exp: pd.DataFrame, gcolors: dict[str, str]) -> None:
    counts = [int((exp["Condition"] == g).sum()) for g in GROUP_ORDER]
    x = np.arange(len(GROUP_ORDER))
    bars = ax.bar(
        x, counts, color=[gcolors[g] for g in GROUP_ORDER], edgecolor="black", width=0.7
    )
    for rect, n in zip(bars, counts, strict=True):
        ax.text(
            rect.get_x() + rect.get_width() / 2,
            n + 0.15,
            str(n),
            ha="center",
            va="bottom",
            fontsize=10,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(GROUP_ORDER)
    ax.set_ylabel("Samples (n)")
    ax.set_ylim(0, max(counts) + 1.5)
    ax.set_title("(a) Samples per group")


def _panel_sex(ax: plt.Axes, exp: pd.DataFrame, scolors: dict[str, str]) -> None:
    x = np.arange(len(GROUP_ORDER))

    def _sex_counts(sex: str) -> np.ndarray:
        return np.array(
            [
                int(((exp["Condition"] == g) & (exp["Sex"] == sex)).sum())
                for g in GROUP_ORDER
            ]
        )

    female = _sex_counts("Female")
    male = _sex_counts("Male")
    ax.bar(
        x, female, color=scolors["Female"], edgecolor="black", width=0.7, label="Female"
    )
    ax.bar(
        x,
        male,
        bottom=female,
        color=scolors["Male"],
        edgecolor="black",
        width=0.7,
        label="Male",
    )
    ax.set_xticks(x)
    ax.set_xticklabels(GROUP_ORDER)
    ax.set_ylabel("Samples (n)")
    ax.set_ylim(0, int((female + male).max()) + 1.5)
    ax.set_title("(b) Sex composition")


def _box_strip(
    ax: plt.Axes,
    exp: pd.DataFrame,
    value_col: str,
    gcolors: dict[str, str],
    rng: np.random.Generator,
) -> None:
    series_by_group: list[np.ndarray] = []
    for g in GROUP_ORDER:
        vals = pd.to_numeric(exp.loc[exp["Condition"] == g, value_col], errors="coerce")
        n_missing = int(vals.isna().sum())
        if n_missing:
            logger.info(
                "%s: dropping %d NA value(s) for group %s", value_col, n_missing, g
            )
        series_by_group.append(vals.dropna().to_numpy(dtype=float))

    positions = np.arange(1, len(GROUP_ORDER) + 1)
    bp = ax.boxplot(
        series_by_group,
        positions=positions,
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.5},
    )
    for patch, g in zip(bp["boxes"], GROUP_ORDER, strict=True):
        patch.set_facecolor(gcolors[g])
        patch.set_alpha(0.45)
        patch.set_edgecolor("black")
    for pos, vals, g in zip(positions, series_by_group, GROUP_ORDER, strict=True):
        jitter = rng.uniform(-0.16, 0.16, size=vals.shape[0])
        ax.scatter(
            np.full(vals.shape[0], pos) + jitter,
            vals,
            color=gcolors[g],
            edgecolor="black",
            linewidth=0.5,
            s=28,
            zorder=3,
        )
    ax.set_xticks(positions)
    ax.set_xticklabels(GROUP_ORDER)


def build_figure(exp: pd.DataFrame) -> Figure:
    gcolors = _group_colors()
    scolors = _sex_colors()
    rng = np.random.default_rng(JITTER_SEED)

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 8.0))
    _panel_counts(axes[0, 0], exp, gcolors)
    _panel_sex(axes[0, 1], exp, scolors)

    _box_strip(axes[1, 0], exp, "Age", gcolors, rng)
    axes[1, 0].set_ylabel("Age (years)")
    axes[1, 0].set_title("(c) Age by group")

    _box_strip(axes[1, 1], exp, "Education", gcolors, rng)
    axes[1, 1].set_ylabel("Education (years)")
    axes[1, 1].set_title("(d) Education by group")

    fig.suptitle(
        "Cohort characterization — plasma-EV proteomics (n=40 experimental)",
        fontsize=15,
        fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig


def build_legend_figure() -> Figure:
    gcolors = _group_colors()
    scolors = _sex_colors()
    fig, (ax_g, ax_s) = plt.subplots(1, 2, figsize=(5.2, 1.7))
    for ax in (ax_g, ax_s):
        ax.axis("off")
    g_handles = [
        Patch(facecolor=gcolors[g], edgecolor="black", label=g) for g in GROUP_ORDER
    ]
    ax_g.legend(
        handles=g_handles,
        title="Diagnostic group",
        loc="center",
        frameon=False,
        handlelength=1.2,
    )
    s_handles = [
        Patch(facecolor=scolors[s], edgecolor="black", label=s) for s in SEX_ORDER
    ]
    ax_s.legend(
        handles=s_handles,
        title="Sex",
        loc="center",
        frameon=False,
        handlelength=1.2,
    )
    return fig


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--stem", type=str, default=DEFAULT_STEM)
    args = parser.parse_args()

    exp = load_experimental(args.data)
    with publication_style():
        fig = build_figure(exp)
        legend_fig = build_legend_figure()
    artifacts = save_figure(fig, args.output_dir, args.stem, legend_fig=legend_fig)
    logger.info("wrote %s", artifacts.png)
    logger.info("wrote %s", artifacts.legend_png)


if __name__ == "__main__":
    main()
