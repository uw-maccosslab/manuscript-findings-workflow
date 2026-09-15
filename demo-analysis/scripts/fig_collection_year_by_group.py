"""Figure: plasma collection year distributed across diagnostic groups.

Stage-1 metadata figure (family: distributions). Documents a covariate skew —
plasma collection year travels with diagnostic group (a secondary confound
alongside source cohort). ADD plasma was collected earlier than the other
groups, and PDD carries a wide spread of late collection years.

Produces a strip plot (jittered points, exact year on y) over a light box
summary, one panel, groups on x in fixed order. Colors come from the project
``group`` dimension of ``state/color_registry.json`` — never invented here.

Dual export: <stem>.svg + <stem>.png (300 DPI), plus a separate legend image
(<stem>.legend.svg / .png). The plot carries no baked-in legend.

Run:
    ./.venv/bin/python scripts/fig_collection_year_by_group.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Final

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.patches import Patch
from scipy import stats

# --- Fixed configuration -----------------------------------------------------

EXPERIMENTAL_GROUPS: Final[tuple[str, ...]] = ("ADD", "PDD", "PDCN", "HCN")
COLOR_DIMENSION: Final[str] = "group"
JITTER_HALF_WIDTH: Final[float] = 0.18
RNG_SEED: Final[int] = 20260915

DEFAULT_DATA: Final[Path] = Path("data/metadata.csv")
DEFAULT_REGISTRY: Final[Path] = Path("state/color_registry.json")
DEFAULT_OUTDIR: Final[Path] = Path("figures/metadata/distributions")
DEFAULT_STEM: Final[str] = "collection-year-by-group"

# --- Publication defaults ----------------------------------------------------


def apply_pub_style() -> None:
    """Legible-at-print-scale defaults; no chartjunk."""
    plt.rcParams.update(
        {
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 12,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 11,
            "ytick.labelsize": 10,
            "legend.frameon": False,
            "svg.fonttype": "none",
        }
    )


# --- Data --------------------------------------------------------------------


def load_group_colors(registry_path: Path, dimension: str) -> dict[str, str]:
    """Read the category->color map for one dimension; never invent colors."""
    registry = json.loads(registry_path.read_text())
    values = registry[dimension]["values"]
    return {str(k): str(v) for k, v in values.items()}


def load_years(data_path: Path, groups: tuple[str, ...]) -> pd.DataFrame:
    """Load metadata, keep experimental groups, parse collection year."""
    df = pd.read_csv(data_path)
    df = df[df["Condition"].isin(groups)].copy()
    parsed = pd.to_datetime(df["Collection date"], format="%m/%d/%Y", errors="raise")
    df["year"] = parsed.dt.year.astype(int)
    return df[["Condition", "year"]].reset_index(drop=True)


def kruskal_year_p(df: pd.DataFrame, groups: tuple[str, ...]) -> float:
    """Kruskal-Wallis p for collection year across the groups."""
    samples = [df.loc[df["Condition"] == g, "year"].to_numpy() for g in groups]
    result = stats.kruskal(*samples)
    return float(result.pvalue)


# --- Plot --------------------------------------------------------------------


def build_plot(
    df: pd.DataFrame,
    groups: tuple[str, ...],
    colors: dict[str, str],
    kw_p: float,
) -> Figure:
    """Strip plot (jittered points, exact year) over a light box summary."""
    rng = np.random.default_rng(RNG_SEED)
    fig, ax = plt.subplots(figsize=(6.4, 4.6))

    positions = list(range(len(groups)))
    per_group_years = [df.loc[df["Condition"] == g, "year"].to_numpy() for g in groups]

    # Light neutral box summary behind the points (median + IQR context).
    bp = ax.boxplot(
        per_group_years,
        positions=positions,
        widths=0.55,
        showfliers=False,
        patch_artist=True,
        medianprops={"color": "#333333", "linewidth": 1.6},
        whiskerprops={"color": "#9a9a9a", "linewidth": 1.0},
        capprops={"color": "#9a9a9a", "linewidth": 1.0},
        boxprops={"edgecolor": "#9a9a9a", "linewidth": 1.0},
        zorder=1,
    )
    for box in bp["boxes"]:
        box.set(facecolor="#f2f2f2", alpha=0.9)

    # Jittered points, exact integer year on y (no vertical distortion).
    for pos, group in enumerate(groups):
        years = df.loc[df["Condition"] == group, "year"].to_numpy()
        x = pos + rng.uniform(-JITTER_HALF_WIDTH, JITTER_HALF_WIDTH, size=years.size)
        ax.scatter(
            x,
            years,
            s=52,
            color=colors[group],
            edgecolor="white",
            linewidth=0.6,
            alpha=0.95,
            zorder=3,
        )

    ax.set_xticks(positions)
    ax.set_xticklabels(list(groups))
    ax.set_xlabel("Diagnostic group")
    ax.set_ylabel("Plasma collection year")
    ax.set_title("Plasma collection year by diagnostic group")

    all_years = df["year"].to_numpy()
    ax.set_ylim(all_years.min() - 0.8, all_years.max() + 0.8)
    ax.set_yticks(range(int(all_years.min()), int(all_years.max()) + 1))
    ax.margins(x=0.06)

    en_dash = "\N{EN DASH}"
    ax.annotate(
        f"Kruskal{en_dash}Wallis p = {kw_p:.3f}  (n = 10/group)",
        xy=(0.99, 0.01),
        xycoords="axes fraction",
        ha="right",
        va="bottom",
        fontsize=9,
        color="#555555",
    )

    fig.tight_layout()
    return fig


def build_legend(groups: tuple[str, ...], colors: dict[str, str]) -> Figure:
    """Legend as its own image — a swatch key, never overlapping the plot."""
    fig, ax = plt.subplots(figsize=(2.2, 1.5))
    ax.axis("off")
    handles = [
        Patch(facecolor=colors[g], edgecolor="white", label=g) for g in groups
    ]
    ax.legend(
        handles=handles,
        title="Diagnostic group",
        loc="center",
        frameon=False,
        handlelength=1.2,
        borderpad=0.4,
        labelspacing=0.7,
    )
    fig.tight_layout()
    return fig


# --- Export ------------------------------------------------------------------


def export(fig: Figure, outdir: Path, stem: str) -> tuple[Path, Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    svg = outdir / f"{stem}.svg"
    png = outdir / f"{stem}.png"
    fig.savefig(svg)
    fig.savefig(png, dpi=300)
    return svg, png


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--stem", type=str, default=DEFAULT_STEM)
    args = parser.parse_args()

    apply_pub_style()

    colors = load_group_colors(args.registry, COLOR_DIMENSION)
    missing = [g for g in EXPERIMENTAL_GROUPS if g not in colors]
    if missing:
        raise SystemExit(f"registry '{COLOR_DIMENSION}' missing colors for: {missing}")
    # 4 categories, well under the 8-color ceiling — direct color encoding is fine.

    df = load_years(args.data, EXPERIMENTAL_GROUPS)
    kw_p = kruskal_year_p(df, EXPERIMENTAL_GROUPS)

    fig = build_plot(df, EXPERIMENTAL_GROUPS, colors, kw_p)
    svg, png = export(fig, args.outdir, args.stem)
    plt.close(fig)

    legend_fig = build_legend(EXPERIMENTAL_GROUPS, colors)
    legend_svg, legend_png = export(legend_fig, args.outdir, f"{args.stem}.legend")
    plt.close(legend_fig)

    print(f"Kruskal-Wallis p (collection year across groups) = {kw_p:.4f}")
    for path in (svg, png, legend_svg, legend_png):
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
