"""Boruta all-relevant importance figure for the PRIMARY ADD-vs-Other selection.

Renders the CACHED Boruta result (it recomputes NOTHING — conventions/results-cache.md;
a *load -> plot* job):

  * ``results/boruta/762e16f294e3/`` — Boruta all-relevant selection, ADD vs Other,
    random-forest importance vs its built-in permuted-shadow noise floor on the whole
    40 x 2340 delivered log2 matrix (no CV). 40 Confirmed / 15 Tentative /
    2285 Rejected; median max-shadow threshold 0.0053; converged in 99 iterations.

Produced by ``scripts/scratch/boruta_add_vs_other.py`` (see ``results/manifest.md``).
Seeded from ``lib/figures/boruta_importance`` (``plot_boruta_importance``) — this reuses
that template's horizontal per-iteration box-plot skeleton (features ordered by Boruta
ranking; a shadow-noise-floor reference line) but DIVERGES from its look in two ways the
dispatch asked for:

  1. **Colored by DE direction, not viridis.** The template colors accepted features by
     median importance on a continuous ``viridis`` scale. Here the 40 Confirmed are
     colored by their **differential-abundance direction** — up-in-ADD vs down-in-ADD —
     from the moderated DE contrast (``results/differential_abundance/add-vs-other-
     moderated/contrast_table.csv``; effect sign of the ``ADD vs Other`` term), via
     the ``primary_contrast`` registry colors (ADD = up, Other = down). Tentative and
     Rejected are the Boruta **decision** states (house-style, not metadata categories —
     template treats decision colors the same way). Because the coloring is now
     **categorical**, this figure ships a **separate legend image** (the template's
     on-axes / colorbar exception applied only to its continuous viridis key).
  2. **Zoomed to the selected head.** Of 2340 features only the informative head is
     drawn legibly — all 40 Confirmed + all 15 Tentative + the top ``top_rejected``
     Rejected (default 10) as an anchor at the noise floor, with gene-symbol labels;
     the remaining ~2275 Rejected are summarized in the title / legend, not drawn.

Encoding (colors read from ``state/color_registry.json``; registry read, not mutated):
  * Confirmed, up in ADD (18)   -> primary_contrast ADD  (#D55E00)
  * Confirmed, down in ADD (22) -> primary_contrast Other (#999999)
  * Tentative (15)              -> house-style Boruta yellow
  * Rejected (top 10 of 2285)   -> house-style light gray
  * max-shadow noise floor (median over 99 iters) = 0.0053 -> dashed line
  * shadow-max IQR (per-iteration spread of the noise floor) -> shaded band

LOAD-BEARING HONESTY (caption; findings 0001/0004): the shadow (permuted-copy) features
are the BUILT-IN null (relevance vs noise) — one pass over the whole 40x2340 matrix, no
CV. EXPLORATORY (n=10 ADD positives). HARD CONFOUND (UNRESOLVED): ADD is 100% Kerchner
cohort (finding 0001), so "all-relevant to ADD" is all-relevant to the ADD-vs-site
partition, not necessarily to disease; the shadow null does not resolve this. The full
narrative lives in the finding's caption, not on the canvas — the figure carries only a
terse caveat marker.

Run:  ./.venv/bin/python scripts/figures/plot_boruta_importance_add_vs_other.py
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.patches import Patch

# Reuse the project shared figure-IO + color-registry modules, and the scratch analysis
# module that defines the cached BorutaResult dataclass (type-directed reload).
_SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_SCRIPTS))
sys.path.insert(0, str(_SCRIPTS / "scratch"))
from boruta_selection import BorutaResult  # noqa: E402
from common.colors import get_color  # noqa: E402
from common.figure_io import publication_style, save_figure  # noqa: E402
from common.result_io import load_cached_result  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

LOGGER = logging.getLogger("plot_boruta_importance_add_vs_other")

__script_meta__: dict[str, object] = {
    "task": "plot-boruta-importance-add-vs-other",
    "kind": "figure",
    "provides": [
        "build_boruta_figure",
        "build_legend_figure",
        "confirmed_directions",
        "main",
    ],
    "uses": [
        "common.figure_io",
        "common.colors",
        "common.result_io",
        "boruta_selection",
    ],
    "seeded_from": {"template": "boruta-importance", "version": "0.1"},
    "description": (
        "Boruta all-relevant importance figure for the cached ADD-vs-Other selection "
        "(results/boruta/762e16f294e3). Horizontal per-iteration box plots for the "
        "selected head (all 40 Confirmed + 15 Tentative + top-N Rejected) ordered by "
        "Boruta ranking, gene-symbol labels; the 40 Confirmed colored by DE direction "
        "(up/down in ADD via primary_contrast) rather than the template's viridis, "
        "Tentative/Rejected by decision (house-style); the median max-shadow noise "
        "floor as a dashed line and the shadow-max IQR as a band. Categorical coloring "
        "-> a SEPARATE legend image. Load -> plot; recomputes nothing. EXPLORATORY; "
        "ADD == Kerchner cohort confound UNRESOLVED."
    ),
}

# --- Cached-result identity (results/manifest.md; provenance for the finding). ---
DEFAULT_CACHE_ROOT = Path("results")
BORUTA_ANALYSIS = "boruta"
BORUTA_FINGERPRINT = "762e16f294e3"
DATA_VERSION = "protein-data.csv sha256:90a0462eb03a02bd"

DEFAULT_DE_TABLE = (
    Path("results")
    / "differential_abundance"
    / "add-vs-other-moderated"
    / "contrast_table.csv"
)
DEFAULT_REGISTRY = Path("state") / "color_registry.json"
DEFAULT_OUTPUT_DIR = Path("figures") / "analysis" / "boruta" / "add-vs-other"
DEFAULT_BASE_NAME = "boruta-importance-add-vs-other"
DEFAULT_TOP_REJECTED = 10

# House-style (non-categorical) colors. The Boruta DECISION states (Tentative/Rejected)
# and the shadow reference encode no metadata category, so they are house-style
# constants, not registry slots (mirrors lib/figures/boruta_importance's rationale).
_TENTATIVE_COLOR = "#F0E442"  # Okabe-Ito yellow — the canonical Boruta "tentative"
_REJECTED_COLOR = "#DDDDDD"  # light gray — deliberately lighter than Other's #999999
_SHADOW_LINE_COLOR = "#000000"  # black dashed — the max-shadow noise floor (median)
_SHADOW_BAND_COLOR = "#bdbdbd"  # per-iteration shadow-max IQR band
_MAX_LABEL_LEN = 22


def gene_label(name: str) -> str:
    """Collapse a UniProt ``sp|ACC|GENE_SPECIES`` id to its gene symbol (else raw)."""
    parts = str(name).split("|")
    if len(parts) >= 3 and parts[2]:
        symbol = parts[2].rsplit("_", 1)[0] if "_" in parts[2] else parts[2]
        text = symbol or str(name)
    else:
        text = str(name)
    if len(text) > _MAX_LABEL_LEN:
        return text[: _MAX_LABEL_LEN - 1] + "…"
    return text


def confirmed_directions(
    result: BorutaResult, de_table_path: str | Path = DEFAULT_DE_TABLE
) -> dict[int, str]:
    """Map each Confirmed feature index -> "up"/"down" in ADD from the DE effect sign.

    The moderated contrast's ``effect`` is the ``ADD vs Other`` log2 fold change, so a
    positive effect is up-in-ADD, negative is down-in-ADD. Fails loud if a Confirmed
    feature is absent from the DE table (the two analyses must share the feature space).
    """
    de = pd.read_csv(Path(de_table_path))
    effect = dict(
        zip(de["feature"].astype(str), de["effect"].astype(float), strict=True)
    )
    names = np.asarray(result.feature_names)
    decisions = np.asarray(result.decision)
    directions: dict[int, str] = {}
    missing: list[str] = []
    for idx in np.where(decisions == "Confirmed")[0]:
        feat = str(names[idx])
        if feat not in effect:
            missing.append(feat)
            continue
        directions[int(idx)] = "up" if effect[feat] > 0.0 else "down"
    if missing:
        raise ValueError(
            f"{len(missing)} Confirmed feature(s) absent from the DE table "
            f"{de_table_path}: {missing[:5]}{'...' if len(missing) > 5 else ''}."
        )
    return directions


def _selected_order(result: BorutaResult, top_rejected: int) -> np.ndarray:
    """Feature indices to draw, top (most important) to bottom.

    Boruta ranking ascending, then median importance descending (the template's order),
    which places the 40 Confirmed, then the 15 Tentative, then the Rejected in a clean
    block sequence. Keeps every Confirmed + Tentative and the first ``top_rejected``
    Rejected as a noise-floor anchor.
    """
    medians = np.asarray(result.importance, dtype=float)
    ranking = np.asarray(result.ranking, dtype=float)
    decisions = np.asarray(result.decision)
    order = np.lexsort((-medians, ranking))
    keep: list[int] = []
    rejected_kept = 0
    for idx in order:
        if decisions[idx] == "Rejected":
            if rejected_kept >= top_rejected:
                continue
            rejected_kept += 1
        keep.append(int(idx))
    return np.asarray(keep, dtype=int)


def _row_color(
    idx: int, decision: str, directions: dict[int, str], c_up: str, c_down: str
) -> str:
    if decision == "Confirmed":
        return c_up if directions.get(idx) == "up" else c_down
    if decision == "Tentative":
        return _TENTATIVE_COLOR
    return _REJECTED_COLOR


def build_boruta_figure(
    result: BorutaResult,
    directions: dict[int, str],
    *,
    top_rejected: int = DEFAULT_TOP_REJECTED,
    c_up: str,
    c_down: str,
) -> Figure:
    """Horizontal per-iteration box plots for the selected head, by DE direction."""
    if top_rejected < 0:
        raise ValueError(f"top_rejected must be >= 0; got {top_rejected}.")
    order = _selected_order(result, top_rejected)
    n_show = order.size
    history = np.asarray(result.importance_history, dtype=float)  # (n_iter, n_features)
    decisions = np.asarray(result.decision)
    names = np.asarray(result.feature_names)
    shadow = np.asarray(result.shadow_max_history, dtype=float)
    shadow_lo, shadow_hi = np.nanpercentile(shadow, [25.0, 75.0])

    n_up = sum(1 for v in directions.values() if v == "up")
    n_down = sum(1 for v in directions.values() if v == "down")

    with publication_style():
        fig, ax = plt.subplots(figsize=(11.0, max(8.0, 0.30 * n_show)))
        try:
            # Shadow-max reference: the per-iteration IQR band + the median noise floor.
            ax.axvspan(
                float(shadow_lo), float(shadow_hi), color=_SHADOW_BAND_COLOR,
                alpha=0.30, zorder=0,
            )
            ax.axvline(
                result.shadow_threshold, color=_SHADOW_LINE_COLOR, linestyle="--",
                linewidth=1.4, zorder=5,
            )

            for pos, feat in enumerate(order):
                decision = str(decisions[feat])
                color = _row_color(feat, decision, directions, c_up, c_down)
                alpha = 0.80 if decision == "Confirmed" else (
                    0.75 if decision == "Tentative" else 0.55
                )
                valid = history[:, feat][~np.isnan(history[:, feat])]
                if valid.size == 0:
                    continue
                box = ax.boxplot(
                    valid, positions=[pos], orientation="horizontal", widths=0.62,
                    patch_artist=True, showfliers=False,
                    medianprops={"color": "black", "linewidth": 1.1},
                    whiskerprops={"color": "#555555", "linewidth": 0.7},
                    capprops={"color": "#555555", "linewidth": 0.7},
                    boxprops={"edgecolor": "#333333", "linewidth": 0.7},
                    zorder=3,
                )
                for patch in box["boxes"]:
                    patch.set_facecolor(color)
                    patch.set_alpha(alpha)

            _draw_region_dividers(ax, order, decisions)

            ax.set_yticks(range(n_show))
            ax.set_yticklabels([gene_label(str(names[f])) for f in order], fontsize=8)
            ax.set_ylim(n_show - 0.5, -0.5)
            ax.set_xlim(left=0.0)
            ax.set_xlabel("Boruta feature importance (per iteration, 99 iterations)")
            ax.grid(axis="x", color="#e5e5e5", linewidth=0.6, zorder=0)
            # Threshold value labelled at the top, in the whitespace left of the boxes
            # (blended transform: x in data at the line, y in axes fraction).
            ax.text(
                result.shadow_threshold, 0.995,
                f"  max-shadow floor (median) = {result.shadow_threshold:.4f}",
                transform=ax.get_xaxis_transform(),
                fontsize=7.5, color="#333333", va="top", ha="left", zorder=6,
            )
            _apply_title(fig, result, n_up, n_down, top_rejected)
            fig.tight_layout(rect=(0, 0.02, 1, 0.96))
        except BaseException:
            plt.close(fig)
            raise
    return fig


def _draw_region_dividers(ax: object, order: np.ndarray, decisions: np.ndarray) -> None:
    """Faint separators + right-margin labels between the decision blocks."""
    labels = [str(decisions[f]) for f in order]
    boundaries: list[tuple[int, str]] = []
    for pos in range(1, len(labels)):
        if labels[pos] != labels[pos - 1]:
            boundaries.append((pos, labels[pos]))
    # Block spans for the region labels (start pos, end pos, name).
    spans: list[tuple[int, int, str]] = []
    start = 0
    for pos in range(1, len(labels) + 1):
        if pos == len(labels) or labels[pos] != labels[start]:
            spans.append((start, pos - 1, labels[start]))
            start = pos
    for pos, _ in boundaries:
        ax.axhline(pos - 0.5, color="#bbbbbb", linewidth=0.7, zorder=1)  # type: ignore[attr-defined]
    for lo, hi, name in spans:
        ax.text(  # type: ignore[attr-defined]
            1.005,
            1.0 - (lo + hi) / 2.0 / (len(labels) - 1) if len(labels) > 1 else 0.5,
            name, transform=ax.transAxes, rotation=90, va="center", ha="left",  # type: ignore[attr-defined]
            fontsize=8, color="#444444", fontweight="bold",
        )


def _apply_title(
    fig: Figure, result: BorutaResult, n_up: int, n_down: int, top_rejected: int
) -> None:
    counts = (
        f"Confirmed {result.n_confirmed} ({n_up} up / {n_down} down in ADD)  |  "
        f"Tentative {result.n_tentative}  |  Rejected {result.n_rejected} "
        f"(top {top_rejected} shown)"
    )
    fig.suptitle(
        f"Boruta all-relevant selection: ADD vs Other\n{counts}",
        fontsize=13, fontweight="bold",
    )
    # Terse caveat MARKERS (hit-counts / caveats, not a paragraph — the full narrative
    # is the caption). Two lines to keep each within the figure width.
    fig.text(
        0.5, 0.005,
        "EXPLORATORY (n=10 ADD)  ·  shadow = built-in permuted-copy null, no CV  ·  "
        "HARD CONFOUND: ADD = 100% Kerchner cohort (finding 0001)\n"
        "7 Confirmed are redundant-but-real (EN-zeroed)",
        ha="center", va="bottom", fontsize=8.5, color="#8a3d00",
    )


def build_legend_figure(
    result: BorutaResult,
    directions: dict[int, str],
    *,
    top_rejected: int,
    c_up: str,
    c_down: str,
) -> Figure:
    """Standalone categorical key (swatches + shadow references) — the legend image."""
    n_up = sum(1 for v in directions.values() if v == "up")
    n_down = sum(1 for v in directions.values() if v == "down")
    handles = [
        Patch(facecolor=c_up, edgecolor="#333333", alpha=0.80,
              label=f"Confirmed — up in ADD ({n_up})"),
        Patch(facecolor=c_down, edgecolor="#333333", alpha=0.80,
              label=f"Confirmed — down in ADD ({n_down})"),
        Patch(facecolor=_TENTATIVE_COLOR, edgecolor="#333333", alpha=0.75,
              label=f"Tentative ({result.n_tentative})"),
        Patch(facecolor=_REJECTED_COLOR, edgecolor="#333333", alpha=0.55,
              label=f"Rejected (top {top_rejected} of {result.n_rejected} shown)"),
        Line2D(
            [0], [0], color=_SHADOW_LINE_COLOR, linestyle="--", linewidth=1.4,
            label=f"max-shadow noise floor (median) = {result.shadow_threshold:.4f}",
        ),
        Patch(facecolor=_SHADOW_BAND_COLOR, alpha=0.30,
              label="shadow-max IQR (per-iteration spread)"),
    ]
    with publication_style():
        legend_fig = plt.figure(figsize=(5.4, 2.1))
        try:
            legend_fig.legend(
                handles=handles, loc="center", frameon=False, ncol=1, fontsize=9,
                title="Boruta decision · DE direction", title_fontsize=10,
            )
        except BaseException:
            plt.close(legend_fig)
            raise
    return legend_fig


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(name)s: %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-root", type=Path, default=DEFAULT_CACHE_ROOT)
    parser.add_argument("--de-table", type=Path, default=DEFAULT_DE_TABLE)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--base-name", type=str, default=DEFAULT_BASE_NAME)
    parser.add_argument("--top-rejected", type=int, default=DEFAULT_TOP_REJECTED)
    args = parser.parse_args(argv)

    result = load_cached_result(
        BorutaResult,
        cache_root=args.cache_root,
        analysis=BORUTA_ANALYSIS,
        fingerprint=BORUTA_FINGERPRINT,
    )
    LOGGER.info(
        "loaded cached Boruta %s: %d Confirmed / %d Tentative / %d Rejected; "
        "shadow_threshold=%.5f; %d iters; data_version=%s",
        BORUTA_FINGERPRINT, result.n_confirmed, result.n_tentative, result.n_rejected,
        result.shadow_threshold, result.n_iter, DATA_VERSION,
    )

    # Registry colors for the DE-direction categories (read, not mutated).
    c_up = get_color("primary_contrast", "ADD", registry_path=args.registry,
                     persist=False)
    c_down = get_color("primary_contrast", "Other", registry_path=args.registry,
                       persist=False)
    directions = confirmed_directions(result, args.de_table)
    LOGGER.info(
        "Confirmed DE direction: %d up / %d down in ADD (colors %s / %s)",
        sum(1 for v in directions.values() if v == "up"),
        sum(1 for v in directions.values() if v == "down"), c_up, c_down,
    )

    fig = build_boruta_figure(
        result, directions, top_rejected=args.top_rejected, c_up=c_up, c_down=c_down,
    )
    legend_fig = build_legend_figure(
        result, directions, top_rejected=args.top_rejected, c_up=c_up, c_down=c_down,
    )
    artifacts = save_figure(
        fig, args.output_dir, args.base_name, legend_fig=legend_fig, dpi=300,
    )
    for path in (
        artifacts.svg, artifacts.png, artifacts.legend_svg, artifacts.legend_png,
    ):
        LOGGER.info("wrote %s", path)


if __name__ == "__main__":
    main()
