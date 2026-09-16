"""Within-arm cohort degeneracy - a data-integrity diagnostic for the four-arm probe.

Makes ONE point in two panels: in the delivered (batch-corrected) matrix, *within* a
diagnostic arm the two collection cohorts have IDENTICAL per-protein means (to machine
precision), even though individual samples genuinely differ. So a within-arm site test
is structurally UNINFORMATIVE here - cohort is orthogonal to abundance given Condition,
and a null within-arm site probe is not evidence of "no site effect".

  * LEFT  - within HCN, a jittered strip of the 10 per-sample complement/coagulation
    summary scores (KEGG-27), split Kerchner (n=7) vs Poston (n=3), each cohort's MEAN
    drawn as a horizontal segment. The two mean segments visibly COINCIDE despite the
    scattered points. Annotated with the max abs per-protein mean diff (2340 proteins).
  * RIGHT - a histogram of the per-protein within-HCN Kerchner-minus-Poston mean
    difference across ALL 2340 proteins: all mass sits at ~0 (values ~1e-14), plotted on
    a symmetric +/-1e-13 axis so the degeneracy reads as universal, not just the score.

DATA / PROVENANCE
-----------------
The 27-protein per-sample summary scores + each sample's cohort are read STRAIGHT from
the cached four-arm result (``results/four_arm_complement/dd884570c4fc/``, fingerprint
``dd884570c4fc``, data version ``protein-data.csv sha256:90a0462eb03a02bd``,
``sample_scores.parquet``) - no recomputation of that result.

The per-protein histogram needs the SIGNED per-protein cohort mean difference across all
2340 proteins, which the cache does NOT store (``cohort_orthogonality.parquet`` holds
only the max/mean abs-diff SUMMARY). That signed vector is therefore recomputed here
directly from the delivered matrix via the Stage-3 loader
(``load_analysis_subset(load_paired_dataset())``, scripts/scratch/project_loader.py) -
the exact ``X[HCN&Kerchner].mean(0) - X[HCN&Poston].mean(0)`` that
``four_arm_complement.py::_cohort_orthogonality`` reduced to its cached summary. This is
a trivial deterministic group-mean subtraction (NOT the CPU result), and the script
FAILS LOUD if the recomputed max abs-diff does not reproduce the cached
``max_abs_protein_mean_diff`` - so a data/loader drift can never silently mis-draw it.

Encoding - colors read from ``state/color_registry.json`` under the ``cohort`` category
(registry read, never mutated): point fill = collection cohort (Kerchner #56B4E9, Poston
#CC79A7). The cohort mean is a black horizontal segment; a faint full-width grey line
marks the common mean so the coincidence is unmistakable. Two cohort colors - well
within the 8-color budget.

Run::

    ./.venv/bin/python scripts/figures/plot_within_arm_cohort_degeneracy.py
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure

# Reuse the project shared figure-IO + color-registry modules and the Stage-3 loader.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scratch"))
from common.colors import get_color
from common.figure_io import publication_style, save_figure

LOGGER = logging.getLogger("plot_within_arm_cohort_degeneracy")

DATA_VERSION = "protein-data.csv sha256:90a0462eb03a02bd"
FINGERPRINT = "dd884570c4fc"

DEFAULT_CACHE_DIR = Path("results") / "four_arm_complement" / FINGERPRINT
DEFAULT_OUTPUT_DIR = (
    Path("figures") / "analysis" / "four_arm_complement" / "complement-coagulation"
)
DEFAULT_BASE_NAME = "within-arm-cohort-degeneracy"
DEFAULT_REGISTRY = Path("state") / "color_registry.json"
DEFAULT_DATA = Path("data") / "protein-data.csv"
DEFAULT_METADATA = Path("data") / "metadata.csv"

# The HCN within-arm site probe: Kerchner (majority) vs Poston.
ARM = "HCN"
SITE_A = "Kerchner"
SITE_B = "Poston"
SCORE_COL = "complement_score"

# Deterministic strip jitter (reproducible render).
_RNG_SEED = 20260915
# Tolerance when asserting the recomputed per-protein max |Δ| reproduces the cache.
_MAX_DELTA_ATOL = 1e-16


def _per_protein_delta(data: Path, metadata: Path) -> np.ndarray:
    """Signed per-protein within-HCN (Kerchner minus Poston) mean diff, 2340 proteins.

    Recomputes the exact vector that ``four_arm_complement.py::_cohort_orthogonality``
    reduced to its cached max/mean summary — from the delivered (batch-corrected) matrix
    via the Stage-3 loader. Deterministic; no RNG.
    """
    from project_loader import load_analysis_subset, load_paired_dataset

    dataset = load_analysis_subset(load_paired_dataset(data, metadata))
    if dataset.scale != "log2":
        raise ValueError(f"expected log2-scale matrix, got scale={dataset.scale!r}")
    meta = dataset.metadata
    arm = meta["Condition"].astype(str).to_numpy()
    coh = meta["Cohort"].astype(str).to_numpy()
    x = dataset.abundances
    m_a = (arm == ARM) & (coh == SITE_A)
    m_b = (arm == ARM) & (coh == SITE_B)
    if m_a.sum() < 2 or m_b.sum() < 2:
        raise ValueError(
            f"within-{ARM} cells too small: {SITE_A} n={int(m_a.sum())}, "
            f"{SITE_B} n={int(m_b.sum())}"
        )
    LOGGER.info(
        "within-%s: %s n=%d, %s n=%d, %d proteins",
        ARM,
        SITE_A,
        int(m_a.sum()),
        SITE_B,
        int(m_b.sum()),
        x.shape[1],
    )
    return np.asarray(x[m_a].mean(axis=0) - x[m_b].mean(axis=0), dtype=float)


def _load_scores(cache_dir: Path) -> pd.DataFrame:
    """Per-sample HCN summary scores + cohort, straight from the cached result."""
    scores = pd.read_parquet(cache_dir / "sample_scores.parquet")
    hcn = scores.loc[scores["arm"] == ARM].copy()
    n_a = int((hcn["cohort"] == SITE_A).sum())
    n_b = int((hcn["cohort"] == SITE_B).sum())
    if n_a == 0 or n_b == 0:
        raise ValueError(f"missing HCN cohort cells: {SITE_A}={n_a}, {SITE_B}={n_b}")
    return hcn


def make_figure(
    hcn: pd.DataFrame,
    delta: np.ndarray,
    *,
    registry_path: Path,
) -> tuple[Figure, Figure]:
    """Build the two-panel figure + a standalone legend figure."""
    color_a = get_color("cohort", SITE_A, registry_path=registry_path, persist=False)
    color_b = get_color("cohort", SITE_B, registry_path=registry_path, persist=False)
    cohort_color = {SITE_A: color_a, SITE_B: color_b}

    max_abs = float(np.abs(delta).max())
    mean_abs = float(np.abs(delta).mean())

    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(10.6, 4.6))

    # ---- LEFT: within-HCN per-sample strip, means coincide -------------------------
    rng = np.random.default_rng(_RNG_SEED)
    x_pos = {SITE_A: 0.0, SITE_B: 1.0}
    means: dict[str, float] = {}
    for site in (SITE_A, SITE_B):
        vals = hcn.loc[hcn["cohort"] == site, SCORE_COL].to_numpy(dtype=float)
        means[site] = float(vals.mean())
        jitter = rng.uniform(-0.11, 0.11, size=vals.size)
        ax_l.scatter(
            np.full(vals.size, x_pos[site]) + jitter,
            vals,
            s=64,
            facecolor=cohort_color[site],
            edgecolor="black",
            linewidth=0.6,
            alpha=0.95,
            zorder=3,
        )

    # The two cohort means coincide to machine precision: one faint full-width line at
    # the common mean, plus a solid black mean segment over each cohort's x-band.
    common_mean = float(np.mean(list(means.values())))
    ax_l.axhline(
        common_mean, color="0.55", linestyle=(0, (4, 3)), linewidth=1.0, zorder=1
    )
    for site in (SITE_A, SITE_B):
        ax_l.plot(
            [x_pos[site] - 0.28, x_pos[site] + 0.28],
            [means[site], means[site]],
            color="black",
            linewidth=2.4,
            solid_capstyle="butt",
            zorder=4,
        )
    mean_gap = abs(means[SITE_A] - means[SITE_B])

    ax_l.set_xlim(-0.6, 1.6)
    ax_l.set_xticks([x_pos[SITE_A], x_pos[SITE_B]])
    ax_l.set_xticklabels([f"{SITE_A}\n(n=7)", f"{SITE_B}\n(n=3)"])
    ax_l.set_ylabel("Complement/coag. summary score\n(mean z, KEGG-27)")
    ax_l.set_title("Within HCN: points differ, means coincide", fontsize=12)
    ax_l.annotate(
        f"cohort mean gap = {mean_gap:.0e}\nper-protein max |Δ| = {max_abs:.1e}"
        f"\n(all 2340 proteins)",
        xy=(0.5, 0.02),
        xycoords="axes fraction",
        ha="center",
        va="bottom",
        fontsize=9,
        color="0.25",
    )

    # ---- RIGHT: per-protein Δ histogram, all mass at ~0 ----------------------------
    xlim = 1e-13
    bins = np.linspace(-xlim, xlim, 61)
    ax_r.hist(delta, bins=bins, color="#0072B2", edgecolor="white", linewidth=0.3)
    ax_r.axvline(0.0, color="black", linewidth=1.0, linestyle=(0, (4, 3)), zorder=3)
    ax_r.set_xlim(-xlim, xlim)
    # Plot labels use the typographic minus (U+2212) and Greek Delta / approx glyphs on
    # purpose - they render correctly in the published figure; silence ruff's ambiguous-
    # Unicode check for these display strings only.
    ax_r.set_xlabel("Per-protein cohort mean Δ\n(Kerchner − Poston, log2)")  # noqa: RUF001
    ax_r.set_ylabel("Proteins (count)")
    ax_r.set_title("All 2340 proteins: cohort mean Δ ≈ 0", fontsize=12)
    n_prot = int(delta.size)
    ax_r.annotate(
        f"n = {n_prot} proteins\nmax |Δ| = {max_abs:.1e}\nmean |Δ| = {mean_abs:.1e}",
        xy=(0.03, 0.97),
        xycoords="axes fraction",
        ha="left",
        va="top",
        fontsize=9,
        color="0.25",
    )

    fig.tight_layout(w_pad=3.0)

    # ---- separate legend image ------------------------------------------------------
    legend_fig = plt.figure(figsize=(3.4, 1.2))
    handles = [
        mpatches.Patch(facecolor=color_a, edgecolor="black", linewidth=0.6,
                       label=f"{SITE_A} (n=7)"),
        mpatches.Patch(facecolor=color_b, edgecolor="black", linewidth=0.6,
                       label=f"{SITE_B} (n=3)"),
        mlines.Line2D([], [], color="black", linewidth=2.4, label="cohort mean"),
    ]
    legend_fig.legend(
        handles=handles,
        loc="center",
        frameon=False,
        title="Collection cohort (within HCN)",
        ncol=1,
    )
    legend_fig.tight_layout()

    return fig, legend_fig


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--base-name", type=str, default=DEFAULT_BASE_NAME)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    args = parser.parse_args(argv)

    hcn = _load_scores(args.cache_dir)
    delta = _per_protein_delta(args.data, args.metadata)

    # Fidelity gate: the recomputed per-protein Δ MUST reproduce the cached summary, or
    # the loader/data has drifted and the panel would be wrong. Fail loud.
    orth = pd.read_parquet(args.cache_dir / "cohort_orthogonality.parquet")
    cached_max = float(
        orth.loc[orth["arm"] == ARM, "max_abs_protein_mean_diff"].iloc[0]
    )
    recomputed_max = float(np.abs(delta).max())
    if not np.isclose(recomputed_max, cached_max, rtol=1e-6, atol=_MAX_DELTA_ATOL):
        raise ValueError(
            f"recomputed within-{ARM} max |Δ| ({recomputed_max:.6e}) does not match "
            f"cached max_abs_protein_mean_diff ({cached_max:.6e}); data/loader drift?"
        )
    LOGGER.info(
        "fidelity OK: recomputed max |Δ|=%.3e matches cached %.3e",
        recomputed_max,
        cached_max,
    )

    with publication_style():
        fig, legend_fig = make_figure(hcn, delta, registry_path=args.registry)
        artifacts = save_figure(
            fig, args.output_dir, args.base_name, legend_fig=legend_fig
        )

    LOGGER.info(
        "wrote:\n  %s\n  %s\n  %s\n  %s",
        artifacts.svg,
        artifacts.png,
        artifacts.legend_svg,
        artifacts.legend_png,
    )
    provenance = {
        "script": "scripts/figures/plot_within_arm_cohort_degeneracy.py",
        "data_version": DATA_VERSION,
        "fingerprint": FINGERPRINT,
        "cache_dir": str(args.cache_dir),
        "recomputed_max_abs_delta": recomputed_max,
        "cached_max_abs_delta": cached_max,
    }
    LOGGER.info("PROVENANCE: %s", json.dumps(provenance))


if __name__ == "__main__":
    main()
