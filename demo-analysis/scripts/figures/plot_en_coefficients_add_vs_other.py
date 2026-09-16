"""Coefficient / selection-frequency panel for the NULL-PASSED ADD-vs-Other elastic net.

Renders the CACHED, CONVERGED, NULL-PASSED elastic-net logistic ClassificationResult
(``results/classification/902c1e98c018/``, fingerprint ``902c1e98c018``, produced by
``scripts/scratch/classify_add_vs_other.py`` with ``run_null=True`` and
``max_iter=20000``) — a *load -> plot* figure job; it recomputes NOTHING
(conventions/results-cache.md). The nested-CV / stability loop is read from the cache
via ``common.result_io.load_cached_result``.

This SUPERSEDES the first-pass panel (fingerprint ``cbd809ac6168``, stem
``en-coefficients-add-vs-other-firstpass``), whose on-canvas caveats
("not tested against a null", "magnitudes approximate / saga not converged") are now
outdated. The selected set is unchanged (same 257 features, best_C=1.0, l1_ratio=0.25,
coefficients byte-identical to the first pass); only the licensing caveats change.

The panel: the top 25 selected proteins (by |all-data standardized coefficient|; these
are also the most stable — selection_frequency >= 0.88, sign_consistency = 1.00), each a
diamond at its SIGNED final coefficient over a bar spanning the resample IQR
(``coef_q25``..``coef_q75``). Direction is the read: positive = toward ADD, negative =
toward Other. Most of the top set (16/25) is DOWN in ADD (negative).

Encoding (colors read from ``state/color_registry.json``; registry read, not mutated):
  * direction via the ``primary_contrast`` category — positive coef (toward ADD) =
    ADD vermillion #D55E00; negative coef (toward Other) = Other gray #999999. Direction
    is ALSO encoded by position (left/right of the x=0 line), so color only reinforces.
  * the resample IQR bar (q25..q75) is drawn in the same direction color, lighter.
  * selection frequency is shown as a terse right-margin number per protein (all >= 0.88
    here — the top set is essentially always selected).
Two categorical colors only (well within the >8-category budget).

Proteins are labeled by gene symbol from a curated accession->symbol map (consistent
with the volcano's top-hit labels); an unmapped accession falls back to the UniProt
entry-name mnemonic present in the feature id (never fabricated).

LOAD-BEARING HONESTY (title markers + caption; state/DATA_DESCRIPTION.md, findings
0001/0004/0006):
  * NULL PASSED — the label-shuffle null HAS been run and passed (observed AUC 0.956 vs
    null mean ~0.50, permutation p = 0.001): the discrimination is real, not
    leakage/overfitting, so these coefficients ARE now licensed as the model's
    discriminative basis (no longer "not tested against a null").
  * CONVERGED — the saga solver converged (max_iter=20000), so coefficient MAGNITUDES
    are now trustworthy (the "magnitudes directional only" caveat is dropped); SIGNS
    remain fully stable (sign_consistency = 1.00).
  * STILL a MINIMAL-OPTIMAL set — a low selection frequency != unimportant; Boruta is
    the all-relevant complement (not yet run).
  * STILL CONFOUND-BOUNDED — ADD = 100% Kerchner (finding 0001); above-null does NOT
    establish disease (the null cannot separate ADD from collection site), so these may
    be site- as much as disease-discriminating.
  * 23/25 of the top set are also DE hits (finding 0006) with 100% coef-sign vs
    DE-log2FC direction agreement.

Run:  ./.venv/bin/python scripts/figures/plot_en_coefficients_add_vs_other.py
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

# Reuse the project shared figure-IO + color-registry modules, and the cached-result
# loader + the ClassificationResult dataclass (in scripts/scratch/).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scratch"))
from classification import ClassificationResult
from common.figure_io import publication_style, save_figure
from common.result_io import load_cached_result

LOGGER = logging.getLogger("plot_en_coefficients_add_vs_other")

DEFAULT_CACHE_ROOT = Path("results")
DEFAULT_ANALYSIS = "classification"
DEFAULT_FINGERPRINT = "902c1e98c018"
DEFAULT_OUTPUT_DIR = Path("figures") / "analysis" / "classification" / "add-vs-other"
DEFAULT_BASE_NAME = "en-coefficients-add-vs-other-nullpass"
DEFAULT_REGISTRY = Path("state") / "color_registry.json"
DEFAULT_TOP_N = 25

# primary_contrast direction -> which class the coefficient sign points toward.
DIR_ADD = "toward ADD (+)"
DIR_OTHER = "toward Other (−)"  # noqa: RUF001 (minus sign is intentional)

# Curated accession -> gene symbol for the top selected proteins. The 6 that also appear
# in the volcano's top-hit labels reuse its exact strings (MEK2/SHP1/USO1/VPS18/IVD/
# NDUFV2) for cross-figure consistency; the rest are HGNC symbols. An accession absent
# here falls back to the UniProt entry-name mnemonic (see _gene_label) — never invented.
GENE_LABELS: dict[str, str] = {
    "Q16822": "PCK2",
    "P29350": "SHP1",       # PTPN6 (volcano label)
    "Q9NRR3": "CDC42SE2",
    "Q6ZU65": "UBN2",
    "O60763": "USO1",
    "Q9HC98": "NEK6",       # shared-peptide group NEK6/NEK7 (see _gene_label)
    "Q8TDX7": "NEK7",
    "P36507": "MEK2",       # MAP2K2 (volcano label)
    "Q6NYC1": "JMJD6",
    "P62333": "PSMC6",
    "O14653": "GOSR2",
    "O43242": "PSMD3",
    "Q9P253": "VPS18",
    "P26440": "IVD",
    "Q9Y3C8": "UFC1",
    "Q8N543": "OGFOD1",
    "P19404": "NDUFV2",
    "P29692": "EEF1D",
    "Q7KZ85": "SUPT6H",
    "Q3MIN7": "RGL3",
    "P40939": "HADHA",
    "O75844": "ZMPSTE24",
    "P54819": "AK2",
    "Q13618": "CUL3",
    "P43686": "PSMC4",
    "Q02108": "GUCY1A1",
}


def _mnemonic(token: str) -> str:
    """UniProt entry-name mnemonic from one ``db|ACC|NAME_ORG`` id (else the token)."""
    parts = token.split("|")
    if len(parts) >= 3:
        name = parts[2]
        return name.rsplit("_", 1)[0] if "_" in name else name
    return token


def _gene_label(feature_id: str) -> str:
    """Gene-symbol label for a feature id, curated by accession; fail-soft to mnemonic.

    Handles shared-peptide groups (``id_a / id_b``) by joining the per-member labels
    with a slash (e.g. ``NEK6/NEK7``), so an ambiguous protein group is never silently
    collapsed to one member.
    """
    labels: list[str] = []
    for token in feature_id.split(" / "):
        parts = token.split("|")
        accession = parts[1] if len(parts) >= 3 else token
        labels.append(GENE_LABELS.get(accession, _mnemonic(token)))
    return "/".join(labels)


def _read_registry_color(category: str, value: str, registry_path: Path) -> str:
    """Hex color for ``(category, value)`` from the registry; read-only."""
    from common.colors import get_color

    return get_color(category, value, registry_path=registry_path, persist=False)


def build_coefficient_panel(
    result: ClassificationResult, registry_path: Path, *, top_n: int
) -> tuple[Figure, Figure, dict[str, str]]:
    """Build the coefficient panel + its separate legend figure. Registry read-only."""
    if top_n <= 0:
        raise ValueError(f"top_n must be positive; got {top_n}.")
    table = result.coefficients
    if len(table) == 0:
        raise ValueError("no features were selected (all coefficients zero).")

    # Largest |coef| at the top: take the head, then reverse so y increases upward.
    sub = table.head(top_n).iloc[::-1].reset_index(drop=True)
    n = len(sub)
    coef = sub["coef"].to_numpy(dtype=float)
    q25 = sub["coef_q25"].to_numpy(dtype=float)
    q75 = sub["coef_q75"].to_numpy(dtype=float)
    freq = sub["selection_frequency"].to_numpy(dtype=float)
    genes = [_gene_label(str(f)) for f in sub["feature"]]

    c_add = _read_registry_color("primary_contrast", "ADD", registry_path)
    c_other = _read_registry_color("primary_contrast", "Other", registry_path)
    color_map = {DIR_ADD: c_add, DIR_OTHER: c_other}
    row_color = [c_add if v > 0 else c_other for v in coef]

    with publication_style():
        fig, ax = plt.subplots(figsize=(7.8, max(3.0, 0.34 * n + 2.0)))

        for i in range(n):
            col = row_color[i]
            ax.plot([q25[i], q75[i]], [i, i], color=col, lw=2.4, alpha=0.45, zorder=2)
            ax.scatter(
                [coef[i]], [i], marker="D", s=52, color=col,
                edgecolor="black", linewidths=0.5, zorder=5,
            )
        ax.axvline(0.0, color="black", lw=0.9, zorder=1)

        ax.set_yticks(range(n))
        ax.set_yticklabels(genes, fontsize=8)
        ax.set_ylim(-0.7, n - 0.3)

        # x-range with right-margin room for the per-row selection-frequency numbers.
        xmin = float(min(q25.min(), coef.min()))
        xmax = float(max(q75.max(), coef.max()))
        span = xmax - xmin
        ax.set_xlim(xmin - 0.06 * span, xmax + 0.30 * span)
        freq_x = xmax + 0.10 * span
        for i in range(n):
            ax.text(freq_x, i, f"{freq[i]:.2f}", va="center", ha="left",
                    fontsize=7, color="#555555")

        ax.set_xlabel(
            "standardized elastic-net coefficient\n"
            "(←  toward Other        |        toward ADD  →)",
            fontsize=10,
        )

        title = (
            f"Elastic-net coefficients — ADD vs Other (top {n} of "
            f"{len(table)} selected)  ·  NULL-PASSED\n"
            "label-shuffle null passed (AUC 0.956 vs ~0.50, p=0.001) → licensed as the "
            "model's discriminative basis; saga converged, signs & magnitudes trusted\n"
            "HARD CONFOUND: ADD = 100% Kerchner (finding 0001) — above-null does not "
            "establish disease; may be site- not disease-discriminating"
        )
        ax.set_title(title, fontsize=10, fontweight="bold", loc="left")

        # One terse stability line (the numbers, not the story) at the foot of the axes.
        ax.annotate(
            f"top {n}: selection frequency ≥ {freq.min():.2f} · "
            "sign consistency = 1.00 · diamond = final coef, bar = resample IQR · "
            "right-margin number = selection frequency",
            xy=(0.0, -0.11), xycoords="axes fraction", ha="left", va="top",
            fontsize=8, color="#333333",
        )
        # Single terse minimal-optimal caveat marker (Boruta all-relevant not yet run) —
        # a low selection frequency is NOT unimportance; these 25/257 are not the
        # complete signature (LOAD-BEARING; findings 0006).
        ax.annotate(
            "minimal-optimal set; low sel-freq ≠ unimportance; "
            "Boruta (all-relevant) not run",
            xy=(0.0, -0.155), xycoords="axes fraction", ha="left", va="top",
            fontsize=8, color="#333333", style="italic",
        )

    legend_fig = _build_legend(c_add, c_other)
    return fig, legend_fig, color_map


def _build_legend(c_add: str, c_other: str) -> Figure:
    """Standalone swatch/line key for the coefficient panel, as its own figure image."""
    with publication_style():
        lf, la = plt.subplots(figsize=(4.8, 2.4))
        la.axis("off")
        handles = [
            mlines.Line2D([], [], marker="D", ls="", ms=9, mfc=c_add, mec="black",
                          mew=0.5, label="positive coef — toward ADD"),
            mlines.Line2D([], [], marker="D", ls="", ms=9, mfc=c_other, mec="black",
                          mew=0.5, label="negative coef — toward Other"),
            mlines.Line2D([], [], color="#888888", lw=2.4, alpha=0.6,
                          label="bar = resample IQR (q25..q75)"),
            mlines.Line2D([], [], marker="D", ls="", ms=9, mfc="white", mec="black",
                          mew=0.5, label="diamond = final all-data coef"),
            mlines.Line2D([], [], marker="$0.98$", ls="", ms=14, mfc="#555555",
                          mec="none", label="right number = selection frequency"),
        ]
        la.legend(handles=handles, loc="center", frameon=False, fontsize=10,
                  handletextpad=0.7, labelspacing=0.9, title="coefficient panel key",
                  title_fontsize=11)
    return lf


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(name)s: %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-root", type=Path, default=DEFAULT_CACHE_ROOT)
    parser.add_argument("--analysis", type=str, default=DEFAULT_ANALYSIS)
    parser.add_argument("--fingerprint", type=str, default=DEFAULT_FINGERPRINT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--base-name", type=str, default=DEFAULT_BASE_NAME)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--top-n", type=int, default=DEFAULT_TOP_N)
    args = parser.parse_args()

    result = load_cached_result(
        ClassificationResult,
        cache_root=args.cache_root,
        analysis=args.analysis,
        fingerprint=args.fingerprint,
    )
    LOGGER.info(
        "loaded cached %s/%s: %d selected features, %s vs %s, nested AUC=%.3f",
        args.analysis, args.fingerprint, len(result.coefficients),
        result.positive_label, result.negative_label, result.cv_auc,
    )

    fig, legend_fig, color_map = build_coefficient_panel(
        result, args.registry, top_n=args.top_n
    )
    artifacts = save_figure(
        fig, args.output_dir, args.base_name, legend_fig=legend_fig, dpi=300
    )
    LOGGER.info("color map: %s", color_map)
    for p in (artifacts.svg, artifacts.png, artifacts.legend_svg, artifacts.legend_png):
        LOGGER.info("wrote %s", p)


if __name__ == "__main__":
    main()
