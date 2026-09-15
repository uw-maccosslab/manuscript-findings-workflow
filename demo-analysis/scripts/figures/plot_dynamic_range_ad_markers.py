"""Rank-abundance (dynamic-range) figure with known AD / EV markers highlighted.

A SECOND dynamic-range render (the plain QC curve produced by ``dynamic_range.py`` is
left untouched). Same base curve, same input state
(``results/qc_states/delivered_linear`` -- the de-logged image of the delivered
normalized + batch-corrected matrix, 40 x 2340, scale=linear; NOT raw intensity), but it
marks 11 detected proteins, color- AND shape-coded by three categories, so the figure
doubles as an EV-enrichment QC check and answers: are known AD / EV markers even
detectable in this untargeted plasma-EV DIA matrix?

The three highlight categories (registry colors + a distinct marker as a second channel):
  * EV-QC anchor (MISEV2018) -- canonical EV-identity/recovery markers. Detecting them is
    the POSITIVE result: 5 of the 7 canonical markers (CD9, syntenin/SDCBP, ALIX/PDCD6IP,
    TSG101, flotillin-1/FLOT1) are present -> consistent with successful EV enrichment;
    the tetraspanins CD63 and CD81 are NOT detected.
  * AD-marker parent protein -- GFAP and APP. HONESTY: only the PARENT protein is detected;
    this is NOT the diagnostic form. APP (sp|P05067|A4_HUMAN) detected != the Abeta42/40
    proteolytic fragment (which needs immunoprecipitation-MS, not untargeted DIA); GFAP
    protein-level ID != a validated glial biomarker readout.
  * EV-cargo / candidate (some contested) -- L1CAM and NCAM1 are the contested
    neuron-derived-EV capture markers (Norman 2021: most plasma L1CAM is soluble, not
    EV-bound); F13A1/F13B (coagulation factor XIII) come from a single small non-replicated
    plasma-EV DDA study. None is an established plasma-EV AD-differential marker.

The load-bearing QC message (stated in the caption, verified against the matrix): the
classic ultralow-abundance ATN / glial / axonal / synaptic plasma AD markers are BELOW the
untargeted plasma DIA-MS detection floor here -- Tau (MAPT), NfL (NEFL), NRGN, CHI3L1,
GAP43, SNAP25, SYP, NPTX2, SYT1, S100A8, and EV-tetraspanins CD63 & CD81 are NOT detected.

Marker list + categorization is from the verified research corpus (research/): findings
`ad-plasma-atn-core-biomarkers`, `ad-plasma-glial-axonal-biomarkers`,
`ad-plasma-ev-proteomic-markers`.

Run:  ./.venv/bin/python scripts/figures/plot_dynamic_range_ad_markers.py
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

# Reuse the project dynamic-range module (compute + plot + save + color/marker legend).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.dataset_io import load_dataset  # noqa: E402
from common.figure_io import save_figure  # noqa: E402
from figures.dynamic_range import plot_dynamic_range  # noqa: E402

LOGGER = logging.getLogger("plot_dynamic_range_ad_markers")

DEFAULT_STATE_DIR = Path("results") / "qc_states" / "delivered_linear"
DEFAULT_OUTPUT_DIR = Path("figures") / "qc" / "dynamic-range"
DEFAULT_BASE_NAME = "dynamic-range-delivered-linear-ad-markers"
DEFAULT_REGISTRY = Path("state") / "color_registry.json"

# Three highlight categories (order fixes registry color assignment).
G_EV = "EV-QC anchor (MISEV2018)"
G_PARENT = "AD-marker parent protein"
G_CARGO = "EV-cargo / candidate (some contested)"

# feature id -> on-canvas label (terse; caveat glyphs "(parent)" / "*" keyed in a footnote)
HIGHLIGHT_FEATURES: dict[str, str] = {
    # EV-QC anchors (positive controls for EV enrichment)
    "sp|P21926|CD9_HUMAN": "CD9",
    "sp|O00560|SDCB1_HUMAN": "syntenin",
    "sp|Q8WUM4|PDC6I_HUMAN": "ALIX",
    "sp|Q99816|TS101_HUMAN": "TSG101",
    "sp|O75955|FLOT1_HUMAN": "FLOT1",
    # AD disease-marker parent proteins (parent detected only)
    "sp|P14136|GFAP_HUMAN": "GFAP (parent)",
    "sp|P05067|A4_HUMAN": "APP (parent)",
    # EV-cargo / candidates (L1CAM/NCAM1 contested)
    "sp|P32004|L1CAM_HUMAN": "L1CAM*",
    "sp|P13591|NCAM1_HUMAN": "NCAM1*",
    "sp|P00488|F13A_HUMAN": "F13A1",
    "sp|P05160|F13B_HUMAN": "F13B",
}

HIGHLIGHT_GROUPS: dict[str, str] = {
    "sp|P21926|CD9_HUMAN": G_EV,
    "sp|O00560|SDCB1_HUMAN": G_EV,
    "sp|Q8WUM4|PDC6I_HUMAN": G_EV,
    "sp|Q99816|TS101_HUMAN": G_EV,
    "sp|O75955|FLOT1_HUMAN": G_EV,
    "sp|P14136|GFAP_HUMAN": G_PARENT,
    "sp|P05067|A4_HUMAN": G_PARENT,
    "sp|P32004|L1CAM_HUMAN": G_CARGO,
    "sp|P13591|NCAM1_HUMAN": G_CARGO,
    "sp|P00488|F13A_HUMAN": G_CARGO,
    "sp|P05160|F13B_HUMAN": G_CARGO,
}

# Second channel: a distinct marker shape per category (colorblind-safe redundancy).
HIGHLIGHT_MARKERS: dict[str, str] = {G_EV: "o", G_PARENT: "s", G_CARGO: "^"}

# The two on-canvas caveat glyphs, defined as a terse footnote key (load-bearing honesty).
_CAVEAT = (
    "(parent) = parent protein detected, NOT the diagnostic Abeta42/40 fragment or "
    "p-tau phosphopeptide (needs IP-MS)      * = contested EV-capture marker (L1CAM / NCAM1)"
)


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
        "loaded %s: %d x %d, scale=%s",
        args.state_dir, n_samples, int(dataset.abundances.shape[1]), dataset.scale,
    )

    title = (
        "Known AD / EV markers on the rank-abundance curve\n"
        "delivered (normalized, batch-corrected) matrix, de-logged to linear "
        f"— {n_samples} experimental samples (processed data, not raw intensity)"
    )

    plot = plot_dynamic_range(
        dataset,
        highlight_features=HIGHLIGHT_FEATURES,
        highlight_groups=HIGHLIGHT_GROUPS,
        highlight_markers=HIGHLIGHT_MARKERS,
        min_intensity=0.0,
        title=title,
        legend_title="marker class",
        registry_path=args.registry,
        persist_colors=True,
    )

    # Terse caveat key for the two on-canvas glyphs ("(parent)" and "*"), placed under
    # the axes so it reads as a key, not interpretation.
    ax = plot.figure.axes[0]
    ax.annotate(
        _CAVEAT,
        xy=(0.0, -0.14),
        xycoords="axes fraction",
        ha="left",
        va="top",
        fontsize=8.5,
        color="#555555",
    )

    artifacts = save_figure(
        plot.figure,
        args.output_dir,
        args.base_name,
        legend_fig=plot.legend_figure,
        dpi=300,
    )

    rank_of = {
        name: i + 1 for i, name in enumerate(plot.result.feature_names_ranked)
    }
    for fid in HIGHLIGHT_FEATURES:
        LOGGER.info("  highlight %-26s rank %d / %d",
                    fid, rank_of[fid], plot.result.n_features_detected)
    LOGGER.info("color map: %s", plot.color_map)
    LOGGER.info("wrote %s", artifacts.svg)
    LOGGER.info("wrote %s", artifacts.png)
    LOGGER.info("wrote %s", artifacts.legend_svg)
    LOGGER.info("wrote %s", artifacts.legend_png)


if __name__ == "__main__":
    main()
