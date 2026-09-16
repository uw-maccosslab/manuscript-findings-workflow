"""GO:CC compartment composition: ADD-vs-Other discriminative set vs background.

The "show" for a NEGATIVE enrichment result. For a panel of informative GO:CC cellular
compartments, plot the FRACTION of proteins annotated to each compartment in three
groups side by side:

  * up-in-ADD query   (intersection / query_size; g:Profiler-mapped)
  * down-in-ADD query (intersection / query_size; g:Profiler-mapped)
  * detected plasma-EV proteome background (term_size / effective-domain)

The message: the query fractions TRACK the background fractions — no compartment is
selectively over-represented (consistent with all g:SCS-corrected p = 1.0). The apparent
"intracellular / organellar + exosomal signature" of the discriminative set reflects the
composition of the detected proteome, not selective enrichment.

Recomputes NO enrichment. Reads the already-computed g:Profiler term tables (custom
detected-proteome background):
  results/enrichment/add-vs-other-boruta-confirmed/up_in_ADD_terms.csv
  results/enrichment/add-vs-other-boruta-confirmed/down_in_ADD_terms.csv
  state/color_registry.json  (primary_contrast: ADD / Other — read, NOT mutated)

Denominators (stated on-figure + in provenance):
  * query fractions use g:Profiler's `query_size` column (mapped genes: 17 of 18 up,
    20 of 22 down) — the denominator `intersection_size` is counted against.
  * background fraction uses the g:Profiler effective domain size (2593 genes with any
    annotation under the custom background), NOT the 2340 raw background proteins.

Run:  ./.venv/bin/python scripts/figures/plot_gocc_composition_add_vs_other.py
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import logging
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.colors import BACKGROUND_COLOR, get_color
from common.figure_io import publication_style, save_figure

LOGGER = logging.getLogger("plot_gocc_composition_add_vs_other")

DEFAULT_ENRICH_DIR = (
    Path("results") / "enrichment" / "add-vs-other-boruta-confirmed"
)
DEFAULT_OUTPUT_DIR = (
    Path("figures") / "analysis" / "enrichment" / "add-vs-other-boruta-confirmed"
)
DEFAULT_BASE_NAME = "gocc-composition-add-vs-other-boruta-confirmed"
DEFAULT_REGISTRY = Path("state") / "color_registry.json"

# g:Profiler effective domain size for the custom detected-proteome background
# (genes with any annotation). Background fraction denominator.
DEFAULT_EFFECTIVE_DOMAIN = 2593

# Submitted query sizes (for the caption / legend); g:Profiler maps 17 of 18, 20 of 22.
SUBMITTED_N_UP = 18
SUBMITTED_N_DOWN = 22

# The informative GO:CC compartments to show, (term_id, display name). Chosen to span
# extracellular / vesicular / cytosolic / organellar / nuclear compartments without
# redundant nested parents; ordered later by background fraction (descending).
COMPARTMENTS: list[tuple[str, str]] = [
    ("GO:0005829", "cytosol"),
    ("GO:0031982", "vesicle"),
    ("GO:0005576", "extracellular region"),
    ("GO:0005634", "nucleus"),
    ("GO:0005886", "plasma membrane"),
    ("GO:0070062", "extracellular exosome"),
    ("GO:0005783", "endoplasmic reticulum"),
    ("GO:0005794", "Golgi apparatus"),
    ("GO:0005739", "mitochondrion"),
    ("GO:0000502", "proteasome complex"),
]

TITLE = "GO:CC compartment composition — ADD-vs-Other discriminative set vs background"
# Load-bearing caveat marker: the negative-result headline, kept to one line.
CAVEAT = (
    "all g:SCS-corrected p = 1.0 — no GO:CC term over-represented "
    "(query tracks background)"
)


def _load_gocc(path: Path) -> tuple[dict[str, dict[str, str]], int]:
    """Return {term_id: row} for GO:CC rows and the query_size (mapped genes)."""
    with path.open() as fh:
        rows = [r for r in csv.DictReader(fh) if r["source"] == "GO:CC"]
    if not rows:
        raise ValueError(f"{path} has no GO:CC rows.")
    query_sizes = {int(r["query_size"]) for r in rows}
    if len(query_sizes) != 1:
        raise ValueError(f"{path} has inconsistent query_size values: {query_sizes}.")
    return {r["term_id"]: r for r in rows}, query_sizes.pop()


def data_version(*paths: Path) -> str:
    """Short combined content hash of the input tables, for provenance."""
    h = hashlib.sha256()
    for p in sorted(paths):
        h.update(p.read_bytes())
    return h.hexdigest()[:16]


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(name)s: %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--enrich-dir", type=Path, default=DEFAULT_ENRICH_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--base-name", type=str, default=DEFAULT_BASE_NAME)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument(
        "--effective-domain", type=int, default=DEFAULT_EFFECTIVE_DOMAIN
    )
    args = parser.parse_args()

    up_path = args.enrich_dir / "up_in_ADD_terms.csv"
    down_path = args.enrich_dir / "down_in_ADD_terms.csv"
    dv = data_version(up_path, down_path)
    LOGGER.info("data_version(up+down terms)=%s", dv)

    up, q_up = _load_gocc(up_path)
    down, q_down = _load_gocc(down_path)
    dom = args.effective_domain
    LOGGER.info(
        "query_size(mapped): up=%d (of %d submitted), down=%d (of %d submitted); "
        "effective domain=%d",
        q_up,
        SUBMITTED_N_UP,
        q_down,
        SUBMITTED_N_DOWN,
        dom,
    )

    # Assemble fractions per compartment. Absent term => intersection 0. term_size is
    # the number of BACKGROUND proteins in the term (same in both tables where present).
    records = []
    for tid, name in COMPARTMENTS:
        u = up.get(tid)
        d = down.get(tid)
        src = u or d
        if src is None:
            raise ValueError(f"{tid} ({name}) not found in either GO:CC table.")
        term_size = int(src["term_size"])
        up_int = int(u["intersection_size"]) if u else 0
        dn_int = int(d["intersection_size"]) if d else 0
        records.append(
            {
                "term_id": tid,
                "name": name,
                "term_size": term_size,
                "up_int": up_int,
                "dn_int": dn_int,
                "up_frac": up_int / q_up,
                "dn_frac": dn_int / q_down,
                "bg_frac": term_size / dom,
            }
        )
    # Order by background fraction descending so the reference baseline reads cleanly.
    records.sort(key=lambda r: r["bg_frac"], reverse=True)

    for r in records:
        LOGGER.info(
            "%-22s bg=%.3f (%d/%d)  up=%.3f (%d/%d)  down=%.3f (%d/%d)",
            r["name"],
            r["bg_frac"],
            r["term_size"],
            dom,
            r["up_frac"],
            r["up_int"],
            q_up,
            r["dn_frac"],
            r["dn_int"],
            q_down,
        )

    # Colors: primary_contrast ADD (up) / Other (down); neutral gray for background.
    # (read from registry; NOT persisted — these categories already exist.)
    c_up = get_color("primary_contrast", "ADD", registry_path=args.registry,
                     persist=False)
    c_down = get_color("primary_contrast", "Other", registry_path=args.registry,
                       persist=False)
    c_bg = BACKGROUND_COLOR
    LOGGER.info("colors: up(ADD)=%s down(Other)=%s background=%s", c_up, c_down, c_bg)

    names = [r["name"] for r in records]
    up_fr = np.array([r["up_frac"] for r in records])
    dn_fr = np.array([r["dn_frac"] for r in records])
    bg_fr = np.array([r["bg_frac"] for r in records])

    y = np.arange(len(records))[::-1]  # top item = highest background fraction
    h = 0.26

    with publication_style():
        fig, ax = plt.subplots(figsize=(8.2, 6.4))
        # Background reference drawn first, hatched, so it reads clearly as the ref.
        ax.barh(y + h, bg_fr, height=h, color=c_bg, edgecolor="#555555",
                linewidth=0.6, hatch="///", zorder=2)
        ax.barh(y, up_fr, height=h, color=c_up, edgecolor="black", linewidth=0.4,
                zorder=2)
        ax.barh(y - h, dn_fr, height=h, color=c_down, edgecolor="black",
                linewidth=0.4, zorder=2)

        ax.set_yticks(y)
        ax.set_yticklabels(names)
        ax.set_xlabel("Fraction of proteins annotated to compartment")
        ax.set_xlim(0, max(bg_fr.max(), up_fr.max(), dn_fr.max()) * 1.12)
        ax.margins(y=0.02)
        ax.tick_params(axis="y", length=0)
        ax.xaxis.grid(True, linewidth=0.4, color="#dddddd", zorder=0)
        ax.set_axisbelow(True)

        fig.subplots_adjust(top=0.85, left=0.28)
        fig.suptitle(TITLE, y=0.985, fontsize=12.5, weight="bold")
        fig.text(0.5, 0.905, CAVEAT, ha="center", va="top", fontsize=9,
                 color="#444444", style="italic")

        # Separate legend image (swatch key). Denominators live in the labels so the
        # bar heights are self-documenting off-canvas.
        legend_fig = plt.figure(figsize=(4.4, 1.15))
        handles = [
            Patch(facecolor=c_up, edgecolor="black", linewidth=0.4,
                  label=f"up in ADD  (query, {q_up} of {SUBMITTED_N_UP} mapped)"),
            Patch(facecolor=c_down, edgecolor="black", linewidth=0.4,
                  label=f"down in ADD  (query, {q_down} of {SUBMITTED_N_DOWN} mapped)"),
            Patch(facecolor=c_bg, edgecolor="#555555", linewidth=0.6, hatch="///",
                  label=f"detected proteome background  (n={dom} eff. domain)"),
        ]
        legend_fig.legend(handles=handles, loc="center", frameon=False,
                          fontsize=9, title="fraction of ...", title_fontsize=9,
                          alignment="left")

    artifacts = save_figure(
        fig, args.output_dir, args.base_name, legend_fig=legend_fig
    )
    LOGGER.info("wrote %s", artifacts.svg)
    LOGGER.info("wrote %s", artifacts.png)
    LOGGER.info("wrote %s", artifacts.legend_svg)
    LOGGER.info("wrote %s", artifacts.legend_png)


if __name__ == "__main__":
    main()
