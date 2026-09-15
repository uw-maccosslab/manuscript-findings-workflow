"""Stage-3 integrity-gate loader for the plasma proteomics dataset (this project).

Reconciles and loads ``data/protein-data.csv`` (a wide, comma-separated feature x
sample matrix whose FIRST column header is the empty string ``""`` and holds the
UniProt feature id) against ``data/metadata.csv`` (48 sample rows, of which 40 are
experimental and 8 are metadata-only pooled-reference controls), returning the
project's standard :class:`~common.data_loading.Dataset` container (seeded from the
``wide-data-loader`` lib/ template — see ``scripts/common/data_loading.py``).

Two obligations, kept separate so verification can reconcile the full load against
the raw source before any feature is dropped (conventions/correctness.md, Stage-3
integrity gate):

  1. :func:`load_paired_dataset` — the FAITHFUL full load. Classifies every metadata
     row as experimental or control via an EXPLICIT allow-list on ``Condition``
     (failing loud on any unrecognized level), reconciles the classification BOTH
     ways against the data file's sample columns (an orphan on either side is an
     error), materializes two identifier-preserving intermediates under
     ``results/prepared/`` (``data/`` is read-only), and calls
     :func:`common.data_loading.load_wide_data` on the 40-experimental subset.
     Returns a ``(40 samples x 2343 features)`` ``Dataset`` on the recorded
     ``scale="log2"``, ordered by acquisition (``Run Order``).
  2. :func:`load_analysis_subset` — the Stage-4 ANALYSIS TRIMMING. Drops exactly the
     three features the scientist decided to exclude (state/DATA_DESCRIPTION.md):
     the yeast spike-in ``sp|P00924|ENO1_YEAST`` and the two >30%-floored
     (not-detected sentinel ``0.373``) features ``sp|Q9NX62|IMPA3_HUMAN`` and
     ``sp|Q5BLP8|CD048_HUMAN``. Fails loud if any target id is absent (so a data
     change can't silently skip a planned drop). Never drops a sample.

Run as a script for a manual sanity check (loads + logs shapes; writes nothing beyond
the ``results/prepared/`` intermediates)::

    ./.venv/bin/python scripts/scratch/project_loader.py
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# `common.data_loading` lives at scripts/common/data_loading.py; make the `scripts/`
# directory importable regardless of how this file is invoked (direct script, `-m`,
# or imported by pytest) — the same pattern used by
# scripts/figures/plot_group_x_cohort_confound.py.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.data_loading import Dataset as Dataset
from common.data_loading import Scale as Scale
from common.data_loading import load_wide_data as load_wide_data

__all__ = [
    "ANALYSIS_DROP_FEATURE_IDS",
    "CONTROL_CONDITIONS",
    "EXPERIMENTAL_CONDITIONS",
    "KNOWN_CONDITIONS",
    "Dataset",
    "Scale",
    "classify_conditions",
    "load_analysis_subset",
    "load_paired_dataset",
    "load_wide_data",
    "main",
    "prepare_feature_id_header",
    "read_data_sample_columns",
    "read_full_metadata",
    "reconcile_sample_partition",
]

__script_meta__: dict[str, object] = {
    "task": None,
    "kind": "module",
    "provides": [
        "EXPERIMENTAL_CONDITIONS",
        "CONTROL_CONDITIONS",
        "KNOWN_CONDITIONS",
        "ANALYSIS_DROP_FEATURE_IDS",
        "read_full_metadata",
        "classify_conditions",
        "read_data_sample_columns",
        "reconcile_sample_partition",
        "prepare_feature_id_header",
        "load_paired_dataset",
        "load_analysis_subset",
    ],
    "uses": ["common.data_loading"],
    "seeded_from": {"template": "wide-data-loader", "version": "0.5"},
    "description": (
        "Stage-3 integrity-gate loader for the plasma proteomics protein matrix: "
        "reconciles the 40-experimental/8-control Condition-based split against the "
        "data file's sample columns (both directions, fail-loud on any orphan), "
        "materializes an identifier-preserving intermediate for the empty-header "
        "feature-id column under results/prepared/, and loads the paired 40x2343 "
        "log2 Dataset via common.data_loading.load_wide_data. A second function "
        "returns the 2340-feature Stage-4 analysis subset (drops the yeast spike-in "
        "+ two >30%-floored sentinel features)."
    ),
}

logger = logging.getLogger(__name__)

# --- Project-specific constants (this study's Condition vocabulary; see
# state/METADATA.md) -----------------------------------------------------------

JOIN_KEY = "SampleID"
CONDITION_COLUMN = "Condition"
FEATURE_ID_COLUMN = "protein"
ORDER_BY_COLUMN = "Run Order"

#: The four experimental diagnostic groups (state/METADATA.md). Explicit allow-list —
#: NOT inferred — so a new or typo'd Condition value fails loud instead of silently
#: falling through to "control" or being dropped.
EXPERIMENTAL_CONDITIONS: frozenset[str] = frozenset({"ADD", "PDD", "PDCN", "HCN"})

#: The two pooled-reference control conditions (metadata-only; no data column).
CONTROL_CONDITIONS: frozenset[str] = frozenset(
    {"Plasma Reference", "TPAD plasma pooled reference"}
)

#: The only six Condition levels this loader knows how to classify.
KNOWN_CONDITIONS: frozenset[str] = EXPERIMENTAL_CONDITIONS | CONTROL_CONDITIONS

EXPECTED_N_EXPERIMENTAL = 40
EXPECTED_N_CONTROL = 8
EXPECTED_N_FEATURES = 2343

#: Stage-4 analysis-subset drops (state/DATA_DESCRIPTION.md, "Preprocessing
#: decisions confirmed with scientist, 2026-09-15"): the yeast spike-in (a loading
#: control, not human biology) and the two features carrying the not-detected
#: sentinel value 0.373 in >30% of samples.
ANALYSIS_DROP_FEATURE_IDS: tuple[str, ...] = (
    "sp|P00924|ENO1_YEAST",
    "sp|Q9NX62|IMPA3_HUMAN",
    "sp|Q5BLP8|CD048_HUMAN",
)

DEFAULT_DATA_FILE = Path("data/protein-data.csv")
DEFAULT_METADATA_FILE = Path("data/metadata.csv")
DEFAULT_PREPARED_DIR = Path("results/prepared")


# ---------------------------------------------------------------------------
# Metadata classification + reconciliation
# ---------------------------------------------------------------------------


def read_full_metadata(metadata_file: Path, *, sep: str = ",") -> pd.DataFrame:
    """Read the FULL metadata table (48 rows for this study), every column as text.

    No numeric coercion happens here — this table feeds only the experimental/
    control classification and reconciliation below; the typed, ordered load of
    the experimental subset happens downstream in
    :func:`common.data_loading.load_wide_data`.

    Fails loud if ``join_key``/``Condition`` are absent, ``join_key`` has
    duplicates, or any ``join_key`` cell is missing/NA (including a ``SampleID``
    that happens to spell a pandas NA token such as ``"None"`` or ``"NA"`` —
    ``dtype=str`` stops type *coercion* but not NA *detection*, so such a value
    would otherwise silently become ``NaN`` and surface as a confusing "orphan"
    downstream instead of a clear error here).
    """
    metadata_file = Path(metadata_file)
    meta = pd.read_csv(metadata_file, sep=sep, dtype=str)
    for col in (JOIN_KEY, CONDITION_COLUMN):
        if col not in meta.columns:
            raise ValueError(
                f"Metadata file {metadata_file} is missing required column {col!r}"
            )
    blank_join_key = meta[JOIN_KEY].isna()
    if bool(blank_join_key.any()):
        rows = (meta.index[blank_join_key][:5] + 2).tolist()  # +2: 1-based, past header
        raise ValueError(
            f"Missing/NA-token {JOIN_KEY!r} value(s) in {metadata_file} at row(s) "
            f"{rows} (of {int(blank_join_key.sum())} total). Note a {JOIN_KEY} "
            f"spelling an NA token (NA, NULL, None, N/A) reads as missing — rename "
            f"it in the source file if it is a real id."
        )
    duplicated = meta[JOIN_KEY].duplicated(keep=False)
    if bool(duplicated.any()):
        dup_values = sorted(meta.loc[duplicated, JOIN_KEY].unique())
        raise ValueError(
            f"Duplicate {JOIN_KEY!r} values in {metadata_file}: {dup_values}"
        )
    logger.info("read %d metadata row(s) from %s", len(meta), metadata_file)
    return meta


def classify_conditions(metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split ``metadata`` into ``(experimental, control)`` via an EXPLICIT allow-list.

    FAILS LOUD if any ``Condition`` value falls outside the six known levels
    (:data:`EXPERIMENTAL_CONDITIONS` | :data:`CONTROL_CONDITIONS`) — guards against a
    new or typo'd group being silently mis-classified or dropped.
    """
    if CONDITION_COLUMN not in metadata.columns:
        raise ValueError(f"metadata is missing required column {CONDITION_COLUMN!r}")
    conditions = metadata[CONDITION_COLUMN]
    unknown = sorted(set(conditions.dropna().tolist()) - KNOWN_CONDITIONS)
    if conditions.isna().any() or unknown:
        if conditions.isna().any():
            unknown = [*unknown, "<NA>"]
        raise ValueError(
            f"Unknown {CONDITION_COLUMN!r} value(s) {unknown} in metadata; known "
            f"experimental levels are {sorted(EXPERIMENTAL_CONDITIONS)}, known "
            f"control levels are {sorted(CONTROL_CONDITIONS)}. Refusing to guess how "
            f"to classify an unrecognized group."
        )
    experimental = metadata[conditions.isin(EXPERIMENTAL_CONDITIONS)].reset_index(
        drop=True
    )
    control = metadata[conditions.isin(CONTROL_CONDITIONS)].reset_index(drop=True)
    logger.info(
        "classified metadata: %d experimental, %d control (of %d total rows)",
        len(experimental),
        len(control),
        len(metadata),
    )
    return experimental, control


def reconcile_sample_partition(
    experimental: pd.DataFrame,
    control: pd.DataFrame,
    data_sample_columns: list[str],
    *,
    data_file: Path,
) -> None:
    """Reconcile the classified metadata partition against the data file's sample
    columns, in BOTH directions. FAILS LOUD on any orphan:

      * an experimental ``SampleID`` with no matching data column (a true orphan —
        NOT a control, so it cannot legitimately be metadata-only);
      * a data column that is not an experimental ``SampleID`` (an unexpected/
        unclassified column, INCLUDING a control that leaked a data column — a
        control id is by construction not an experimental id, so it is caught here
        with the same "does not match any experimental SampleID" diagnosis).

    Given both checks pass, every metadata row absent from the data is necessarily
    exactly the control set (the third documented reconciliation direction) — this
    is asserted defensively below rather than trusted as an invariant.
    """
    data_cols = set(data_sample_columns)
    exp_ids = set(experimental[JOIN_KEY].tolist())
    ctrl_ids = set(control[JOIN_KEY].tolist())

    exp_missing_from_data = sorted(exp_ids - data_cols)
    if exp_missing_from_data:
        raise ValueError(
            f"{len(exp_missing_from_data)} experimental SampleID(s) are in metadata "
            f"but have no matching column in {data_file} (true orphan(s) -- not a "
            f"known control): {exp_missing_from_data}"
        )
    data_not_experimental = sorted(data_cols - exp_ids)
    if data_not_experimental:
        raise ValueError(
            f"{len(data_not_experimental)} column(s) of {data_file} do not match any "
            f"experimental SampleID in metadata (a control SampleID is expected to "
            f"be metadata-only, so it is never a valid data column either): "
            f"{data_not_experimental}"
        )
    absent_from_data = (exp_ids | ctrl_ids) - data_cols
    unclassified_absent = sorted(absent_from_data - ctrl_ids)
    if unclassified_absent:  # pragma: no cover - implied unreachable by the checks
        # above (exp_ids is already a subset of data_cols at this point), kept as a
        # defensive fail-loud backstop rather than a silent invariant.
        raise ValueError(
            f"{len(unclassified_absent)} metadata SampleID(s) are absent from "
            f"{data_file} but are not classified as control: {unclassified_absent}"
        )
    logger.info(
        "reconciled %d experimental SampleID(s) against %d data column(s) (exact "
        "match); %d control SampleID(s) confirmed metadata-only",
        len(exp_ids),
        len(data_cols),
        len(ctrl_ids),
    )


# ---------------------------------------------------------------------------
# Empty-header feature-id column
# ---------------------------------------------------------------------------


def read_data_sample_columns(data_file: Path, *, sep: str = ",") -> list[str]:
    """Return the data file's sample-column headers (every column but the first).

    Verifies the first header field is the empty string ``""`` (this study's
    feature-id column carries no header) — fails loud if the format has changed,
    rather than silently treating some other column as the id column.

    Parses the raw header line with :mod:`csv` rather than ``pandas.read_csv``:
    pandas silently renames an empty/unlabeled header to ``"Unnamed: 0"`` (so
    checking for the literal empty string against a pandas-parsed header would
    never see it), and the manual parse also correctly handles both this study's
    quoted-empty (``""``) convention and a plain unquoted-empty one.
    """
    data_file = Path(data_file)
    with data_file.open("r", newline="") as fh:
        reader = csv.reader(fh, delimiter=sep)
        try:
            raw_columns = next(reader)
        except StopIteration as exc:
            raise ValueError(f"{data_file} is empty (no header row).") from exc
    if not raw_columns or raw_columns[0] != "":
        raise ValueError(
            f"Expected {data_file} to have an empty ('') header on its first "
            f"column (the feature-id column); got {raw_columns[:1]!r}. Data format "
            f"may have changed."
        )
    sample_columns = [str(c) for c in raw_columns[1:]]
    if not sample_columns:
        raise ValueError(f"{data_file} has no sample columns after the id column.")
    return sample_columns


def prepare_feature_id_header(
    raw_data_file: Path,
    prepared_data_file: Path,
    *,
    data_sep: str = ",",
    feature_id_column: str = FEATURE_ID_COLUMN,
) -> Path:
    """Materialize a byte-faithful copy of ``raw_data_file`` with its empty first
    header field (the feature-id column) renamed to ``feature_id_column``.

    ``data/`` is read-only (the project's ``guard_readonly_data`` hook), and
    :func:`common.data_loading.load_wide_data` needs a *named* id column, but this
    study's data file ships an EMPTY header (``""``) on the leading column. This
    splices ONLY the raw text of that first header field out of the header line and
    replaces it with the bare ``feature_id_column`` name (no delimiter/special
    character, so it never needs CSV quoting) — every other byte of the header
    line, and every byte of every subsequent line (every id, every abundance value),
    is streamed through completely unmodified. No value can be reformatted,
    re-typed, or truncated by this step; it is a pure rename. Written under
    ``results/prepared/`` (never into ``data/``).

    Raises if the first header field is not exactly the expected empty value —
    refuses to guess/blindly rename a header that has actually changed.
    """
    if data_sep in feature_id_column or '"' in feature_id_column:
        raise ValueError(
            f"feature_id_column {feature_id_column!r} contains the separator or a "
            f"quote character and would need CSV quoting; pass a plain identifier."
        )
    raw_data_file = Path(raw_data_file)
    prepared_data_file = Path(prepared_data_file)
    with raw_data_file.open("r", newline="") as fh:
        first_line = fh.readline()
        if not first_line:
            raise ValueError(f"{raw_data_file} is empty (no header row).")
        rest = fh.read()

    line_ending = "\r\n" if first_line.endswith("\r\n") else "\n"
    has_ending = first_line.endswith(line_ending)
    body = first_line[: -len(line_ending)] if has_ending else first_line

    empty_quoted = f'""{data_sep}'
    if body.startswith(empty_quoted):
        remainder = body[len(empty_quoted) - len(data_sep) :]  # keep the leading sep
    elif body.startswith(data_sep):
        remainder = body
    else:
        raise ValueError(
            f"Could not find an empty first header field in {raw_data_file}'s header "
            f"line (expected it to start with {empty_quoted!r} or {data_sep!r}); got "
            f"{body[:20]!r}. Refusing to guess how to rename it."
        )
    new_header_line = feature_id_column + remainder
    if has_ending:
        new_header_line += line_ending

    prepared_data_file.parent.mkdir(parents=True, exist_ok=True)
    with prepared_data_file.open("w", newline="") as out:
        out.write(new_header_line)
        out.write(rest)
    logger.info(
        "prepared %s -> %s (renamed the empty first header field to %r; all other "
        "bytes copied unchanged)",
        raw_data_file,
        prepared_data_file,
        feature_id_column,
    )
    return prepared_data_file


# ---------------------------------------------------------------------------
# Public loaders
# ---------------------------------------------------------------------------


def load_paired_dataset(
    data_file: Path = DEFAULT_DATA_FILE,
    metadata_file: Path = DEFAULT_METADATA_FILE,
    *,
    prepared_dir: Path = DEFAULT_PREPARED_DIR,
    data_sep: str = ",",
    metadata_sep: str = ",",
    order_by: str | None = ORDER_BY_COLUMN,
    expected_n_samples: int = EXPECTED_N_EXPERIMENTAL,
    expected_n_features: int = EXPECTED_N_FEATURES,
    expected_n_control: int | None = None,
    scale: Scale = "log2",
) -> Dataset:
    """Reconcile + load the paired experimental protein Dataset (Stage-3 loader).

    1. Reads the full metadata and classifies experimental (40) vs control (8) via
       an explicit ``Condition`` allow-list (:func:`classify_conditions`).
    2. Reconciles that classification BOTH ways against the data file's sample
       columns (:func:`reconcile_sample_partition`) — an orphan on either side is a
       hard error.
    3. Materializes two identifier-preserving intermediates under ``prepared_dir``
       (``results/prepared/`` by default; ``data/`` is read-only): the data matrix
       with its empty feature-id header renamed to ``"protein"``
       (:func:`prepare_feature_id_header`), and the 40-row experimental-only
       metadata subset — so the pairing check inside
       :func:`common.data_loading.load_wide_data` certifies an exact 40<->40 match.
    4. Loads via :func:`common.data_loading.load_wide_data` with ``scale="log2"``
       (recorded, not inferred), ``missing="error"`` (the matrix has no NaN),
       ``require_unique_features=True``, and the expected-shape sanity checks,
       ordered by acquisition (``Run Order``, coerced numeric).
    5. Optionally re-checks the classified control count (``expected_n_control``,
       ``None`` by default — the reconciliation in step 2 already guarantees every
       metadata row absent from the data is exactly the classified control set;
       this is an additional opt-in numeric sanity check, e.g. ``8`` for the real
       project files).

    Returns
    -------
    Dataset
        ``(expected_n_samples, expected_n_features)`` abundances, ``scale="log2"``.
    """
    data_file = Path(data_file)
    metadata_file = Path(metadata_file)
    prepared_dir = Path(prepared_dir)

    logger.info(
        "load_paired_dataset: data=%s metadata=%s scale=%s order_by=%s",
        data_file,
        metadata_file,
        scale,
        order_by,
    )

    full_metadata = read_full_metadata(metadata_file, sep=metadata_sep)
    experimental, control = classify_conditions(full_metadata)

    data_sample_columns = read_data_sample_columns(data_file, sep=data_sep)
    logger.info(
        "data file %s: %d sample column(s)", data_file, len(data_sample_columns)
    )
    reconcile_sample_partition(
        experimental, control, data_sample_columns, data_file=data_file
    )

    prepared_dir.mkdir(parents=True, exist_ok=True)
    prepared_data_file = prepare_feature_id_header(
        data_file,
        prepared_dir / f"{data_file.stem}.protein-id{data_file.suffix}",
        data_sep=data_sep,
    )
    prepared_metadata_file = (
        prepared_dir / f"{metadata_file.stem}.experimental{metadata_file.suffix}"
    )
    experimental.to_csv(prepared_metadata_file, sep=metadata_sep, index=False)
    logger.info(
        "wrote experimental-subset metadata (%d row(s)) to %s",
        len(experimental),
        prepared_metadata_file,
    )

    numeric_columns = (order_by,) if order_by is not None else ()
    dataset = load_wide_data(
        prepared_data_file,
        prepared_metadata_file,
        join_key=JOIN_KEY,
        id_columns=(FEATURE_ID_COLUMN,),
        feature_id_column=FEATURE_ID_COLUMN,
        numeric_columns=numeric_columns,
        order_by=order_by,
        missing="error",
        require_unique_features=True,
        expected_n_samples=expected_n_samples,
        expected_n_features=expected_n_features,
        data_sep=data_sep,
        metadata_sep=metadata_sep,
        scale=scale,
    )
    logger.info(
        "loaded paired dataset: %d sample(s) x %d feature(s), scale=%s",
        dataset.abundances.shape[0],
        dataset.abundances.shape[1],
        dataset.scale,
    )
    if expected_n_control is not None and len(control) != expected_n_control:
        raise ValueError(
            f"Expected {expected_n_control} control SampleID(s) but classified "
            f"{len(control)}: {sorted(control[JOIN_KEY])}"
        )
    return dataset


def load_analysis_subset(
    dataset: Dataset,
    *,
    drop_feature_ids: tuple[str, ...] = ANALYSIS_DROP_FEATURE_IDS,
) -> Dataset:
    """Return the Stage-4 ANALYSIS SUBSET: ``dataset`` with ``drop_feature_ids``
    removed (the yeast spike-in + the two >30%-floored sentinel features by default).

    Samples/metadata are never touched — only features are dropped. FAILS LOUD if
    ``dataset.scale`` is not ``"log2"``, or if any id in ``drop_feature_ids`` is not
    present in ``dataset.feature_names`` (so a data change can't silently skip a
    planned drop).
    """
    if dataset.scale != "log2":
        raise ValueError(
            f"load_analysis_subset expects a log2-scale Dataset; got scale="
            f"{dataset.scale!r}."
        )
    feature_names = dataset.feature_names
    present = set(feature_names.tolist())
    missing = [fid for fid in drop_feature_ids if fid not in present]
    if missing:
        raise ValueError(
            f"Target feature id(s) to drop are not present in the dataset: "
            f"{missing}. Refusing to silently skip a planned drop -- confirm the "
            f"id(s) still match the current data (state/DATA_DESCRIPTION.md)."
        )

    drop_set = set(drop_feature_ids)
    keep_mask = np.array([name not in drop_set for name in feature_names], dtype=bool)
    n_dropped = int((~keep_mask).sum())
    if n_dropped != len(drop_set):
        raise ValueError(
            f"Expected to drop exactly {len(drop_set)} feature(s) but the mask "
            f"removes {n_dropped}; feature ids may not be unique."
        )

    trimmed = Dataset(
        abundances=dataset.abundances[:, keep_mask],
        feature_names=feature_names[keep_mask],
        feature_metadata=dataset.feature_metadata.loc[keep_mask].reset_index(drop=True),
        metadata=dataset.metadata,
        scale=dataset.scale,
    )
    if trimmed.abundances.shape[0] != dataset.abundances.shape[0]:
        raise ValueError(  # pragma: no cover - defensive; keep_mask never touches rows
            "load_analysis_subset must not change the sample count."
        )
    logger.info(
        "analysis subset: dropped %d feature(s) %s -> %d sample(s) x %d feature(s)",
        n_dropped,
        sorted(drop_set),
        trimmed.abundances.shape[0],
        trimmed.abundances.shape[1],
    )
    return trimmed


# ---------------------------------------------------------------------------
# Manual sanity-check entry point
# ---------------------------------------------------------------------------


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA_FILE)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA_FILE)
    parser.add_argument("--prepared-dir", type=Path, default=DEFAULT_PREPARED_DIR)
    args = parser.parse_args()

    full = load_paired_dataset(
        args.data,
        args.metadata,
        prepared_dir=args.prepared_dir,
        expected_n_control=EXPECTED_N_CONTROL,
    )
    logger.info("full dataset: %d x %d, scale=%s", *full.abundances.shape, full.scale)
    trimmed = load_analysis_subset(full)
    logger.info(
        "analysis subset: %d x %d, scale=%s",
        *trimmed.abundances.shape,
        trimmed.scale,
    )


if __name__ == "__main__":
    main()
