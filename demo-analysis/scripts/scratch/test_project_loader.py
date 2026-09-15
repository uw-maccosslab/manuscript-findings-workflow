"""Tests for ``scripts/scratch/project_loader.py`` (the Stage-3 integrity-gate loader).

Both obligations of conventions/correctness.md's "Data loading is two obligations":

  A. Test the loader — unit tests on hand-built fixtures, property/invariant tests,
     a planted-truth fixture, and edge cases (below).
  B. Verify the loaded data on the real file — ``test_real_data_reconciliation``
     loads ``data/protein-data.csv`` + ``data/metadata.csv`` and checks shape, exact
     experimental/control counts, and a spot-checked raw value.

Run:
    ./.venv/bin/python -m pytest scripts/scratch/test_project_loader.py -v
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from hypothesis import given
from hypothesis import strategies as st
from project_loader import (
    ANALYSIS_DROP_FEATURE_IDS,
    Dataset,
    classify_conditions,
    load_analysis_subset,
    load_paired_dataset,
    prepare_feature_id_header,
    read_data_sample_columns,
    read_full_metadata,
    reconcile_sample_partition,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
REAL_DATA_FILE = PROJECT_ROOT / "data" / "protein-data.csv"
REAL_METADATA_FILE = PROJECT_ROOT / "data" / "metadata.csv"

# --------------------------------------------------------------------------
# Fixture builders
# --------------------------------------------------------------------------


def _write_metadata_csv(path: Path, rows: list[dict[str, str]]) -> None:
    """Write a metadata CSV from a list of row-dicts (column union across rows)."""
    pd.DataFrame(rows).to_csv(path, index=False)


def _write_data_csv(
    path: Path,
    feature_ids: Sequence[str],
    sample_ids: Sequence[str],
    matrix: np.ndarray,
) -> None:
    """Write a wide data CSV: an UNQUOTED-empty first header + one column per sample.

    ``matrix`` has shape ``(n_features, n_samples)``.
    """
    df = pd.DataFrame(matrix, columns=list(sample_ids))
    df.insert(0, "", list(feature_ids))
    df.to_csv(path, index=False)


def _small_fixture(
    tmp_path: Path,
) -> tuple[Path, Path, list[str], list[str], np.ndarray]:
    """A small, hand-verified 5-experimental / 2-control fixture.

    Feature ids deliberately include a pipe-delimited UniProt-style id and a
    date-like id (the classic spreadsheet-corruption trap) to prove no coercion.
    Sample columns are written in a SCRAMBLED (non-Run-Order) order to prove the
    loader reorders by ``Run Order``, not by file column order.
    """
    feature_ids = ["sp|P00001|TEST_HUMAN", "2-Mar", "sp|Q99999|OTHR_HUMAN", "1/2/2023"]
    # (sample_id, condition, run_order)
    experimental = [
        ("TPAD1001", "ADD", "30"),
        ("TPAD1002", "PDD", "10"),
        ("TPAD1003", "PDCN", "50"),
        ("TPAD1004", "HCN", "20"),
        ("TPAD1005", "ADD", "40"),
    ]
    control = [
        ("CP01", "Plasma Reference", ""),
        ("KP01", "TPAD plasma pooled reference", ""),
    ]

    metadata_rows = [
        {"SampleID": sid, "Condition": cond, "Run Order": run}
        for sid, cond, run in [*experimental, *control]
    ]
    _write_metadata_csv(tmp_path / "metadata.csv", metadata_rows)

    rng = np.random.default_rng(20260915)
    # File column order intentionally NOT sorted by Run Order (scrambled).
    scrambled_order = [4, 1, 3, 0, 2]  # index into `experimental`
    sample_ids = [experimental[i][0] for i in scrambled_order]
    matrix = rng.uniform(15.0, 30.0, size=(len(feature_ids), len(sample_ids)))
    _write_data_csv(tmp_path / "protein-data.csv", feature_ids, sample_ids, matrix)

    return (
        tmp_path / "protein-data.csv",
        tmp_path / "metadata.csv",
        feature_ids,
        sample_ids,
        matrix,
    )


# --------------------------------------------------------------------------
# Unit tests: hand-built fixture — orientation, pairing, scale, id fidelity
# --------------------------------------------------------------------------


def test_load_paired_dataset_hand_built_fixture(tmp_path: Path) -> None:
    data_file, metadata_file, feature_ids, sample_ids, matrix = _small_fixture(tmp_path)

    dataset = load_paired_dataset(
        data_file,
        metadata_file,
        prepared_dir=tmp_path / "prepared",
        expected_n_samples=5,
        expected_n_features=4,
    )

    # Orientation: rows=samples, cols=features.
    assert dataset.abundances.shape == (5, 4)
    # Scale recorded (not inferred).
    assert dataset.scale == "log2"
    # Feature ids preserved verbatim, in source column order, including the
    # pipe-delimited and date-like ids (no type coercion).
    assert list(dataset.feature_names) == feature_ids
    assert "2-Mar" in dataset.feature_names
    assert "1/2/2023" in dataset.feature_names
    assert "sp|P00001|TEST_HUMAN" in dataset.feature_names

    # Pairing: exactly the 5 experimental SampleIDs, indexed by SampleID.
    assert set(dataset.metadata.index) == {
        "TPAD1001",
        "TPAD1002",
        "TPAD1003",
        "TPAD1004",
        "TPAD1005",
    }
    assert dataset.metadata.index.name == "SampleID"

    # Sample order follows ascending Run Order, not file column order.
    assert dataset.metadata["Run Order"].tolist() == [10, 20, 30, 40, 50]
    assert dataset.metadata.index.tolist() == [
        "TPAD1002",
        "TPAD1004",
        "TPAD1001",
        "TPAD1005",
        "TPAD1003",
    ]

    # Every loaded value equals its source cell (spot check via source lookup).
    for out_row, sid in enumerate(dataset.metadata.index):
        src_col = sample_ids.index(sid)
        for feat_idx in range(4):
            assert dataset.abundances[out_row, feat_idx] == pytest.approx(
                matrix[feat_idx, src_col]
            )


def test_load_analysis_subset_drops_exactly_three(tmp_path: Path) -> None:
    feature_ids = [
        "sp|P00924|ENO1_YEAST",
        "sp|Q9NX62|IMPA3_HUMAN",
        "sp|Q5BLP8|CD048_HUMAN",
        "sp|P99999|KEEP_HUMAN",
    ]
    sample_ids = ["TPAD1001", "TPAD1002"]
    metadata_rows = [
        {"SampleID": "TPAD1001", "Condition": "ADD", "Run Order": "1"},
        {"SampleID": "TPAD1002", "Condition": "HCN", "Run Order": "2"},
    ]
    _write_metadata_csv(tmp_path / "metadata.csv", metadata_rows)
    matrix = np.full((4, 2), 20.0)
    _write_data_csv(tmp_path / "protein-data.csv", feature_ids, sample_ids, matrix)

    dataset = load_paired_dataset(
        tmp_path / "protein-data.csv",
        tmp_path / "metadata.csv",
        prepared_dir=tmp_path / "prepared",
        expected_n_samples=2,
        expected_n_features=4,
    )
    trimmed = load_analysis_subset(dataset)
    assert trimmed.abundances.shape == (2, 1)
    assert list(trimmed.feature_names) == ["sp|P99999|KEEP_HUMAN"]
    # Samples untouched.
    assert trimmed.metadata.index.tolist() == dataset.metadata.index.tolist()
    assert trimmed.scale == "log2"


# --------------------------------------------------------------------------
# prepare_feature_id_header — both quoting conventions
# --------------------------------------------------------------------------


def test_prepare_feature_id_header_quoted_empty(tmp_path: Path) -> None:
    raw = tmp_path / "raw.csv"
    raw.write_text('"","A","B"\n"feat|1",1.0,2.0\n"feat|2",3.0,4.0\n')
    prepared = tmp_path / "prepared" / "raw.csv"
    result = prepare_feature_id_header(raw, prepared, data_sep=",")
    text = result.read_text()
    assert text.startswith('protein,"A","B"\n')
    # Every other byte carried through unchanged.
    assert text == 'protein,"A","B"\n"feat|1",1.0,2.0\n"feat|2",3.0,4.0\n'


def test_prepare_feature_id_header_unquoted_empty(tmp_path: Path) -> None:
    raw = tmp_path / "raw.csv"
    raw.write_text(",A,B\nfeat1,1.0,2.0\n")
    prepared = tmp_path / "prepared" / "raw.csv"
    result = prepare_feature_id_header(raw, prepared, data_sep=",")
    assert result.read_text() == "protein,A,B\nfeat1,1.0,2.0\n"


def test_prepare_feature_id_header_rejects_nonempty_header(tmp_path: Path) -> None:
    raw = tmp_path / "raw.csv"
    raw.write_text('"already_named","A"\n"feat1",1.0\n')
    with pytest.raises(ValueError, match="empty"):
        prepare_feature_id_header(raw, tmp_path / "prepared" / "raw.csv")


# --------------------------------------------------------------------------
# Property / invariant tests (Hypothesis)
# --------------------------------------------------------------------------

# Tokens pandas' CSV reader treats as missing (NA) REGARDLESS of ``dtype=str`` — an
# id strategy must exclude these, or a generated id would silently become NaN on
# read instead of round-tripping (a correctly-handled but different failure mode
# from the one this property test targets; exercised deliberately elsewhere by
# test_read_full_metadata_na_token_sample_id_fails_loud /
# common.data_loading's own ``_check_ids_present`` tests).
_PANDAS_NA_TOKENS = frozenset(
    {
        "",
        "#N/A",
        "#N/A N/A",
        "#NA",
        "-1.#IND",
        "-1.#QNAN",
        "-NaN",
        "-nan",
        "1.#IND",
        "1.#QNAN",
        "<NA>",
        "N/A",
        "NA",
        "NULL",
        "NaN",
        "None",
        "n/a",
        "nan",
        "null",
    }
)
_ID_ALPHABET = st.characters(
    whitelist_categories=("Lu", "Ll", "Nd"), whitelist_characters="_|.-"
)
_id_strategy = st.text(alphabet=_ID_ALPHABET, min_size=1, max_size=10).filter(
    lambda s: s not in _PANDAS_NA_TOKENS
)
# Sample ids additionally must not collide with the reserved feature-id column name.
_sample_id_strategy = st.text(alphabet=_ID_ALPHABET, min_size=1, max_size=8).filter(
    lambda s: s not in _PANDAS_NA_TOKENS and s != "protein"
)
_value_strategy = st.floats(
    min_value=1.0, max_value=100.0, allow_nan=False, allow_infinity=False, width=64
)


@given(
    n_features=st.integers(min_value=1, max_value=5),
    n_samples=st.integers(min_value=2, max_value=5),
    data=st.data(),
)
def test_property_values_preserved_and_ordered(
    n_features: int,
    n_samples: int,
    data: st.DataObject,
    tmp_path_factory: pytest.TempPathFactory,
) -> None:
    feature_ids = data.draw(
        st.lists(_id_strategy, min_size=n_features, max_size=n_features, unique=True)
    )
    sample_ids = data.draw(
        st.lists(
            _sample_id_strategy,
            min_size=n_samples,
            max_size=n_samples,
            unique=True,
        )
    )
    run_orders = data.draw(st.permutations(list(range(1, n_samples + 1))))
    values = data.draw(
        st.lists(
            st.lists(_value_strategy, min_size=n_samples, max_size=n_samples),
            min_size=n_features,
            max_size=n_features,
        )
    )
    matrix = np.array(values, dtype=float)

    tmp_path = tmp_path_factory.mktemp("prop")
    metadata_rows = [
        {"SampleID": sid, "Condition": "ADD", "Run Order": str(run)}
        for sid, run in zip(sample_ids, run_orders, strict=True)
    ]
    _write_metadata_csv(tmp_path / "metadata.csv", metadata_rows)
    _write_data_csv(tmp_path / "protein-data.csv", feature_ids, sample_ids, matrix)

    dataset = load_paired_dataset(
        tmp_path / "protein-data.csv",
        tmp_path / "metadata.csv",
        prepared_dir=tmp_path / "prepared",
        expected_n_samples=n_samples,
        expected_n_features=n_features,
    )

    # Loaded value count == source non-id cells.
    assert dataset.abundances.size == n_features * n_samples

    # Sample order follows Run Order (ascending), not file column order.
    assert dataset.metadata["Run Order"].tolist() == sorted(run_orders)

    # Every loaded value equals the corresponding source cell (no invented values).
    for out_row, sid in enumerate(dataset.metadata.index):
        src_col = sample_ids.index(sid)
        for feat_idx in range(n_features):
            assert dataset.abundances[out_row, feat_idx] == pytest.approx(
                matrix[feat_idx, src_col]
            )
    assert list(dataset.feature_names) == feature_ids


# --------------------------------------------------------------------------
# Planted-truth fixture: a known mean difference must be recovered
# --------------------------------------------------------------------------


def test_planted_truth_mean_difference_recovered(tmp_path: Path) -> None:
    rng = np.random.default_rng(20260915)
    n_per_group = 6
    group_a_ids = [f"TPAD{1000 + i}" for i in range(n_per_group)]
    group_b_ids = [f"TPAD{1000 + n_per_group + i}" for i in range(n_per_group)]
    sample_ids = group_a_ids + group_b_ids
    conditions = ["ADD"] * n_per_group + ["PDD"] * n_per_group
    run_orders = rng.permutation(np.arange(1, 2 * n_per_group + 1))

    feature_ids = ["sp|P00001|EFFECT_HUMAN", "sp|P00002|NULL_HUMAN"]
    known_delta = 2.5
    base = 20.0
    noise = rng.normal(0.0, 0.01, size=(2, 2 * n_per_group))
    matrix = np.empty((2, 2 * n_per_group))
    matrix[0, :n_per_group] = base + noise[0, :n_per_group]
    matrix[0, n_per_group:] = base + known_delta + noise[0, n_per_group:]
    matrix[1, :] = base + noise[1, :]  # no planted effect

    metadata_rows = [
        {"SampleID": sid, "Condition": cond, "Run Order": str(run)}
        for sid, cond, run in zip(sample_ids, conditions, run_orders, strict=True)
    ]
    _write_metadata_csv(tmp_path / "metadata.csv", metadata_rows)
    _write_data_csv(tmp_path / "protein-data.csv", feature_ids, sample_ids, matrix)

    dataset = load_paired_dataset(
        tmp_path / "protein-data.csv",
        tmp_path / "metadata.csv",
        prepared_dir=tmp_path / "prepared",
        expected_n_samples=2 * n_per_group,
        expected_n_features=2,
    )

    is_b = (dataset.metadata["Condition"] == "PDD").to_numpy()
    effect_col = list(dataset.feature_names).index("sp|P00001|EFFECT_HUMAN")
    null_col = list(dataset.feature_names).index("sp|P00002|NULL_HUMAN")

    recovered_delta = (
        dataset.abundances[is_b, effect_col].mean()
        - dataset.abundances[~is_b, effect_col].mean()
    )
    recovered_null = (
        dataset.abundances[is_b, null_col].mean()
        - dataset.abundances[~is_b, null_col].mean()
    )
    assert recovered_delta == pytest.approx(known_delta, abs=0.05)
    assert recovered_null == pytest.approx(0.0, abs=0.05)


# --------------------------------------------------------------------------
# Edge cases — each must FAIL LOUD
# --------------------------------------------------------------------------


def test_true_orphan_experimental_sample_id_fails_loud(tmp_path: Path) -> None:
    """An experimental SampleID with no matching data column is a true orphan."""
    feature_ids = ["sp|P00001|TEST_HUMAN"]
    sample_ids = ["TPAD1001"]  # data file has only 1 column
    metadata_rows = [
        {"SampleID": "TPAD1001", "Condition": "ADD", "Run Order": "1"},
        {"SampleID": "TPAD1002", "Condition": "HCN", "Run Order": "2"},  # orphan
    ]
    _write_metadata_csv(tmp_path / "metadata.csv", metadata_rows)
    matrix = np.full((1, 1), 20.0)
    _write_data_csv(tmp_path / "protein-data.csv", feature_ids, sample_ids, matrix)

    with pytest.raises(ValueError, match="orphan"):
        load_paired_dataset(
            tmp_path / "protein-data.csv",
            tmp_path / "metadata.csv",
            prepared_dir=tmp_path / "prepared",
        )


def test_data_column_absent_from_metadata_fails_loud(tmp_path: Path) -> None:
    feature_ids = ["sp|P00001|TEST_HUMAN"]
    sample_ids = ["TPAD1001", "TPAD9999"]  # TPAD9999 is not in metadata at all
    metadata_rows = [{"SampleID": "TPAD1001", "Condition": "ADD", "Run Order": "1"}]
    _write_metadata_csv(tmp_path / "metadata.csv", metadata_rows)
    matrix = np.full((1, 2), 20.0)
    _write_data_csv(tmp_path / "protein-data.csv", feature_ids, sample_ids, matrix)

    with pytest.raises(ValueError, match="do not match any experimental SampleID"):
        load_paired_dataset(
            tmp_path / "protein-data.csv",
            tmp_path / "metadata.csv",
            prepared_dir=tmp_path / "prepared",
        )


def test_duplicate_feature_id_fails_loud(tmp_path: Path) -> None:
    feature_ids = ["sp|P00001|TEST_HUMAN", "sp|P00001|TEST_HUMAN"]
    sample_ids = ["TPAD1001"]
    metadata_rows = [{"SampleID": "TPAD1001", "Condition": "ADD", "Run Order": "1"}]
    _write_metadata_csv(tmp_path / "metadata.csv", metadata_rows)
    matrix = np.full((2, 1), 20.0)
    _write_data_csv(tmp_path / "protein-data.csv", feature_ids, sample_ids, matrix)

    with pytest.raises(ValueError, match="Duplicate feature ids"):
        load_paired_dataset(
            tmp_path / "protein-data.csv",
            tmp_path / "metadata.csv",
            prepared_dir=tmp_path / "prepared",
        )


def test_nan_abundance_fails_loud_under_missing_error(tmp_path: Path) -> None:
    metadata_rows = [
        {"SampleID": "TPAD1001", "Condition": "ADD", "Run Order": "1"},
        {"SampleID": "TPAD1002", "Condition": "HCN", "Run Order": "2"},
    ]
    _write_metadata_csv(tmp_path / "metadata.csv", metadata_rows)
    # Blank cell -> NaN on read.
    (tmp_path / "protein-data.csv").write_text(
        ',"TPAD1001","TPAD1002"\n"sp|P00001|TEST_HUMAN",20.0,\n'
    )

    with pytest.raises(ValueError, match="missing value"):
        load_paired_dataset(
            tmp_path / "protein-data.csv",
            tmp_path / "metadata.csv",
            prepared_dir=tmp_path / "prepared",
        )


def test_expected_n_control_mismatch_fails_loud(tmp_path: Path) -> None:
    """An opt-in ``expected_n_control`` check catches an unexpected control count."""
    feature_ids = ["sp|P00001|TEST_HUMAN"]
    sample_ids = ["TPAD1001"]
    metadata_rows = [
        {"SampleID": "TPAD1001", "Condition": "ADD", "Run Order": "1"},
        {"SampleID": "CP01", "Condition": "Plasma Reference", "Run Order": ""},
    ]
    _write_metadata_csv(tmp_path / "metadata.csv", metadata_rows)
    matrix = np.full((1, 1), 20.0)
    _write_data_csv(tmp_path / "protein-data.csv", feature_ids, sample_ids, matrix)

    with pytest.raises(ValueError, match="Expected 2 control SampleID"):
        load_paired_dataset(
            tmp_path / "protein-data.csv",
            tmp_path / "metadata.csv",
            prepared_dir=tmp_path / "prepared",
            expected_n_samples=1,
            expected_n_features=1,
            expected_n_control=2,
        )


def test_unknown_condition_level_fails_loud(tmp_path: Path) -> None:
    metadata_rows = [
        {"SampleID": "TPAD1001", "Condition": "ADD", "Run Order": "1"},
        {"SampleID": "TPAD1002", "Condition": "SomeNewGroup", "Run Order": "2"},
    ]
    metadata = pd.DataFrame(metadata_rows)
    with pytest.raises(ValueError, match="Unknown 'Condition' value"):
        classify_conditions(metadata)


def test_missing_target_feature_at_analysis_drop_fails_loud() -> None:
    dataset = Dataset(
        abundances=np.array([[20.0, 21.0]]),
        feature_names=np.array(["sp|P00924|ENO1_YEAST", "sp|P99999|KEEP_HUMAN"]),
        feature_metadata=pd.DataFrame(
            {"protein": ["sp|P00924|ENO1_YEAST", "sp|P99999|KEEP_HUMAN"]}
        ),
        metadata=pd.DataFrame({"Condition": ["ADD"]}, index=pd.Index(["TPAD1001"])),
        scale="log2",
    )
    # ANALYSIS_DROP_FEATURE_IDS names 3 ids; only 1 is present here -> the other 2
    # are "missing target features" and must fail loud, not be silently skipped.
    with pytest.raises(ValueError, match="not present in the dataset"):
        load_analysis_subset(dataset, drop_feature_ids=ANALYSIS_DROP_FEATURE_IDS)


# --------------------------------------------------------------------------
# Additional edge cases (empty / single-sample / ties)
# --------------------------------------------------------------------------


def test_empty_data_file_fails_loud(tmp_path: Path) -> None:
    metadata_rows = [{"SampleID": "TPAD1001", "Condition": "ADD", "Run Order": "1"}]
    _write_metadata_csv(tmp_path / "metadata.csv", metadata_rows)
    (tmp_path / "protein-data.csv").write_text(',"TPAD1001"\n')  # header only

    with pytest.raises(ValueError, match="no feature rows"):
        load_paired_dataset(
            tmp_path / "protein-data.csv",
            tmp_path / "metadata.csv",
            prepared_dir=tmp_path / "prepared",
        )


def test_single_sample_loads_successfully(tmp_path: Path) -> None:
    feature_ids = ["sp|P00001|TEST_HUMAN", "sp|P00002|OTHER_HUMAN"]
    sample_ids = ["TPAD1001"]
    metadata_rows = [{"SampleID": "TPAD1001", "Condition": "ADD", "Run Order": "1"}]
    _write_metadata_csv(tmp_path / "metadata.csv", metadata_rows)
    matrix = np.array([[20.0], [21.0]])
    _write_data_csv(tmp_path / "protein-data.csv", feature_ids, sample_ids, matrix)

    dataset = load_paired_dataset(
        tmp_path / "protein-data.csv",
        tmp_path / "metadata.csv",
        prepared_dir=tmp_path / "prepared",
        expected_n_samples=1,
        expected_n_features=2,
    )
    assert dataset.abundances.shape == (1, 2)
    assert dataset.metadata.index.tolist() == ["TPAD1001"]


def test_run_order_ties_break_stably(tmp_path: Path) -> None:
    """Two samples with an identical Run Order keep their original file order."""
    feature_ids = ["sp|P00001|TEST_HUMAN"]
    sample_ids = ["TPAD1001", "TPAD1002", "TPAD1003"]
    metadata_rows = [
        {"SampleID": "TPAD1001", "Condition": "ADD", "Run Order": "5"},
        {"SampleID": "TPAD1002", "Condition": "HCN", "Run Order": "5"},  # tie
        {"SampleID": "TPAD1003", "Condition": "PDD", "Run Order": "1"},
    ]
    _write_metadata_csv(tmp_path / "metadata.csv", metadata_rows)
    matrix = np.array([[20.0, 21.0, 22.0]])
    _write_data_csv(tmp_path / "protein-data.csv", feature_ids, sample_ids, matrix)

    dataset = load_paired_dataset(
        tmp_path / "protein-data.csv",
        tmp_path / "metadata.csv",
        prepared_dir=tmp_path / "prepared",
        expected_n_samples=3,
        expected_n_features=1,
    )
    # TPAD1003 (Run Order 1) first; then the tie at 5 keeps metadata-file order:
    # TPAD1001 before TPAD1002.
    assert dataset.metadata.index.tolist() == ["TPAD1003", "TPAD1001", "TPAD1002"]


# --------------------------------------------------------------------------
# read_full_metadata / classify_conditions / reconcile_sample_partition — units
# --------------------------------------------------------------------------


def test_classify_conditions_counts(tmp_path: Path) -> None:
    metadata_rows = [
        {"SampleID": "TPAD1001", "Condition": "ADD", "Run Order": "1"},
        {"SampleID": "TPAD1002", "Condition": "PDD", "Run Order": "2"},
        {"SampleID": "CP01", "Condition": "Plasma Reference", "Run Order": ""},
    ]
    metadata = pd.DataFrame(metadata_rows)
    experimental, control = classify_conditions(metadata)
    assert len(experimental) == 2
    assert len(control) == 1
    assert control["SampleID"].tolist() == ["CP01"]


def test_reconcile_sample_partition_control_with_data_column_fails_loud() -> None:
    """A control SampleID that leaked a data column is caught as "not an
    experimental SampleID" (a control id is, by construction, never one)."""
    experimental = pd.DataFrame({"SampleID": ["TPAD1001"]})
    control = pd.DataFrame({"SampleID": ["CP01"]})
    with pytest.raises(ValueError, match="do not match any experimental SampleID"):
        reconcile_sample_partition(
            experimental, control, ["TPAD1001", "CP01"], data_file=Path("x.csv")
        )


def test_read_full_metadata_missing_join_key_fails_loud(tmp_path: Path) -> None:
    path = tmp_path / "metadata.csv"
    path.write_text("NotSampleID,Condition\nA,ADD\n")
    with pytest.raises(ValueError, match="missing required column"):
        read_full_metadata(path)


def test_read_full_metadata_na_token_sample_id_fails_loud(tmp_path: Path) -> None:
    """A SampleID that spells a pandas NA token (e.g. "None") must fail loud, not
    silently become NaN and surface later as a confusing orphan."""
    path = tmp_path / "metadata.csv"
    path.write_text("SampleID,Condition\nNone,ADD\nTPAD1002,HCN\n")
    with pytest.raises(ValueError, match="Missing/NA-token"):
        read_full_metadata(path)


def test_read_data_sample_columns_requires_empty_first_header(tmp_path: Path) -> None:
    path = tmp_path / "data.csv"
    path.write_text("protein,A,B\nfeat1,1.0,2.0\n")
    with pytest.raises(ValueError, match="empty"):
        read_data_sample_columns(path)


# --------------------------------------------------------------------------
# Obligation B: reconciliation on the REAL project files
# --------------------------------------------------------------------------


@pytest.mark.skipif(
    not (REAL_DATA_FILE.exists() and REAL_METADATA_FILE.exists()),
    reason="real project data files not present",
)
def test_real_data_reconciliation(tmp_path: Path) -> None:
    full_metadata = read_full_metadata(REAL_METADATA_FILE)
    assert len(full_metadata) == 48

    experimental, control = classify_conditions(full_metadata)
    assert len(experimental) == 40
    assert len(control) == 8
    for condition in ("ADD", "PDD", "PDCN", "HCN"):
        assert int((experimental["Condition"] == condition).sum()) == 10
    assert set(control["SampleID"]) == {
        "CP01",
        "CP02",
        "CP03",
        "CP04",
        "KP01",
        "KP02",
        "KP03",
        "KP04",
    }

    dataset = load_paired_dataset(
        REAL_DATA_FILE,
        REAL_METADATA_FILE,
        prepared_dir=tmp_path / "prepared",
        expected_n_control=8,
    )
    assert dataset.abundances.shape == (40, 2343)
    assert dataset.scale == "log2"
    assert len(set(dataset.feature_names.tolist())) == 2343

    # Spot-check a raw value against the source file (row 1, sample TPAD1010).
    # (pandas renames the file's empty first header to "Unnamed: 0"; this is only
    # this test's independent re-read of the source for reconciliation, not how
    # project_loader itself parses the id column.)
    raw = pd.read_csv(REAL_DATA_FILE, dtype={"Unnamed: 0": str})
    raw_value = float(
        raw.loc[raw["Unnamed: 0"] == "sp|O15194|CTDSL_HUMAN", "TPAD1010"].iloc[0]
    )
    feat_idx = list(dataset.feature_names).index("sp|O15194|CTDSL_HUMAN")
    row_idx = dataset.metadata.index.tolist().index("TPAD1010")
    assert dataset.abundances[row_idx, feat_idx] == pytest.approx(raw_value)

    trimmed = load_analysis_subset(dataset)
    assert trimmed.abundances.shape == (40, 2343 - len(ANALYSIS_DROP_FEATURE_IDS))
    assert trimmed.scale == "log2"
    remaining = set(trimmed.feature_names.tolist())
    for dropped_id in ANALYSIS_DROP_FEATURE_IDS:
        assert dropped_id not in remaining
    # Samples untouched by the analysis-subset trim.
    assert trimmed.metadata.index.tolist() == dataset.metadata.index.tolist()
