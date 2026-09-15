"""Tests for ``scripts/scratch/qc_prep.py`` (the Stage-3 QC prep-once script).

Layers:
  * unit -- hand-built fixtures for ``build_delivered_states`` (the single-state
    log2/linear derivation) and ``assess_confounding`` (the warning-capturing wrapper),
    including a planted-truth confounding check and edge cases (fail-loud guards,
    single sample, near-degenerate crosstabs);
  * serialization -- ``report_to_jsonable`` / ``write_confounding_report`` /
    ``write_readme`` produce the artifacts the README/report contract promises, and
    round-trip through ``json.loads``;
  * integration on the REAL project files -- ``run()`` end-to-end into a temp output
    directory, verifying shapes/scales and that the confounding report corroborates
    findings 0001/0004 without creating a new finding.

Run:
    ./.venv/bin/python -m pytest scripts/scratch/test_qc_prep.py -v
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# `common.*` lives at scripts/common/; this test imports it directly (to build
# fixtures / assert on its exception types), so it needs the same sys.path setup
# qc_prep.py itself performs -- pytest's own rootdir insertion only adds this file's
# directory (scripts/scratch, no __init__.py there), not scripts/.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.batch_correct import BatchConfoundingWarning
from common.data_loading import Dataset
from common.dataset_io import load_dataset
from qc_prep import (
    CORROBORATES_FINDINGS,
    DeliveredStates,
    _parse_args,
    assess_confounding,
    build_delivered_states,
    report_to_jsonable,
    run,
    sha256_16,
    write_confounding_report,
    write_readme,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
REAL_DATA_FILE = PROJECT_ROOT / "data" / "protein-data.csv"
REAL_METADATA_FILE = PROJECT_ROOT / "data" / "metadata.csv"

# --------------------------------------------------------------------------
# Fixture builders
# --------------------------------------------------------------------------


def _ds(
    abundances: np.ndarray,
    *,
    scale: str = "log2",
    batch: list[str] | None = None,
    covariate: list[str] | None = None,
    extra_meta: dict[str, list[str]] | None = None,
) -> Dataset:
    """A small hand-built Dataset for unit tests."""
    ab = np.asarray(abundances, dtype=float)
    n_samples, n_features = ab.shape
    feature_names = np.array([f"F{j}" for j in range(n_features)], dtype=str)
    feature_metadata = pd.DataFrame({"protein": feature_names})
    data: dict[str, list[str]] = {"sample": [f"s{i}" for i in range(n_samples)]}
    if batch is not None:
        data["batch"] = list(batch)
    if covariate is not None:
        data["covariate"] = list(covariate)
    if extra_meta is not None:
        for key, values in extra_meta.items():
            data[key] = list(values)
    metadata = pd.DataFrame(data).set_index("sample")
    return Dataset(
        abundances=ab,
        feature_names=feature_names,
        feature_metadata=feature_metadata,
        metadata=metadata,
        scale=scale,
    )


# --------------------------------------------------------------------------
# build_delivered_states -- unit + planted-truth + edge cases
# --------------------------------------------------------------------------


def test_build_delivered_states_hand_verified() -> None:
    """A tiny 2x2 log2 matrix; the linear image must equal 2**log2 exactly."""
    log_values = np.array([[1.0, 2.0], [3.0, 0.0]])
    ds = _ds(log_values, scale="log2")
    states = build_delivered_states(ds)

    assert isinstance(states, DeliveredStates)
    assert states.log is ds
    assert states.log.scale == "log2"
    assert states.linear.scale == "linear"
    np.testing.assert_allclose(states.linear.abundances, np.exp2(log_values))
    # Known values: 2**1=2, 2**2=4, 2**3=8, 2**0=1.
    np.testing.assert_allclose(
        states.linear.abundances, np.array([[2.0, 4.0], [8.0, 1.0]])
    )


def test_build_delivered_states_shape_and_feature_names_preserved() -> None:
    rng = np.random.default_rng(0)
    ds = _ds(rng.uniform(10, 30, size=(5, 7)), scale="log2")
    states = build_delivered_states(ds)
    assert states.linear.abundances.shape == ds.abundances.shape
    assert list(states.linear.feature_names) == list(ds.feature_names)
    pd.testing.assert_frame_equal(states.linear.metadata, ds.metadata)


def test_build_delivered_states_linear_is_strictly_positive() -> None:
    """Planted-truth-adjacent: every de-logged value must recover positivity even for
    large-magnitude and near-zero log2 inputs (no accidental zero/negative)."""
    ds = _ds(np.array([[-5.0, 0.0, 30.0, 38.7]]), scale="log2")
    states = build_delivered_states(ds)
    assert bool(np.all(states.linear.abundances > 0))
    np.testing.assert_allclose(
        states.linear.abundances, np.exp2(np.array([[-5.0, 0.0, 30.0, 38.7]]))
    )


def test_build_delivered_states_rejects_non_log2_scale() -> None:
    ds = _ds(np.array([[1.0, 2.0]]), scale="linear")
    with pytest.raises(ValueError, match="scale 'log2'"):
        build_delivered_states(ds)


def test_build_delivered_states_rejects_glog2_scale() -> None:
    """glog2 is a real Scale but has no clean linear inverse -- must fail loud too."""
    ds = _ds(np.array([[1.0, 2.0]]), scale="glog2")
    with pytest.raises(ValueError, match="scale 'log2'"):
        build_delivered_states(ds)


def test_build_delivered_states_rejects_non_finite() -> None:
    ds = _ds(np.array([[1.0, np.nan]]), scale="log2")
    with pytest.raises(ValueError, match="non-finite"):
        build_delivered_states(ds)


def test_build_delivered_states_single_sample_single_feature() -> None:
    """Edge case: the smallest possible Dataset (1x1) still round-trips correctly."""
    ds = _ds(np.array([[5.0]]), scale="log2")
    states = build_delivered_states(ds)
    assert states.linear.abundances.shape == (1, 1)
    assert states.linear.abundances[0, 0] == pytest.approx(32.0)  # 2**5


# --------------------------------------------------------------------------
# assess_confounding -- unit, planted confound, and captured-warning behavior
# --------------------------------------------------------------------------


def test_assess_confounding_perfect_confound_captures_warning() -> None:
    """Planted-truth: a batch/covariate pair nested 1:1 must be flagged PERFECT and
    the emitted BatchConfoundingWarning must be captured, not raised to the caller."""
    ds = _ds(
        np.zeros((4, 3)),
        scale="log2",
        batch=["A", "A", "B", "B"],
        covariate=["x", "x", "y", "y"],
    )
    with warnings.catch_warnings():
        # Any escaping warning would fail the test under -W error; assert none escapes.
        warnings.simplefilter("error", BatchConfoundingWarning)
        reports, captured = assess_confounding(
            ds,
            batch_column="batch",
            covariate_columns=["covariate"],
            warn_threshold=0.5,
        )
    assert len(reports) == 1
    assert reports[0].perfectly_confounded
    assert reports[0].cramers_v == pytest.approx(1.0)
    assert len(captured) == 1
    assert captured[0]["category"] == "BatchConfoundingWarning"
    assert "PERFECTLY confounded" in captured[0]["message"]


def test_assess_confounding_independent_no_warning_captured() -> None:
    ds = _ds(
        np.zeros((4, 3)),
        scale="log2",
        batch=["A", "A", "B", "B"],
        covariate=["x", "y", "x", "y"],
    )
    reports, captured = assess_confounding(
        ds, batch_column="batch", covariate_columns=["covariate"], warn_threshold=0.5
    )
    assert not reports[0].perfectly_confounded
    assert reports[0].cramers_v == pytest.approx(0.0)
    assert captured == []


def test_assess_confounding_multiple_covariates_preserve_order() -> None:
    ds = _ds(
        np.zeros((4, 3)),
        scale="log2",
        batch=["A", "A", "B", "B"],
        covariate=["x", "y", "x", "y"],
        extra_meta={"other": ["p", "p", "q", "q"]},
    )
    reports, _ = assess_confounding(
        ds,
        batch_column="batch",
        covariate_columns=["covariate", "other"],
        warn_threshold=0.5,
    )
    assert [r.covariate_column for r in reports] == ["covariate", "other"]
    assert not reports[0].perfectly_confounded
    assert reports[1].perfectly_confounded


def test_assess_confounding_propagates_missing_column_error() -> None:
    """A programming/config error (unknown column) must NOT be swallowed as a
    warning."""
    ds = _ds(np.zeros((4, 3)), scale="log2", batch=["A", "A", "B", "B"])
    with pytest.raises(ValueError, match="not a metadata column"):
        assess_confounding(
            ds, batch_column="batch", covariate_columns=["nope"], warn_threshold=0.5
        )


# --------------------------------------------------------------------------
# report_to_jsonable / write_confounding_report -- serialization
# --------------------------------------------------------------------------


def test_report_to_jsonable_roundtrips_through_json() -> None:
    ds = _ds(
        np.zeros((4, 3)),
        scale="log2",
        batch=["A", "A", "B", "B"],
        covariate=["x", "x", "y", "y"],
    )
    reports, _ = assess_confounding(
        ds, batch_column="batch", covariate_columns=["covariate"], warn_threshold=0.5
    )
    payload = report_to_jsonable(reports[0])
    # Must be JSON-serializable (no numpy/pandas types leaking through).
    text = json.dumps(payload)
    reloaded = json.loads(text)
    assert reloaded["batch_column"] == "batch"
    assert reloaded["covariate_column"] == "covariate"
    assert reloaded["perfectly_confounded"] is True
    assert reloaded["crosstab"]["A"]["x"] == 2
    assert reloaded["crosstab"]["B"]["y"] == 2
    assert reloaded["crosstab"]["A"]["y"] == 0


def test_write_confounding_report_writes_valid_json_with_corroboration(
    tmp_path: Path,
) -> None:
    ds = _ds(
        np.zeros((4, 3)),
        scale="log2",
        batch=["A", "A", "B", "B"],
        covariate=["x", "x", "y", "y"],
    )
    reports, captured = assess_confounding(
        ds, batch_column="batch", covariate_columns=["covariate"], warn_threshold=0.5
    )
    out = tmp_path / "report.json"
    path = write_confounding_report(
        reports,
        captured,
        output_path=out,
        data_version={"data": "deadbeef", "metadata": "cafef00d"},
    )
    assert path == out
    payload = json.loads(out.read_text())
    assert payload["data_version"] == {"data": "deadbeef", "metadata": "cafef00d"}
    assert payload["corroborates_findings"] == list(CORROBORATES_FINDINGS)
    assert "0001" in payload["corroborates_findings"]
    assert "0004" in payload["corroborates_findings"]
    assert "not itself a new finding" in payload["note"]
    assert len(payload["reports"]) == 1
    assert payload["warnings_captured"][0]["category"] == "BatchConfoundingWarning"


def test_write_confounding_report_empty_reports_edge_case(tmp_path: Path) -> None:
    """Edge case: an empty report/warning list must still write valid JSON."""
    out = tmp_path / "sub" / "report.json"
    path = write_confounding_report(
        [], [], output_path=out, data_version={"data": "x", "metadata": "y"}
    )
    payload = json.loads(path.read_text())
    assert payload["reports"] == []
    assert payload["warnings_captured"] == []


# --------------------------------------------------------------------------
# write_readme -- content contract
# --------------------------------------------------------------------------


def test_write_readme_states_single_delivered_state(tmp_path: Path) -> None:
    path = write_readme(tmp_path)
    assert path == tmp_path / "README.md"
    text = path.read_text()
    assert "only one delivered processing state" in text
    assert "0004" in text
    assert "0003" in text
    assert "NOT raw" in text or "NOT raw / unnormalized intensity" in text
    assert "delivered_log" in text
    assert "delivered_linear" in text
    assert "pooled-reference" in text.lower() or "pooled reference" in text.lower()


# --------------------------------------------------------------------------
# sha256_16
# --------------------------------------------------------------------------


def test_sha256_16_matches_known_hash(tmp_path: Path) -> None:
    f = tmp_path / "x.txt"
    f.write_text("hello world")
    import hashlib

    expected = hashlib.sha256(b"hello world").hexdigest()[:16]
    assert sha256_16(f) == expected
    assert len(sha256_16(f)) == 16


# --------------------------------------------------------------------------
# Argument parsing
# --------------------------------------------------------------------------


def test_parse_args_defaults() -> None:
    args = _parse_args([])
    assert args.batch_column == "Cohort"
    assert args.covariate_columns == ["Condition"]
    assert args.warn_threshold == pytest.approx(0.5)
    assert args.expected_n_samples == 40


def test_parse_args_overrides() -> None:
    args = _parse_args(
        [
            "--batch-column",
            "SiteX",
            "--covariate-columns",
            "Condition",
            "Sex",
            "--warn-threshold",
            "0.7",
        ]
    )
    assert args.batch_column == "SiteX"
    assert args.covariate_columns == ["Condition", "Sex"]
    assert args.warn_threshold == pytest.approx(0.7)


# --------------------------------------------------------------------------
# Integration: run() end-to-end on the REAL project files
# --------------------------------------------------------------------------


@pytest.mark.skipif(
    not (REAL_DATA_FILE.exists() and REAL_METADATA_FILE.exists()),
    reason="real project data files not present",
)
def test_run_end_to_end_on_real_data(tmp_path: Path) -> None:
    output_dir = tmp_path / "qc_states"
    args = argparse.Namespace(
        data=REAL_DATA_FILE,
        metadata=REAL_METADATA_FILE,
        prepared_dir=tmp_path / "prepared",
        output_dir=output_dir,
        batch_column="Cohort",
        covariate_columns=["Condition"],
        warn_threshold=0.5,
        expected_n_samples=40,
        log_level="INFO",
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BatchConfoundingWarning)
        run(args)

    # Both states were written and are independently loadable.
    delivered_log = load_dataset(output_dir / "delivered_log")
    delivered_linear = load_dataset(output_dir / "delivered_linear")

    assert delivered_log.abundances.shape == (40, 2340)
    assert delivered_log.scale == "log2"
    assert delivered_linear.abundances.shape == (40, 2340)
    assert delivered_linear.scale == "linear"
    assert bool(np.all(delivered_linear.abundances > 0))
    # The linear state is the exact de-logged image of the log state.
    np.testing.assert_allclose(
        delivered_linear.abundances, np.exp2(delivered_log.abundances)
    )
    # Same feature order/names and same sample pairing in both states.
    assert list(delivered_log.feature_names) == list(delivered_linear.feature_names)
    assert (
        delivered_log.metadata.index.tolist()
        == delivered_linear.metadata.index.tolist()
    )

    # The confounding report was written and corroborates 0001/0004 (Cohort/Condition
    # is known -- Cramer's V = 0.515 group x cohort per finding 0001/results -- so the
    # real-data run must find nonzero, warned confounding, not independence).
    report_path = output_dir / "batch_confounding_report.json"
    assert report_path.is_file()
    payload = json.loads(report_path.read_text())
    assert payload["reports"][0]["batch_column"] == "Cohort"
    assert payload["reports"][0]["covariate_column"] == "Condition"
    assert payload["reports"][0]["cramers_v"] > 0.0
    assert payload["corroborates_findings"] == ["0001", "0004"]

    readme_path = output_dir / "README.md"
    assert readme_path.is_file()
    assert "only one delivered processing state" in readme_path.read_text()
