"""Regression checks for stale data and incorrectly marked passing evidence."""
import pytest
from evidence_utils import hashes, check_hashes, check_validation


def test_reject_changed_output_even_when_old_success_record_exists(tmp_path):
    (tmp_path / "errors.csv").write_text("0.00125\n")
    recorded = hashes(tmp_path, ("errors.csv",))
    check_hashes(tmp_path, recorded, ("errors.csv",))
    (tmp_path / "errors.csv").write_text("0.0\n")
    with pytest.raises(ValueError, match="stale"):
        check_hashes(tmp_path, recorded, ("errors.csv",))


def test_reject_incomplete_manifest(tmp_path):
    with pytest.raises(ValueError, match="missing"):
        check_hashes(tmp_path, {}, ("test_experiments.py",))


@pytest.mark.parametrize("value", [1.0, float("nan"), float("inf")])
def test_recompute_values_instead_of_trusting_pass_flag(value):
    record = {"passed": True, "number_of_checks": 1,
              "checks": [{"name": "residual", "value": value, "limit": 1e-10,
                           "comparison": "<=", "passed": True}]}
    with pytest.raises(ValueError):
        check_validation(record)
