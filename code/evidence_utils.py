"""Small provenance checks shared by report generation and delivery auditing."""
import hashlib
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path

RUN_SOURCES = ("code/run_experiments.py", "code/reference_errors.py", "code/evidence_utils.py")
RUN_OUTPUTS = ("errors.csv", "diagnostics.csv", "validation.json", "nodal_solutions.npz",
               "figures/convergence.png", "figures/convergence.pdf",
               "figures/two_element_check.png", "figures/two_element_check.pdf")
TEST_SOURCES = (*RUN_SOURCES, "code/test_experiments.py", "code/test_evidence.py", "code/run_tests.py")
BUILD_SOURCES = ("code/build_reports.py", "code/chinese_template.md", "code/evidence_utils.py")
BUILD_INPUTS = ("runs/primary/metadata.json", "evidence/tests.json")
BUILD_OUTPUTS = ("first_weekly_progress_Jiahui_Gu.docx", "中文研究说明.md")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def hashes(base, names):
    return {name: digest(Path(base) / name) for name in names}


def check_hashes(base, recorded, required):
    if set(recorded) != set(required):
        raise ValueError("Evidence manifest has missing or unexpected file entries")
    for name in required:
        path = Path(base) / name
        if not path.is_file() or recorded[name] != digest(path):
            raise ValueError(f"Evidence is missing or stale: {path}")


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def check_validation(record):
    checks = record["checks"]
    if not checks or record["number_of_checks"] != len(checks) or record["passed"] is not True:
        raise ValueError("Numerical validation is empty, inconsistent or failed")
    for check in checks:
        value, limit = check["value"], check["limit"]
        comparison = check["comparison"]
        if comparison not in ("<=", ">=") or not math.isfinite(value) or not math.isfinite(limit):
            raise ValueError("Invalid numerical validation value or comparison")
        passed = value <= limit if comparison == "<=" else value >= limit
        if not passed or check["passed"] is not True:
            raise ValueError(f"Numerical check failed: {check['name']}")


def verify_run(run, root):
    meta = read_json(Path(run) / "metadata.json")
    check_hashes(root, meta["code_sha256"], RUN_SOURCES)
    check_hashes(run, meta["artifact_sha256"], RUN_OUTPUTS)
    validation = read_json(Path(run) / "validation.json")
    check_validation(validation)
    if meta["validation_passed"] is not True:
        raise ValueError("Run metadata does not report successful validation")
    return meta, validation


def junit_counts(path):
    tree = ET.parse(path)
    cases = tree.findall(".//testcase")
    return {"total": len(cases),
            "failures": sum(c.find("failure") is not None for c in cases),
            "errors": sum(c.find("error") is not None for c in cases),
            "skipped": sum(c.find("skipped") is not None for c in cases)}


def verify_tests(root):
    record = read_json(Path(root) / "evidence/tests.json")
    check_hashes(root, record["code_sha256"], TEST_SOURCES)
    check_hashes(root, record["evidence_sha256"], ("evidence/pytest.xml", "evidence/pytest.txt"))
    counts = junit_counts(Path(root) / "evidence/pytest.xml")
    if (record["counts"] != counts or record["exit_code"] != 0 or record["passed"] is not True
            or counts["total"] == 0 or any(counts[k] for k in ("failures", "errors", "skipped"))):
        raise ValueError("Current tests did not all pass")
    return record


def verify_build(root):
    record = read_json(Path(root) / "evidence/report_build.json")
    check_hashes(root, record["code_sha256"], BUILD_SOURCES)
    check_hashes(root, record["input_sha256"], BUILD_INPUTS)
    check_hashes(root, record["output_sha256"], BUILD_OUTPUTS)
    return record
