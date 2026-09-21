"""Run pytest and save an exit code, JUnit result and tested-file hashes."""
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from evidence_utils import TEST_SOURCES, hashes, junit_counts

ROOT = Path(__file__).resolve().parents[1]


def main():
    evidence = ROOT / "evidence"
    evidence.mkdir(exist_ok=True)
    before = hashes(ROOT, TEST_SOURCES)
    command = [sys.executable, "-B", "-X", "utf8", "-m", "pytest",
               str(ROOT / "code/test_experiments.py"), str(ROOT / "code/test_evidence.py"),
               "-q", "-p", "no:cacheprovider", "--junitxml", str(evidence / "pytest.xml")]
    temp_root = evidence / "test_temp"
    temp_root.mkdir(exist_ok=True)
    environment = os.environ.copy()
    environment["PYTEST_DEBUG_TEMPROOT"] = str(temp_root)
    completed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True,
                               encoding="utf-8", env=environment)
    log = completed.stdout + completed.stderr
    (evidence / "pytest.txt").write_text(log, encoding="utf-8")
    counts = junit_counts(evidence / "pytest.xml") if (evidence / "pytest.xml").exists() else {}
    unchanged = before == hashes(ROOT, TEST_SOURCES)
    passed = (completed.returncode == 0 and unchanged and counts.get("total", 0) > 0
              and all(counts.get(k) == 0 for k in ("failures", "errors", "skipped")))
    record = {"recorded_at_utc": datetime.now(timezone.utc).isoformat(), "command": command,
              "python_executable": sys.executable, "exit_code": completed.returncode,
              "counts": counts, "code_sha256": before, "code_unchanged_during_tests": unchanged,
              "evidence_sha256": hashes(ROOT, ("evidence/pytest.xml", "evidence/pytest.txt")),
              "passed": passed}
    (evidence / "tests.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(log, end="")
    print(f"Test evidence current: {passed}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
