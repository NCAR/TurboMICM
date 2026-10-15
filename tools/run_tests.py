"""Print every reference case exactly once; absent implementations fail closed."""
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parents[1]
cases = [line.split()[0] for line in (root / "tests/cases.txt").read_text().splitlines()
         if line and not line.startswith("#")]
result = subprocess.run([root / "bin/test_cpu"], text=True, capture_output=True)
sys.stderr.write(result.stderr)
observed = {}
for line in result.stdout.splitlines():
    key, status, reason = line.split(" ", 2)
    if key in observed or key not in cases or status not in {"PASS", "FAIL", "N/A"}:
        raise RuntimeError(f"Invalid test report: {line}")
    observed[key] = (status, reason)
counts = {"PASS": 0, "FAIL": 0, "N/A": 0}
for key in cases:
    status, reason = observed.get(key, ("FAIL", "case port pending"))
    counts[status] += 1
    print(key, status, reason)
print(" ".join(f"{k}={v}" for k, v in counts.items()), file=sys.stderr)
sys.exit(bool(result.returncode or counts["FAIL"]))
