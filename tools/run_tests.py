"""Print every reference case exactly once; absent implementations fail closed."""
from pathlib import Path
import subprocess
import sys
from test_tutorials import run_tutorials

root = Path(__file__).resolve().parents[1]
cases = [line.split()[0] for line in (root / "tests/cases.txt").read_text().splitlines()
         if line and not line.startswith("#")]
extra = subprocess.run([root / "bin/test_config"], text=True, capture_output=True)
sys.stderr.write(extra.stdout + extra.stderr)
ts1 = subprocess.run([sys.executable, root / "tools/test_ts1.py"], text=True, capture_output=True)
sys.stderr.write(ts1.stdout + ts1.stderr)
result = subprocess.run([root / "bin/test_cpu"], text=True, capture_output=True)
sys.stderr.write(result.stderr)
observed = {}
for line in result.stdout.splitlines():
    key, status, reason = line.split(" ", 2)
    if key in observed or key not in cases or status not in {"PASS", "FAIL", "N/A"}:
        raise RuntimeError(f"Invalid test report: {line}")
    observed[key] = (status, reason)
for key, report in run_tutorials(root).items():
    if key in observed:
        raise RuntimeError(f"Duplicate program report: {key}")
    observed[key] = report
counts = {"PASS": 0, "FAIL": 0, "N/A": 0}
for key in cases:
    status, reason = observed.get(key, ("FAIL", "case port pending"))
    counts[status] += 1
    print(key, status, reason)
print(" ".join(f"{k}={v}" for k, v in counts.items()), file=sys.stderr)
sys.exit(bool(extra.returncode or ts1.returncode or result.returncode or counts["FAIL"]))
