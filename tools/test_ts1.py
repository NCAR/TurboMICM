"""TS1 regression for post-interval nonnegative concentration clamping.

Input is cell 3974 from the reference generator's 10,000-cell seed-2026 input.
Expected data is from MICM, built with -O3 -DNDEBUG, 100 timed steps of 30 s
plus one warm-up step. Without the public solver's interval-end clamp, this
cell diverges by 5e-4 relative although the smaller benchmark passes.
"""
import csv
import math
from pathlib import Path
import subprocess

root = Path(__file__).resolve().parents[1]
run = root / "build/test-ts1"
run.mkdir(parents=True, exist_ok=True)
subprocess.run([root / "bin/ts1_bench_cpu", "--threads", "1", "1", "100", "30",
                root / "tests/data/ts1_overshoot_input.csv"], cwd=run, check=True, capture_output=True)
with (run / "ts1_final.csv").open() as f:
    actual = next(csv.DictReader(f))
with (root / "tests/data/ts1_overshoot_expected.csv").open() as f:
    expected = next(csv.DictReader(f))
assert actual.keys() == expected.keys()
for name, value in actual.items():
    a, e = float(value), float(expected[name])
    assert math.isfinite(a) and abs(a - e) <= 1e-30 + 1e-10 * abs(e), (name, a, e)
print("TS1 interval-clamp regression: PASS")
