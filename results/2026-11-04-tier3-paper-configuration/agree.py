#!/usr/bin/env python3
"""Agreement test: the C library's retaining solver and restarting loop against
the archived Python solver's operation counts and digits on the 64 cells of the
eleven scenarios (two groups of 32, one thread)."""
import csv, os, re, subprocess, sys, tempfile
from pathlib import Path

HOME = Path.home()
BIN = HOME / "NSF-NCAR/TurboMICM-paper/bin/ts1_scenario"
ARCHIVE = HOME / "DavidFillmore/Fast-ODE-Atmospheric-Chem/data/scenarios"
RUNS = HOME / "DavidFillmore/ODE/fastode/results/ts1/runs.csv"
SCORE = HOME / "NSF-NCAR/TurboMICM/tier3/scenarios/score.py"
SCENARIOS = ["tier1-30s", "diurnal-30min", "diurnal-10min", "diurnal-5min", "polluted-10min", "clean-10min",
             "summer-15min", "tropics-20min", "winter60-10min", "polluted-48h-30min", "clean-48h-5min"]
PAPER = ["--option", "quantize=1", "--option", "controller=h211b", "--option", "safety=0.9",
         "--option", "exact_end=1", "--option", "reject_elementary=1", "--option", "strict_min=1"]
LOOP = ["--option", "warm_start=0", "--option", "quantize=0", "--option", "controller=micm"]
CONFIGS = {
    "retaining rtol 1e-2": ("alg", 0.01, ["--method", "ros3", "--option", "rtol=0.01", "--option", "atol=1e-16"] + PAPER),
    "retaining rtol 5e-3": ("alg", 0.005, ["--method", "ros3", "--option", "rtol=0.005", "--option", "atol=5e-17"] + PAPER),
    "restarting rtol 3e-2": ("restart", 0.03, ["--method", "ros3", "--option", "rtol=0.03", "--option", "atol=1e-15"] + LOOP),
    "restarting rtol 3e-4": ("restart", 0.0003, ["--method", "ros3", "--option", "rtol=0.0003", "--option", "atol=1e-15"] + LOOP),
}
# jacobian_updates is compared separately: the library re-evaluates the Jacobian after a rejection
# (its in-place factorization overwrites it), so its count exceeds the Python count by the rejections.
KEYS = [("internal_steps", "steps"), ("accepted", "accepted"), ("rejected", "rejected"), ("function_calls", "rhs"),
        ("decompositions", "lu"), ("solves", "solves")]

targets = {}
for r in csv.DictReader(open(RUNS)):
    if r["cells"] == "64":
        targets[(r["config"], float(r["rtol"]), r["scenario"])] = r

only = sys.argv[1:]  # optional config-name filters
mismatches = 0
for name, (config, rtol, args) in CONFIGS.items():
    if only and not any(o in name for o in only):
        continue
    print(f"== {name}")
    print(f"{'scenario':20s} {'C steps':>8s} {'py':>6s} {'rej':>5s} {'py':>5s} {'rhs':>6s} {'py':>6s} {'solves':>7s} {'py':>6s} {'C digits':>9s} {'py':>7s} counts")
    for s in SCENARIOS:
        scn = ARCHIVE / s / "scenario.txt"
        t = targets[(config, rtol, s)]
        with tempfile.TemporaryDirectory() as d:
            out = subprocess.run([str(BIN), "--threads", "1"] + args + [str(scn)], cwd=d, capture_output=True, text=True)
            if out.returncode != 0:
                print(s, "FAILED:", out.stderr.strip()[:200]); mismatches += 1; continue
            stats = dict(re.findall(r"(\w+)=([-\d.]+)", out.stdout))
            ref = ARCHIVE / s / ("reference_snapshots.csv" if (Path(d) / "ts1_snapshots.csv").exists() else "reference_final.csv")
            res = Path(d) / ("ts1_snapshots.csv" if (Path(d) / "ts1_snapshots.csv").exists() else "ts1_final.csv")
            sc = subprocess.run([sys.executable, str(SCORE), str(res), str(ref), "--json"], capture_output=True, text=True)
            try:
                import json; j = json.loads(sc.stdout); digits = float(j.get("overall_digits", float("nan")))
            except Exception:
                m = re.search(r"overall[^0-9-]*(-?\d+\.\d+)", sc.stdout); digits = float(m.group(1)) if m else float("nan")
        same = all(int(float(stats[c])) == int(t[p]) for c, p in KEYS) and int(float(stats['jacobian_updates'])) == int(t['jac']) + int(t['rejected'])
        mismatches += 0 if same else 1
        print(f"{s:20s} {int(float(stats['internal_steps'])):8d} {int(t['steps']):6d} {int(float(stats['rejected'])):5d} {int(t['rejected']):5d} "
              f"{int(float(stats['function_calls'])):6d} {int(t['rhs']):6d} {int(float(stats['solves'])):7d} {int(t['solves']):6d} "
              f"{digits:9.4f} {float(t['digits']):7.4f} {'same' if same else 'DIFFER'}")
print("count mismatches:", mismatches)
