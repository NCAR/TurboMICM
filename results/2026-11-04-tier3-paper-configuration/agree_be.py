#!/usr/bin/env python3
"""Agreement of the library's backward Euler with the Python prototype: digits at
MICM's default tolerances on the 64 cells of all eleven scenarios, and with
equal sub-steps on the three scenarios of the paper's comparison."""
import csv, json, re, subprocess, sys, tempfile
from pathlib import Path

HOME = Path.home()
BIN = HOME / "NSF-NCAR/TurboMICM-paper/bin/ts1_scenario"
ARCHIVE = HOME / "DavidFillmore/Fast-ODE-Atmospheric-Chem/data/scenarios"
RUNS = HOME / "DavidFillmore/ODE/fastode/results/ts1/runs.csv"
SUBSTEPS = HOME / "DavidFillmore/ODE/fastode/results/ts1/be_substeps.csv"
SCORE = HOME / "NSF-NCAR/TurboMICM/tier3/scenarios/score.py"
SCENARIOS = ["tier1-30s", "diurnal-30min", "diurnal-10min", "diurnal-5min", "polluted-10min", "clean-10min",
             "summer-15min", "tropics-20min", "winter60-10min", "polluted-48h-30min", "clean-48h-5min"]

def run(scenario, substeps):
    scn = ARCHIVE / scenario / "scenario.txt"
    with tempfile.TemporaryDirectory() as d:
        args = [str(BIN), "--threads", "1", "--method", "be", "--option", "rtol=1e-6", "--option", "atol=1e-3",
                "--option", f"substeps={substeps}", str(scn)]
        out = subprocess.run(args, cwd=d, capture_output=True, text=True)
        if out.returncode != 0:
            return None, out.stderr.strip()[:200]
        stats = dict(re.findall(r"(\w+)=([-\d.]+)", out.stdout))
        snaps = (Path(d) / "ts1_snapshots.csv").exists()
        res = Path(d) / ("ts1_snapshots.csv" if snaps else "ts1_final.csv")
        ref = ARCHIVE / scenario / ("reference_snapshots.csv" if snaps else "reference_final.csv")
        sc = subprocess.run([sys.executable, str(SCORE), str(res), str(ref), "--json"], capture_output=True, text=True)
        j = json.loads(sc.stdout)
        return (float(j["overall_digits"]), stats), None

py = {r["scenario"]: r for r in csv.DictReader(open(RUNS)) if r["config"] == "be" and r["cells"] == "64"}
print("== backward Euler at MICM's defaults (rtol 1e-6, atol 1e-3), one step per interval")
print(f"{'scenario':20s} {'C digits':>9s} {'py':>8s}  C stats")
bad = 0
for s in SCENARIOS:
    r, err = run(s, 1)
    if r is None:
        print(s, "FAILED", err); bad += 1; continue
    d, st = r
    flag = "same" if abs(d - round(float(py[s]["digits"]), 3)) <= 0.0015 else "DIFFER"
    bad += flag != "same"
    print(f"{s:20s} {d:9.3f} {float(py[s]['digits']):8.4f}  steps={st.get('internal_steps')} jac={st.get('jacobian_updates')} lu={st.get('decompositions')} solves={st.get('solves')} rhs={st.get('function_calls')} {flag}")
print("== equal sub-steps (the paper's Table of the backward Euler comparison)")
subs = {(r["scenario"], int(r["substeps"])): r for r in csv.DictReader(open(SUBSTEPS))}
for (s, n), r in sorted(subs.items(), key=lambda kv: (SCENARIOS.index(kv[0][0]), kv[0][1])):
    if n == 1:
        continue
    res, err = run(s, n)
    if res is None:
        print(s, n, "FAILED", err); bad += 1; continue
    d, st = res
    flag = "same" if abs(d - round(float(r["digits"]), 3)) <= 0.0015 else "DIFFER"
    bad += flag != "same"
    print(f"{s:20s} x{n:<4d} {d:9.3f} {float(r['digits']):8.4f}  newton/interval py={float(r['newton_per_interval']):.0f} C jac={st.get('jacobian_updates')} {flag}")
print("mismatches:", bad)
