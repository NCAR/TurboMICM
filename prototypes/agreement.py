"""The test that each C method agrees with its prototype.

    python prototypes/agreement.py [--scenario DIR] [--methods rodas3,extrap,...] [--keep]

On a small scenario (prototypes/mini/, made by make_mini.py: the first cells
of a disclosed scenario for a few steps), runs bin/ts1_scenario and the
prototype with the same method and options, one group, one thread, and
compares: the final states (and snapshots) must agree to roundoff, and the
operation counts (internal steps, accepted, rejected, function calls,
Jacobian updates, decompositions, solves, matrix reuses) must be equal.
Prints one line per case, PASS or FAIL with the largest relative difference,
and exits 1 on any failure. Expected differences, and why, are in
prototypes/README.md.
"""
import argparse
import csv
import re
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ENTRY = HERE.parent
PYTHON = sys.executable

CASES = [
    ("micm3", []),
    ("ros2", ["rtol=1e-3", "atol=1e-18"]),
    ("ros3", ["rtol=1e-3", "atol=1e-18", "controller=predictive"]),
    ("rodas3", ["rtol=1e-3", "atol=1e-18"]),
    ("rodas3", ["rtol=1e-3", "atol=1e-18", "norm=cellmax", "clamp=1", "quantize=1"]),
    ("rodas4", ["rtol=1e-3", "atol=1e-18", "controller=h211b"]),
    ("ros2", ["rtol=1e-2", "atol=1e-18", "reuse_age=50", "reuse_change=0.2", "quantize=1", "freeze_band=1.5"]),
    ("ros34pw2", ["rtol=1e-3", "atol=1e-18"]),
    ("ros34pw2", ["rtol=1e-2", "atol=1e-18", "reuse_age=50", "reuse_change=0.2", "quantize=1", "freeze_band=1.5", "controller=h211b", "h211_b=1", "h211_k=1.7", "safety=1.3"]),
    ("extrap", ["rtol=1e-3", "atol=1e-18", "order=3"]),
    ("extrap", ["rtol=1e-3", "atol=1e-18", "order=4", "sequence=shifted"]),
]
COUNTS = ["internal_steps", "accepted", "rejected", "function_calls", "jacobian_updates", "decompositions", "solves", "matrix_reuses"]
PROTO_NAMES = {"internal_steps": "steps"}


def read_table(path):
    with open(path, newline="") as file:
        rows = list(csv.reader(file))
    return rows[0], [[float(x) for x in r] for r in rows[1:] if r]


def compare(a_path, b_path):
    """Largest relative difference between two CSVs with the same rows; species matched by name."""
    ha, ra = read_table(a_path)
    hb, rb = read_table(b_path)
    if len(ra) != len(rb):
        return float("inf")
    index = {name: i for i, name in enumerate(hb)}
    worst = 0.0
    for x, y in zip(ra, rb):
        for name, va in zip(ha, x):
            vb = y[index[name]]
            scale = max(abs(va), abs(vb), 1e-30)
            worst = max(worst, abs(va - vb) / scale)
    return worst


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scenario", type=Path, default=HERE / "mini" / "scenario.txt")
    parser.add_argument("--methods", help="comma-separated subset of the cases' methods")
    parser.add_argument("--tolerance", type=float, default=1e-9, help="largest relative difference allowed")
    parser.add_argument("--keep", action="store_true", help="keep the work directories")
    args = parser.parse_args()
    binary = ENTRY / "bin" / "ts1_scenario"
    if not binary.exists():
        sys.exit("build bin/ts1_scenario first (make bench)")
    if not args.scenario.exists():
        sys.exit(f"{args.scenario} is missing; run prototypes/make_mini.py")
    wanted = set(args.methods.split(",")) if args.methods else None
    failures = 0
    for method, options in CASES:
        if wanted and method not in wanted:
            continue
        with tempfile.TemporaryDirectory(prefix="agreement-", dir=ENTRY / "build", delete=not args.keep) as work:
            work = Path(work)
            c_dir = work / "c"
            c_dir.mkdir()
            cmd = [str(binary), "--threads", "1", "--method", method, *sum((["--option", o] for o in options), []), str(args.scenario)]
            run = subprocess.run(cmd, cwd=c_dir, capture_output=True, text=True)
            if run.returncode != 0:
                print(f"FAIL {method} {' '.join(options)}: ts1_scenario exit {run.returncode}: {run.stderr.strip()[:200]}")
                failures += 1
                continue
            c_counts = dict(re.findall(r"(\w+)=(\d+)", run.stdout))
            output = c_dir / ("ts1_snapshots.csv" if (c_dir / "ts1_snapshots.csv").exists() else "ts1_final.csv")
            p_out = work / "proto.csv"
            pcmd = [PYTHON, str(HERE / "scenario.py"), str(args.scenario), "--method", method, "--cells", "0",
                    *sum((["--option", o] for o in options), []), "--write", str(p_out)]
            prun = subprocess.run(pcmd, cwd=work, capture_output=True, text=True)
            if prun.returncode != 0:
                print(f"FAIL {method} {' '.join(options)}: prototype exit {prun.returncode}: {(prun.stderr or prun.stdout).strip()[-300:]}")
                failures += 1
                continue
            p_counts = dict(re.findall(r"(\w+)=(\S+)", prun.stdout.strip().splitlines()[-1]))
            worst = compare(output, p_out)
            mismatched = [k for k in COUNTS if c_counts.get(k) != p_counts.get(PROTO_NAMES.get(k, k))]
            ok = worst <= args.tolerance and not mismatched
            failures += 0 if ok else 1
            detail = f"largest relative difference {worst:.2e}; counts " + " ".join(
                f"{k}={c_counts.get(k)}" + ("" if k not in mismatched else f"(prototype {p_counts.get(PROTO_NAMES.get(k, k))})") for k in COUNTS)
            print(f"{'PASS' if ok else 'FAIL'} {method} {' '.join(options)}: {detail}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
