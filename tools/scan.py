#!/usr/bin/env python3
"""Run a structure scan: a matrix of ts1_scenario configurations on scenarios.

    tools/scan.py NAME --points POINTS_FILE --scenarios DIR [DIR ...]
                  [--rounds 1] [--threads 10] [--binary bin/ts1_scenario]
                  [--out scans/NAME.csv] [--pause 1]

POINTS_FILE has one configuration per line: LABEL<TAB>METHOD<TAB>OPTIONS,
OPTIONS being space-separated KEY=VALUE pairs for --option ('#' starts a
comment). Every configuration runs on every scenario directory (one with
scenario.txt and a reference), ROUNDS times, in rotating order so that no
configuration always runs first. Each run is scored with scenarios/score.py
against the scenario's reference (once per distinct output, by SHA-256), and
one row goes to the CSV: label, method, options, scenario, cells, round,
elapsed_ms, process_s, the solver counts, digits, worst species and its
digits, gate, and the output's hash. Outputs are deleted after scoring.
Runs in the conda environment's Python; needs only the standard library.
"""
import argparse
import csv
import datetime
import hashlib
import json
import re
import shlex
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ENTRY = Path(__file__).resolve().parent.parent
ROOT = ENTRY.parent
COUNTS = ["internal_steps", "accepted", "rejected", "function_calls", "jacobian_updates", "decompositions", "solves", "matrix_reuses", "rate_updates"]


def read_points(path):
    points = []
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            sys.exit(f"{path}: expected LABEL<TAB>METHOD<TAB>OPTIONS, got {line!r}")
        label, method = parts[0].strip(), parts[1].strip()
        options = parts[2].split() if len(parts) > 2 else []
        points.append((label, method, options))
    return points


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("name")
    parser.add_argument("--points", type=Path, required=True)
    parser.add_argument("--scenarios", type=Path, nargs="+", required=True)
    parser.add_argument("--rounds", type=int, default=1)
    parser.add_argument("--threads", type=int, default=10)
    parser.add_argument("--binary", type=Path, default=ENTRY / "bin" / "ts1_scenario")
    parser.add_argument("--score", type=Path, default=ROOT / "scenarios" / "score.py")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--pause", type=float, default=1.0)
    parser.add_argument("--timeout", type=float, default=1800.0)
    args = parser.parse_args()
    out = args.out or ENTRY / "scans" / f"{args.name}.csv"
    points = read_points(args.points)
    scenarios = [d for d in args.scenarios if (d / "scenario.txt").exists()]
    if not scenarios:
        sys.exit("no scenario directories with scenario.txt")
    fields = ["label", "method", "options", "scenario", "cells", "round", "threads", "elapsed_ms", "process_s", *COUNTS,
              "digits", "worst_species", "worst_digits", "gate", "status", "output_sha256", "started"]
    new = not out.exists()
    out.parent.mkdir(parents=True, exist_ok=True)
    scored = {}  # sha -> (digits, worst, worst_digits, gate)
    with open(out, "a", newline="") as table:
        writer = csv.DictWriter(table, fieldnames=fields)
        if new:
            writer.writeheader()
        for round_number in range(1, args.rounds + 1):
            shift = (round_number - 1) * len(points) // args.rounds % len(points)
            order = points[shift:] + points[:shift]
            for scenario in scenarios:
                reference = scenario / "reference_snapshots.csv"
                if not reference.exists():
                    reference = scenario / "reference_final.csv"
                output_name = "ts1_snapshots.csv" if reference.name.startswith("reference_snapshots") else "ts1_final.csv"
                cells = sum(1 for _ in open(scenario / "cells.csv")) - 1
                for label, method, options in order:
                    with tempfile.TemporaryDirectory(prefix="scan-", dir=ENTRY / "build") as work:
                        work = Path(work)
                        cmd = [str(args.binary), "--threads", str(args.threads), "--method", method,
                               *sum((["--option", o] for o in options), []), str(scenario / "scenario.txt")]
                        started = datetime.datetime.now().isoformat(timespec="seconds")
                        t0 = time.perf_counter()
                        try:
                            run = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=args.timeout)
                            status = f"exit {run.returncode}" if run.returncode else "ok"
                            stdout = run.stdout
                        except subprocess.TimeoutExpired:
                            status, stdout = "timeout", ""
                        process_s = time.perf_counter() - t0
                        values = dict(re.findall(r"(\w+)=([0-9.]+)", stdout))
                        row = {"label": label, "method": method, "options": " ".join(options), "scenario": scenario.name, "cells": cells,
                               "round": round_number, "threads": args.threads, "elapsed_ms": values.get("elapsed_ms", ""),
                               "process_s": f"{process_s:.3f}", "status": status, "started": started,
                               "digits": "", "worst_species": "", "worst_digits": "", "gate": "FAIL", "output_sha256": ""}
                        for k in COUNTS:
                            row[k] = values.get(k, "")
                        result = work / output_name
                        if status == "ok" and result.exists():
                            sha = sha256(result)
                            row["output_sha256"] = sha
                            if sha not in scored:
                                check = subprocess.run([sys.executable, str(args.score), str(result), str(reference), "--json"], capture_output=True, text=True)
                                try:
                                    summary = json.loads(check.stdout.strip().splitlines()[-1])
                                except (IndexError, ValueError):
                                    summary = {}
                                worst = summary.get("worst_species", [["", ""]])[0]
                                scored[sha] = (summary.get("overall_digits", ""), worst[0], worst[1], summary.get("gate", "FAIL"))
                            row["digits"], row["worst_species"], row["worst_digits"], row["gate"] = scored[sha]
                        writer.writerow(row)
                        table.flush()
                        print(f"round {round_number} {scenario.name:<20} {label:<28} elapsed_ms={row['elapsed_ms']:>10} digits={row['digits']} worst={row['worst_species']} {row['worst_digits']} gate={row['gate']} {status if status != 'ok' else ''}", flush=True)
                        time.sleep(args.pause)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
