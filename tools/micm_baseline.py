#!/usr/bin/env python3
"""Time and score MICM's own scenario driver at the experimenter's baseline settings.

    tools/micm_baseline.py --scenarios DIR [DIR ...] [--rounds 3] [--threads 10]
                           [--binary bin/micm_scenario] [--out scans/MICM-baselines.csv]

Runs MICM's three-stage Rosenbrock (32-cell groups) at the two settings the
ranking uses as its baselines, relative 1e-2 and 3e-4 with absolute 1e-15,
on each scenario, ROUNDS times in rotating order, scores each output with
scenarios/score.py against the scenario's reference (once per distinct
output), and writes one CSV row per run in tools/scan.py's spirit, so that
RESULTS.md can set the operating points against MICM on the same machine
and session. Standard library only.
"""
import argparse, csv, datetime, hashlib, json, re, subprocess, sys, tempfile, time
from pathlib import Path

ENTRY = Path(__file__).resolve().parents[1]
ROOT = ENTRY.parent  # the Tier 3 run root: scenarios/ and the scorer live there
SETTINGS = [("micm-L2", "1e-2", "1e-15"), ("micm-L3", "3e-4", "1e-15")]
FIELDS = ["label", "rtol", "atol", "scenario", "cells", "round", "threads", "elapsed_ms", "process_s",
          "internal_steps", "accepted", "rejected", "function_calls", "jacobian_updates", "decompositions", "solves",
          "digits", "worst_species", "worst_digits", "gate", "status", "output_sha256", "started"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scenarios", type=Path, nargs="+", required=True)
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--threads", type=int, default=10)
    parser.add_argument("--binary", type=Path, default=ENTRY / "bin" / "micm_scenario")
    parser.add_argument("--score", type=Path, default=ROOT / "scenarios" / "score.py")
    parser.add_argument("--out", type=Path, default=ENTRY / "scans" / "MICM-baselines.csv")
    parser.add_argument("--pause", type=float, default=1.0)
    args = parser.parse_args()
    scored = {}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as table:
        writer = csv.DictWriter(table, fieldnames=FIELDS)
        writer.writeheader()
        for round_number in range(1, args.rounds + 1):
            for scenario in args.scenarios:
                scenario = scenario.resolve()
                reference = scenario / "reference_snapshots.csv"
                if not reference.exists():
                    reference = scenario / "reference_final.csv"
                output_name = "ts1_snapshots.csv" if reference.name.startswith("reference_snapshots") else "ts1_final.csv"
                cells = sum(1 for _ in open(scenario / "cells.csv")) - 1
                order = SETTINGS[(round_number - 1) % len(SETTINGS):] + SETTINGS[:(round_number - 1) % len(SETTINGS)]
                for label, rtol, atol in order:
                    with tempfile.TemporaryDirectory(prefix="micm-", dir=ENTRY / "build") as work:
                        work = Path(work)
                        cmd = [str(args.binary), "--threads", str(args.threads), "--group-size", "32",
                               "--rtol", rtol, "--atol", atol, str(scenario / "scenario.txt")]
                        started = datetime.datetime.now().isoformat(timespec="seconds")
                        t0 = time.perf_counter()
                        run = subprocess.run(cmd, cwd=work, capture_output=True, text=True)
                        process_s = time.perf_counter() - t0
                        status = f"exit {run.returncode}" if run.returncode else "ok"
                        values = dict(re.findall(r"(\w+)=([0-9.]+)", run.stdout))
                        row = {k: values.get(k, "") for k in FIELDS}
                        row.update({"label": label, "rtol": rtol, "atol": atol, "scenario": scenario.name, "cells": cells,
                                    "round": round_number, "threads": args.threads, "process_s": f"{process_s:.3f}",
                                    "status": status, "started": started, "digits": "", "worst_species": "",
                                    "worst_digits": "", "gate": "FAIL", "output_sha256": ""})
                        result = work / output_name
                        if status == "ok" and result.exists():
                            sha = sha256(result)
                            row["output_sha256"] = sha
                            if sha not in scored:
                                check = subprocess.run([sys.executable, str(args.score), str(result), str(reference), "--json"],
                                                       capture_output=True, text=True)
                                try:
                                    summary = json.loads(check.stdout.strip().splitlines()[-1])
                                except (IndexError, ValueError):
                                    summary = {}
                                worst = summary.get("worst_species", [["", ""]])[0]
                                scored[sha] = (summary.get("overall_digits", ""), worst[0], worst[1], summary.get("gate", "FAIL"))
                            row["digits"], row["worst_species"], row["worst_digits"], row["gate"] = scored[sha]
                        writer.writerow(row)
                        table.flush()
                        print(f"round {round_number} {scenario.name:<20} {label:<8} elapsed_ms={row['elapsed_ms']:>10} "
                              f"digits={row['digits']} worst={row['worst_species']} {row['worst_digits']} gate={row['gate']} "
                              f"{status if status != 'ok' else ''}", flush=True)
                        time.sleep(args.pause)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
