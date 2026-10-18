#!/usr/bin/env python3
"""Time TS1 benchmark programs the same way, in rotating round-robin order.

    admin/time_entries.py --input CELLS_CSV --reference MICM_FINAL --out DIR \\
        NAME=COMMAND [NAME=COMMAND ...]
    admin/time_entries.py --scenario SCENARIO --reference REFERENCE --out DIR \\
        NAME=COMMAND [NAME=COMMAND ...]

Every COMMAND is a TS1 benchmark taking `[--threads N] CELLS STEPS DT CELLS_CSV`,
like reference/micm_ts1_ref and an entry's bin/ts1_bench_cpu, and may carry its
own options, for example "micm_ts1_ref --group-size 32". Each round runs every
program at every thread count with the same arguments and input, rotating the
program order from round to round so that no program always runs first or last:
each round starts an equal share of the list further on (one program further,
when there are as many rounds as programs).
Each run gets a fresh working directory. Its ts1_final.csv is checked against
MICM_FINAL with reference/compare_final.py and then deleted.

With --scenario (Tier 3), every COMMAND takes `[--threads N] SCENARIO`, like
tier3/scenarios/micm_scenario and an entry's bin/ts1_scenario, with its method
and options in the COMMAND. Each run is scored with tier3/scenarios/score.py
against REFERENCE, a scenario's reference_snapshots.csv (scoring the run's
ts1_snapshots.csv) or reference_final.csv (scoring its ts1_final.csv), and the
gate is score.py's. A 10,000-cell scenario can be scored against the 1,000-cell
reference of the same arguments: score.py scores only the reference's cells,
which are the first 1,000, and applies the gate to every cell. Scoring a
10,000-cell output takes much longer than the run, so each program's first
complete run is scored, and a later run only if its output's SHA-256 differs
from that run's (--score-every-run scores them all).

DIR/timings.csv gets one row per run: the program's own elapsed_ms and solver
statistics, the whole process's wall time (which also counts reading the input,
setup, and writing the output), the gate result, and, with --scenario, the
overall digits, the worst species, and the output's SHA-256. DIR/summary.md
gives the median and range of elapsed_ms for each program and thread count, and
with --scenario the median digits; DIR/environment.txt records the machine and
the commands. Needs only the Python standard library.
"""

import argparse
import csv
import datetime
import hashlib
import json
import platform
import re
import shlex
import statistics
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
STATS = ("internal_steps", "accepted", "rejected", "function_calls", "jacobian_updates", "decompositions", "solves")

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("programs", nargs="+", metavar="NAME=COMMAND")
parser.add_argument("--input", type=Path, help="the CELLS_CSV every program reads (Tier 1)")
parser.add_argument("--scenario", type=Path, help="the scenario file every program runs (Tier 3)")
parser.add_argument("--reference", type=Path, required=True,
                    help="MICM's ts1_final.csv for the same input, or the scenario's reference_final.csv or reference_snapshots.csv")
parser.add_argument("--out", type=Path, required=True)
parser.add_argument("--cells", type=int, default=10000)
parser.add_argument("--steps", type=int, default=100)
parser.add_argument("--dt", default="30")
parser.add_argument("--threads", default="10,1", help="comma-separated thread counts (default: 10,1)")
parser.add_argument("--rounds", type=int, default=5)
parser.add_argument("--pause", type=float, default=5.0, help="seconds between runs (default: 5)")
parser.add_argument("--compare", type=Path, default=ROOT / "reference" / "compare_final.py")
parser.add_argument("--score", type=Path, default=ROOT / "tier3" / "scenarios" / "score.py")
parser.add_argument("--score-every-run", action="store_true",
                    help="with --scenario, score every run, not only each program's first and any whose output differs")
args = parser.parse_args()
if (args.input is None) == (args.scenario is None):
    parser.error("give exactly one of --input and --scenario")
if args.scenario:
    with open(args.reference, newline="") as file:
        snapshots = next(csv.reader(file))[:2] == ["step", "cell"]
    output_name = "ts1_snapshots.csv" if snapshots else "ts1_final.csv"

programs = []
for item in args.programs:
    name, _, command = item.partition("=")
    if not name or not command:
        parser.error(f"expected NAME=COMMAND, got {item!r}")
    programs.append((name, shlex.split(command)))
thread_counts = [int(t) for t in args.threads.split(",")]


def environment():
    lines = [
        f"date: {datetime.datetime.now().isoformat(timespec='seconds')}",
        f"host: {platform.node()}",
        f"platform: {platform.platform()}",
    ]
    for command in (
        ["sysctl", "-n", "hw.model", "machdep.cpu.brand_string", "hw.perflevel0.physicalcpu", "hw.perflevel1.physicalcpu"],
        ["pmset", "-g", "batt"],
        ["uptime"],
    ):
        try:
            result = subprocess.run(command, capture_output=True, text=True)
            lines.append(f"{shlex.join(command)}: {' | '.join(result.stdout.split(chr(10))).strip(' |')}")
        except OSError:
            pass
    if args.scenario:
        lines += [f"scenario: {args.scenario.resolve()}", f"reference: {args.reference.resolve()} (scores {output_name})"]
        lines += [f"threads: {args.threads}; rounds: {args.rounds}"]
    else:
        lines += [f"input: {args.input.resolve()}", f"reference: {args.reference.resolve()}"]
        lines += [f"cells, steps, dt: {args.cells}, {args.steps}, {args.dt}; threads: {args.threads}; rounds: {args.rounds}"]
    lines += [f"program {name}: {shlex.join(command)}" for name, command in programs]
    return lines


args.out.mkdir(parents=True, exist_ok=True)
(args.out / "environment.txt").write_text("\n".join(environment()) + "\n")
table = args.out / "timings.csv"
fields = ["round", "threads", "program", "elapsed_ms", *STATS, "process_s", "gate", "max_relative_difference", "started",
          "digits", "worst_species", "output_sha256"]
with open(table, "w", newline="") as file:
    csv.writer(file).writerow(fields)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as file:
        for block in iter(lambda: file.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


first_score = {}  # program -> (SHA-256, gate, digits, worst, run directory) of its first scored run

for round_number in range(1, args.rounds + 1):
    shift = (round_number - 1) * len(programs) // args.rounds % len(programs)
    order = programs[shift:] + programs[:shift]
    for threads in thread_counts:
        for name, command in order:
            work = args.out / "runs" / f"{name}-t{threads}-r{round_number}"
            work.mkdir(parents=True, exist_ok=True)
            started = datetime.datetime.now().isoformat(timespec="seconds")
            if args.scenario:
                argv = [*command, "--threads", str(threads), str(args.scenario.resolve())]
            else:
                argv = [*command, "--threads", str(threads), str(args.cells), str(args.steps), args.dt, str(args.input.resolve())]
            begin = time.perf_counter()
            run = subprocess.run(argv, cwd=work, capture_output=True, text=True)
            process_s = time.perf_counter() - begin
            (work / "output.txt").write_text(run.stdout + run.stderr)
            values = dict(re.findall(r"(\w+)=([0-9.]+)", run.stdout))
            final = work / "ts1_final.csv"
            digits = worst = sha = ""
            if args.scenario:
                output = work / output_name
                if run.returncode == 0 and output.exists():
                    sha = sha256(output)
                    known = first_score.get(name)
                    if known and known[0] == sha and not args.score_every_run:
                        gate, digits, worst = known[1:4]
                        (work / "score.txt").write_text(f"output identical to {known[4]}'s (SHA-256 {sha}); its score is reused\n")
                    else:
                        check = subprocess.run(
                            [sys.executable, str(args.score), str(output), str(args.reference.resolve()), "--json"],
                            capture_output=True,
                            text=True,
                        )
                        (work / "score.txt").write_text(check.stdout + check.stderr)
                        try:
                            scored = json.loads(check.stdout.strip().splitlines()[-1])
                        except (IndexError, ValueError):
                            scored = {}
                        gate = scored.get("gate", "FAIL")
                        digits = scored.get("overall_digits", "")
                        worst = " ".join(f"{species} {value}" for species, value in scored.get("worst_species", [])[:1])
                        first_score.setdefault(name, (sha, gate, digits, worst, work.name))
                else:
                    gate = f"ERROR {run.returncode}"
                largest = ""
                for filename in ("ts1_final.csv", "ts1_snapshots.csv"):
                    (work / filename).unlink(missing_ok=True)
            elif run.returncode == 0 and final.exists():
                check = subprocess.run(
                    [sys.executable, str(args.compare), str(final), str(args.reference.resolve())],
                    capture_output=True,
                    text=True,
                )
                (work / "compare.txt").write_text(check.stdout + check.stderr)
                gate = "PASS" if check.returncode == 0 else "FAIL"
                found = re.search(r"largest relative difference (\S+)", check.stdout)
                largest = found.group(1) if found else ""
                final.unlink()
            else:
                gate, largest = f"ERROR {run.returncode}", ""
            row = [round_number, threads, name, values.get("elapsed_ms", ""), *(values.get(s, "") for s in STATS)]
            row += [f"{process_s:.3f}", gate, largest, started, digits, worst, sha]
            with open(table, "a", newline="") as file:
                csv.writer(file).writerow(row)
            print(
                f"round {round_number} threads {threads:>2} {name:<12} elapsed_ms={values.get('elapsed_ms', '?'):>12}"
                f" process_s={process_s:8.3f} gate={gate}" + (f" digits={digits}" if args.scenario else f" largest={largest}"),
                flush=True,
            )
            time.sleep(args.pause)

with open(table, newline="") as file:
    rows = list(csv.DictReader(file))
lines = [
    "| Program | Threads | Runs passing the gate | Median elapsed_ms | Range | Median process wall time (s) |"
    + (" Median digits |" if args.scenario else ""),
    "|---|---:|---:|---:|---:|---:|" + ("---:|" if args.scenario else ""),
]
for threads in thread_counts:
    for name, _ in programs:
        runs = [row for row in rows if row["program"] == name and int(row["threads"]) == threads]
        good = [row for row in runs if row["gate"] == "PASS" and row["elapsed_ms"]]
        if not good:
            lines.append(f"| {name} | {threads} | 0 of {len(runs)} | | | |" + (" |" if args.scenario else ""))
            continue
        times = [float(row["elapsed_ms"]) for row in good]
        walls = [float(row["process_s"]) for row in good]
        line = (
            f"| {name} | {threads} | {len(good)} of {len(runs)} | {statistics.median(times):,.1f} |"
            f" {min(times):,.1f}–{max(times):,.1f} | {statistics.median(walls):,.1f} |"
        )
        if args.scenario:
            digits = [float(row["digits"]) for row in good if row["digits"]]
            line += f" {statistics.median(digits):.2f} |" if digits else " |"
        lines.append(line)
(args.out / "summary.md").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
