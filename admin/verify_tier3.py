#!/usr/bin/env python3
"""Verify a finished Tier 3 entry: the Tier 1 gate, its operating points, and its build.

    admin/verify_tier3.py RUN BRANCH WORK --points POINTS
                          [--scenarios DIR ...] [--compilers clang,gcc-15]
                          [--threads 10] [--stages build,tests,ts1,points,audit]

RUN is the run name (sol-xhigh) and BRANCH its collected branch in this
repository (tier3/sol/xhigh). WORK is a verification directory, made if missing.
For each compiler, WORK/<cc>/ gets links micm, reference, and scenarios to the
Tier 3 run root's, and a detached worktree of BRANCH at WORK/<cc>/RUN. Stages:

  build   make env, then make cpu bench CC=<cc> OPTFLAGS="-O3 -DNDEBUG"
  tests   make test-cpu; its "<id> PASS|FAIL|N/A" lines are counted against
          reference/micm_cases.txt, every id exactly once
  ts1     the Tier 1 TS1 gate: bin/ts1_bench_cpu on 10,000 cells of seed 2026,
          100 steps of 30 s, checked against MICM's final state with
          reference/compare_final.py (1e-10 relative). MICM's driver, the cells,
          and its final state are made once, in WORK/../tier1-gate/
  points  every declared operating point on every scenario: bin/ts1_scenario
          --threads THREADS OPTIONS SCENARIO, scored with tier3/scenarios/score.py
          against the scenario's reference, at levels 2 and 3
  audit   the compile and link flags (make -n -B), the libraries each program
          links (otool -L), and the lines of C that match a banned construct

POINTS lists the declared operating points, one per line, as NAME<TAB>OPTIONS,
OPTIONS being bin/ts1_scenario's arguments before the scenario; # starts a
comment. --scenarios takes directories of scenario directories (default: the
disclosed set). The points stage runs with the first compiler only, the build
the ranking times, unless --points-all-compilers is given.

WORK/report.md collects every stage's results; the raw logs are in WORK/logs/
and the scenario runs in WORK/runs/. Needs only the Python standard library.
"""

import argparse
import json
import os
import re
import shlex
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
TIER3_ROOT = Path.home() / "NSF-NCAR/ODE/benchmarks/tier3"
OPTFLAGS = "-O3 -DNDEBUG"
STAGES = ("build", "tests", "ts1", "points", "audit")

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("run")
parser.add_argument("branch")
parser.add_argument("work", type=Path)
parser.add_argument("--points", type=Path, help="the declared operating points (needed for the points stage)")
parser.add_argument("--scenarios", type=Path, nargs="+", default=[TIER3_ROOT / "scenarios"])
parser.add_argument("--compilers", default="clang,gcc-15")
parser.add_argument("--threads", type=int, default=10)
parser.add_argument("--stages", default=",".join(STAGES))
parser.add_argument("--points-all-compilers", action="store_true")
parser.add_argument("--timeout", type=float, default=3600.0, help="seconds allowed for one scenario run")
args = parser.parse_args()
stages = args.stages.split(",")
unknown = set(stages) - set(STAGES)
if unknown:
    parser.error(f"unknown stages: {' '.join(sorted(unknown))}")
compilers = args.compilers.split(",")
work = args.work.resolve()
logs = work / "logs"
logs.mkdir(parents=True, exist_ok=True)
report = []


def sdk_env():
    env = dict(os.environ)
    if sys.platform == "darwin" and "SDKROOT" not in env:
        env["SDKROOT"] = subprocess.run(["xcrun", "--sdk", "macosx", "--show-sdk-path"],
                                        capture_output=True, text=True).stdout.strip()
    return env


def run(command, cwd, log, timeout=None):
    """Run COMMAND in CWD, writing its output to LOG; return (exit status, output, seconds)."""
    start = time.perf_counter()
    try:
        done = subprocess.run(command, cwd=cwd, capture_output=True, text=True, timeout=timeout, env=sdk_env(),
                              shell=isinstance(command, str))
        status, output = done.returncode, done.stdout + done.stderr
    except subprocess.TimeoutExpired as expired:
        status, output = "timeout", f"{expired.stdout or ''}{expired.stderr or ''}\nTIMEOUT after {timeout} s"
    seconds = time.perf_counter() - start
    log.write_text(f"$ {command if isinstance(command, str) else shlex.join(map(str, command))}\n"
                   f"# cwd {cwd}; exit {status}; {seconds:.1f} s\n{output}")
    return status, output, seconds


def link(path, target):
    # -n semantics: never follow an existing link (a racing ln -s writes inside the target)
    if path.is_symlink():
        if os.readlink(path) == str(target):
            return
        path.unlink()
    os.symlink(target, path)


def tree(cc):
    return work / cc / args.run


def prepare(cc):
    base = work / cc
    base.mkdir(parents=True, exist_ok=True)
    for name in ("micm", "reference", "scenarios"):
        link(base / name, (TIER3_ROOT / name).resolve())
    if not tree(cc).exists():
        subprocess.run(["git", "-C", str(REPO), "worktree", "add", "--detach", str(tree(cc)), args.branch],
                       check=True, capture_output=True)
    head = subprocess.run(["git", "-C", str(tree(cc)), "log", "-1", "--format=%h %s"],
                          capture_output=True, text=True).stdout.strip()
    return head


def make(cc, target, log, timeout=None):
    return run(["make", target, f"CC={cc}", f"OPTFLAGS={OPTFLAGS}"], tree(cc), log, timeout)


heads = {cc: prepare(cc) for cc in compilers}
report.append(f"# Verification of tier3/{args.run}\n")
report.append(f"Branch `{args.branch}`, commit {heads[compilers[0]]}; compilers {', '.join(compilers)}; "
              f"OPTFLAGS `{OPTFLAGS}`; {time.strftime('%Y-%m-%d %H:%M')}.\n")

# build
if "build" in stages:
    report.append("## Build\n")
    report.append("| Compiler | make env | make cpu bench |")
    report.append("|---|---|---|")
    for cc in compilers:
        env_status, _, env_s = run(["make", "env"], tree(cc), logs / f"{cc}-env.log")
        build = [make(cc, target, logs / f"{cc}-{target}.log") for target in ("cpu", "bench")]
        statuses = " ".join(f"{t} exit {b[0]} ({b[2]:.0f} s)" for t, b in zip(("cpu", "bench"), build))
        warnings = sum(output.count("warning:") for _, output, _ in build)
        report.append(f"| {cc} | exit {env_status} ({env_s:.0f} s) | {statuses}; {warnings} warnings |")
    report.append("")

# tests
if "tests" in stages:
    cases = [line.split()[0] for line in (TIER3_ROOT / "reference" / "micm_cases.txt").read_text().splitlines()
             if line.strip()]
    case_set = set(cases)
    report.append("## MICM's test cases\n")
    report.append("| Compiler | make test-cpu | PASS | FAIL | N/A | Ids missing | Ids repeated | Unknown ids |")
    report.append("|---|---|---:|---:|---:|---:|---:|---:|")
    failures = {}
    for cc in compilers:
        status, output, seconds = make(cc, "test-cpu", logs / f"{cc}-test-cpu.log")
        seen = {}
        for line in output.splitlines():
            found = re.match(r"^\s*(\S+)\s+(PASS|FAIL|N/A)\b", line)
            # a listed id, or one shaped like a GoogleTest case id that the list lacks
            if found and (found.group(1) in case_set or re.fullmatch(r"[\w-]+:[\w/]+\.[\w/]+", found.group(1))):
                seen.setdefault(found.group(1), []).append(found.group(2))
        counts = {verdict: sum(1 for v in seen.values() for x in v if x == verdict) for verdict in ("PASS", "FAIL", "N/A")}
        missing = [c for c in cases if c not in seen]
        repeated = [c for c, v in seen.items() if len(v) > 1]
        unknown_ids = [c for c in seen if c not in case_set]
        failures[cc] = [c for c, v in seen.items() if "FAIL" in v]
        report.append(f"| {cc} | exit {status} ({seconds:.0f} s) | {counts['PASS']} | {counts['FAIL']} | {counts['N/A']} "
                      f"| {len(missing)} | {len(repeated)} | {len(unknown_ids)} |")
    report.append("")
    for cc, failed in failures.items():
        if failed:
            report.append(f"{cc} failures: {' '.join(failed[:40])}\n")

# ts1
if "ts1" in stages:
    gate_dir = work.parent / "tier1-gate"
    gate_dir.mkdir(exist_ok=True)
    cells = gate_dir / "cells_2026.csv"
    micm_final = gate_dir / "micm_final_2026.csv"
    if not cells.exists():
        run([sys.executable, str(TIER3_ROOT / "reference" / "make_cells_csv.py"), "10000", "2026", "--output", str(cells)],
            gate_dir, gate_dir / "make_cells.log")
    if not micm_final.exists():
        driver = gate_dir / "micm_ts1_ref"
        if not driver.exists():
            run(["clang++", "-std=c++20", *OPTFLAGS.split(), "-I", str(TIER3_ROOT / "micm/include"),
                 "-I", str(TIER3_ROOT / "micm/test/benchmark"), str(TIER3_ROOT / "reference/micm_ts1_ref.cpp"),
                 "-o", str(driver)], gate_dir, gate_dir / "build_micm_ts1_ref.log")
        status, _, _ = run([str(driver), "10000", "100", "30", str(cells)], gate_dir, gate_dir / "micm_ts1_ref.log")
        if status == 0:
            (gate_dir / "ts1_final.csv").rename(micm_final)
    report.append("## Tier 1 TS1 gate (10,000 cells, seed 2026, 100 steps of 30 s)\n")
    report.append("| Compiler | ts1_bench_cpu | elapsed_ms | compare_final.py (1e-10) |")
    report.append("|---|---|---:|---|")
    for cc in compilers:
        out = work / "runs" / cc / "tier1-ts1"
        out.mkdir(parents=True, exist_ok=True)
        status, output, _ = run([str(tree(cc) / "bin/ts1_bench_cpu"), "10000", "100", "30", str(cells)], out,
                                logs / f"{cc}-ts1_bench_cpu.log")
        elapsed = re.search(r"elapsed_ms=([0-9.]+)", output)
        verdict = "no result"
        if status == 0 and (out / "ts1_final.csv").exists():
            check_status, check, _ = run([sys.executable, str(TIER3_ROOT / "reference/compare_final.py"),
                                          str(out / "ts1_final.csv"), str(micm_final)], out, logs / f"{cc}-ts1-compare.log")
            verdict = ("PASS" if check_status == 0 else "FAIL") + ": " + " ".join(check.split())[:160]
        report.append(f"| {cc} | exit {status} | {elapsed.group(1) if elapsed else ''} | {verdict} |")
    report.append("")

# points
if "points" in stages:
    if args.points is None:
        parser.error("the points stage needs --points")
    points = []
    for line in args.points.read_text().splitlines():
        line = line.split("#", 1)[0].rstrip()
        if line.strip():
            name, _, options = line.partition("\t")
            points.append((re.sub(r"[^\w.+-]", "_", name.strip()), shlex.split(options)))
    scenarios = sorted(d for base in args.scenarios for d in base.iterdir() if (d / "scenario.txt").exists())
    report.append(f"## Operating points on the scenarios (--threads {args.threads}, 1,000 cells)\n")
    for name, options in points:
        report.append(f"- `{name}`: `bin/ts1_scenario {shlex.join(options)} SCENARIO`")
    report.append("")
    report.append("| Compiler | Point | Scenario | Exit | elapsed_ms | Digits | Worst species | Gate | Level 2 | Level 3 |")
    report.append("|---|---|---|---|---:|---:|---|---|---|---|")
    rows = []
    for cc in (compilers if args.points_all_compilers else compilers[:1]):
        for name, options in points:
            for scenario in scenarios:
                out = work / "runs" / cc / name / scenario.name
                out.mkdir(parents=True, exist_ok=True)
                for stale in ("ts1_final.csv", "ts1_snapshots.csv"):
                    (out / stale).unlink(missing_ok=True)
                status, output, _ = run([str(tree(cc) / "bin/ts1_scenario"), "--threads", str(args.threads), *options,
                                         str(scenario / "scenario.txt")], out, out / "run.log", args.timeout)
                elapsed = re.search(r"elapsed_ms=([0-9.]+)", output)
                reference = scenario / "reference_snapshots.csv"
                result = out / "ts1_snapshots.csv"
                if not reference.exists():
                    reference, result = scenario / "reference_final.csv", out / "ts1_final.csv"
                row = {"compiler": cc, "point": name, "scenario": scenario.name, "exit": status,
                       "elapsed_ms": elapsed.group(1) if elapsed else ""}
                if status == 0 and result.exists():
                    for level in (2, 3):
                        _, scored, _ = run([sys.executable, str(REPO / "tier3/scenarios/score.py"), str(result),
                                            str(reference), "--level", str(level), "--json"], out, out / f"score-{level}.log")
                        try:
                            summary = json.loads(scored.strip().splitlines()[-1])
                        except (IndexError, ValueError):
                            summary = {"gate": "FAIL", "level_result": "FAIL", "error": scored.strip()[-200:]}
                        row["digits"] = summary.get("overall_digits", "")
                        row["worst"] = " ".join(f"{n} {v}" for n, v in summary.get("worst_species", [])[:1])
                        row["gate"] = summary.get("gate", "FAIL")
                        row[f"level{level}"] = summary.get("level_result", "FAIL")
                else:
                    row.update(digits="", worst="", gate="FAIL", level2="FAIL", level3="FAIL")
                rows.append(row)
                report.append(f"| {cc} | {name} | {scenario.name} | {status} | {row['elapsed_ms']} | {row['digits']} "
                              f"| {row['worst']} | {row['gate']} | {row['level2']} | {row['level3']} |")
                print(f"{cc} {name} {scenario.name}: exit {status} digits {row['digits']} gate {row['gate']}", flush=True)
    (work / "points.json").write_text(json.dumps(rows, indent=1))
    report.append("")
    for name, _ in points:
        mine = [r for r in rows if r["point"] == name and r["compiler"] == compilers[0]]
        verdicts = {
            "gate on every scenario": all(r["gate"] == "PASS" for r in mine),
            "level 2 on every scenario": all(r["level2"] == "PASS" for r in mine),
            "level 3 on every scenario": all(r["level3"] == "PASS" for r in mine),
        }
        report.append(f"- `{name}`: " + "; ".join(f"{k}: {'yes' if v else 'NO'}" for k, v in verdicts.items()))
    report.append("")

# audit
if "audit" in stages:
    cc = compilers[0]
    report.append("## Build audit\n")
    status, output, _ = run(["make", "-n", "-B", "cpu", "bench", f"CC={cc}", f"OPTFLAGS={OPTFLAGS}"], tree(cc),
                            logs / "make-n.log")
    def tokens(line):
        try:
            return shlex.split(line)
        except ValueError:
            return line.split()
    flags = sorted({token for line in output.splitlines() for token in tokens(line)
                    if token.startswith("-") and token != "-o"} if status == 0 else [])
    suspicious = [f for f in flags if re.match(r"-(march|mcpu|mtune|ffast-math|funsafe|fno-math|ffp-|fassociative|"
                                               r"freciprocal|ffinite|fno-signed|fno-trapping|flto|fprofile|fauto-profile|"
                                               r"Ofast|O[0124sz]?$|fvectorize|fslp|funroll|mllvm)", f)]
    report.append(f"`make -n -B cpu bench CC={cc}`: exit {status}; flags used: `{' '.join(flags)}`\n")
    report.append(f"Flags outside the rules (review): {', '.join(suspicious) or 'none'}\n")
    for program in sorted((tree(cc) / "bin").glob("*")):
        if program.is_file() and os.access(program, os.X_OK):
            libs = subprocess.run(["otool", "-L", str(program)], capture_output=True, text=True).stdout.split("\n")[1:]
            report.append(f"- `bin/{program.name}` links: " + ", ".join(l.split()[0] for l in libs if l.strip()))
    report.append("")
    envs = {p.parent for p in tree(cc).rglob("conda-meta")}
    sources = [p for p in tree(cc).rglob("*") if p.suffix in (".c", ".h") and p.is_file()
               and not any(env in p.parents for env in envs)]
    banned = re.compile(r"__attribute__|__builtin_|\basm\b|__asm|vector_size|ext_vector_type|long\s+double|"
                        r"#\s*include\s*<(arm_|immintrin|xmmintrin|emmintrin|x86intrin|Accelerate|dispatch|mach/)|"
                        r"#\s*pragma\s+(clang|GCC)|omp\s+target|__APPLE__|__aarch64__|__x86_64__|__ARM_NEON|"
                        r"_Float16|__fp16|__bf16")
    matches = []
    for source in sources:
        for number, line in enumerate(source.read_text(errors="replace").splitlines(), 1):
            if banned.search(line):
                matches.append(f"{source.relative_to(tree(cc))}:{number}: {line.strip()[:140]}")
    (logs / "banned.txt").write_text("\n".join(matches) + "\n")
    report.append(f"C sources scanned: {len(sources)} (tracked and generated); lines matching a banned construct: "
                  f"{len(matches)} (all in `logs/banned.txt`)\n")
    for match in matches[:30]:
        report.append(f"    {match}")
    report.append("")

(work / "report.md").write_text("\n".join(report) + "\n")
print("\n".join(report))
