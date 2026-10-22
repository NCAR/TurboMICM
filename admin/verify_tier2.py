#!/usr/bin/env python3
"""Verify a finished Tier 2 entry: the Tier 1 gate, the specialized solver, and the build.

    admin/verify_tier2.py RUN BRANCH WORK [--compilers clang,gcc-15] [--threads 18]
                          [--scenarios DIR ...] [--ts1-input CELLS MICM_FINAL ...]
                          [--stages build,tests,ts1,specialized,generator,daynight,audit]

RUN is the run name (astra-max) and BRANCH its collected branch in this
repository (tier2/astra/max). WORK is a verification directory, made if missing.
For each compiler, WORK/<cc>/ gets links micm, reference, and scenarios to the
Tier 2 run root's, and a detached worktree of BRANCH at WORK/<cc>/RUN. Stages:

  build        make env, then make cpu bench CC=<cc> OPTFLAGS="-O3 -DNDEBUG" in the
               fresh worktree, timed: the clean-build time. The C files the build
               creates or rewrites are the generated C, whose size is reported
  tests        make test-cpu; its "<id> PASS|FAIL|N/A" lines are counted against
               reference/micm_cases.txt, every id exactly once
  ts1          the Tier 1 TS1 gate: bin/ts1_bench_cpu on 10,000 cells of seed
               2026, 100 steps of 30 s, against MICM's final state with
               reference/compare_final.py (1e-10 relative). MICM's driver, the
               cells, and its final state are made once, in WORK/../tier1-gate/
  specialized  the same for bin/ts1_bench_specialized, on the seed-2026 cells and
               on each --ts1-input pair (a CELLS_CSV and MICM's final state for it)
  generator    make test-specialized: the entry's tests of its generator
  daynight     bin/ts1_scenario_specialized --threads THREADS --rtol 1e-6 --atol
               1e-15 --group-size 32 on every scenario (if the program rejects
               --group-size, from a second build made with GROUP_SIZE=32), scored
               with scenarios/score.py against MICM's result at the same settings
               (--level 10: 10 digits overall, 9 for every species, and the gate)
               and, for scale, against the reference solution
  audit        the compile and link flags (make -n -B), the libraries each program
               links (otool -L), and the lines of C, generated C included, that
               match a banned construct or start a process

--scenarios takes directories of scenario directories (default: the disclosed
set, ROOT/scenarios). The ts1, specialized, and daynight stages use every
compiler. WORK/report.md collects the stages' results (WORK/report-<stages>.md
for a run of only some stages), and WORK/results.json the numbers, merged across
runs; the raw logs are in WORK/logs/ and the runs in WORK/runs/. Needs only the
Python standard library.
"""

import argparse
import hashlib
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
TIER2_ROOT = Path.home() / "NSF-NCAR/ODE/benchmarks/tier2"
OPTFLAGS = "-O3 -DNDEBUG"
STAGES = ("build", "tests", "ts1", "specialized", "generator", "daynight", "audit")
DAYNIGHT = ["--rtol", "1e-6", "--atol", "1e-15"]

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("run")
parser.add_argument("branch")
parser.add_argument("work", type=Path)
parser.add_argument("--scenarios", type=Path, nargs="+", default=[TIER2_ROOT / "scenarios"])
parser.add_argument("--ts1-input", type=Path, nargs=2, action="append", default=[], metavar=("CELLS", "MICM_FINAL"))
parser.add_argument("--compilers", default="clang,gcc-15")
parser.add_argument("--threads", type=int, default=18)
parser.add_argument("--stages", default=",".join(STAGES))
parser.add_argument("--timeout", type=float, default=7200.0, help="seconds allowed for one make or one run")
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
results_path = work / "results.json"
results = json.loads(results_path.read_text()) if results_path.exists() else {}


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
        done = subprocess.run(command, cwd=cwd, capture_output=True, text=True, timeout=timeout or args.timeout,
                              env=sdk_env(), shell=isinstance(command, str))
        status, output = done.returncode, done.stdout + done.stderr
    except subprocess.TimeoutExpired as expired:
        status, output = "timeout", f"{expired.stdout or ''}{expired.stderr or ''}\nTIMEOUT after {timeout} s"
    seconds = time.perf_counter() - start
    log.parent.mkdir(parents=True, exist_ok=True)
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


def tree(cc, variant=""):
    return work / f"{cc}{variant}" / args.run


def prepare(cc, variant=""):
    base = work / f"{cc}{variant}"
    base.mkdir(parents=True, exist_ok=True)
    for name in ("micm", "reference", "scenarios"):
        link(base / name, (TIER2_ROOT / name).resolve())
    if not tree(cc, variant).exists():
        subprocess.run(["git", "-C", str(REPO), "worktree", "add", "--detach", str(tree(cc, variant)), args.branch],
                       check=True, capture_output=True)
    return subprocess.run(["git", "-C", str(tree(cc, variant)), "log", "-1", "--format=%h %s"],
                          capture_output=True, text=True).stdout.strip()


def make(cc, target, log, variant="", extra=()):
    return run(["make", target, f"CC={cc}", f"OPTFLAGS={OPTFLAGS}", *extra], tree(cc, variant), log)


def conda_roots(root):
    """Conda environments (conda-meta/) and extracted packages of a package cache (info/index.json)."""
    return {p.parent for p in root.rglob("conda-meta")} | {p.parent.parent for p in root.rglob("info/index.json")}


def c_files(root):
    skip = conda_roots(root)
    return {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in root.rglob("*")
            if p.suffix in (".c", ".h") and p.is_file() and not any(s in p.parents for s in skip)}


def elapsed_and_stats(output):
    values = dict(re.findall(r"(\w+)=([0-9.]+)", output))
    return values.get("elapsed_ms", ""), values


def compare(final, reference, log):
    status, check, _ = run([sys.executable, str(TIER2_ROOT / "reference/compare_final.py"), str(final), str(reference)],
                           final.parent, log)
    largest = re.search(r"largest relative difference (\S+)", check)
    return ("PASS" if status == 0 else "FAIL"), (largest.group(1) if largest else ""), " ".join(check.split())[:160]


heads = {cc: prepare(cc) for cc in compilers}
report.append(f"# Verification of {args.branch}\n")
report.append(f"Branch `{args.branch}`, commit {heads[compilers[0]]}; compilers {', '.join(compilers)}; "
              f"OPTFLAGS `{OPTFLAGS}`; {args.threads} threads; {time.strftime('%Y-%m-%d %H:%M')}.\n")
results["commit"] = heads[compilers[0]]

# build
if "build" in stages:
    report.append("## Build (fresh worktree)\n")
    report.append("| Compiler | make env | make cpu | make bench | Clean build | Warnings | Generated C |")
    report.append("|---|---|---|---|---:|---:|---|")
    for cc in compilers:
        env_status, _, env_s = run(["make", "env"], tree(cc), logs / f"{cc}-env.log")
        before = c_files(tree(cc))
        build = [make(cc, target, logs / f"{cc}-{target}.log") for target in ("cpu", "bench")]
        after = c_files(tree(cc))
        generated = sorted(p for p, digest in after.items() if before.get(p) != digest)
        lines = sum(len(p.read_text(errors="replace").splitlines()) for p in generated)
        size = sum(p.stat().st_size for p in generated)
        warnings = sum(output.count("warning:") for _, output, _ in build)
        seconds = sum(b[2] for b in build)
        report.append(f"| {cc} | exit {env_status} ({env_s:.0f} s) | exit {build[0][0]} ({build[0][2]:.1f} s) "
                      f"| exit {build[1][0]} ({build[1][2]:.1f} s) | {seconds:.1f} s | {warnings} "
                      f"| {len(generated)} files, {lines:,} lines, {size / 1e6:.2f} MB |")
        results.setdefault("build", {})[cc] = {
            "env_exit": env_status, "cpu_exit": build[0][0], "bench_exit": build[1][0], "clean_build_s": round(seconds, 1),
            "warnings": warnings, "generated_files": [str(p.relative_to(tree(cc))) for p in generated],
            "generated_lines": lines, "generated_bytes": size}
    report.append("")
    for cc in compilers:
        files = results["build"][cc]["generated_files"]
        report.append(f"{cc}: generated C: " + (", ".join(f"`{f}`" for f in files[:20]) or "none") +
                      (f", and {len(files) - 20} more" if len(files) > 20 else "") + "\n")

# tests
if "tests" in stages:
    cases = [line.split()[0] for line in (TIER2_ROOT / "reference" / "micm_cases.txt").read_text().splitlines()
             if line.strip() and not line.startswith("#")]
    case_set = set(cases)
    report.append("## MICM's test cases (the Tier 1 gate)\n")
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
        counts = {v: sum(1 for x in seen.values() for y in x if y == v) for v in ("PASS", "FAIL", "N/A")}
        missing = [c for c in cases if c not in seen]
        repeated = [c for c, v in seen.items() if len(v) > 1]
        unknown_ids = [c for c in seen if c not in case_set]
        failures[cc] = [c for c, v in seen.items() if "FAIL" in v]
        report.append(f"| {cc} | exit {status} ({seconds:.0f} s) | {counts['PASS']} | {counts['FAIL']} | {counts['N/A']} "
                      f"| {len(missing)} | {len(repeated)} | {len(unknown_ids)} |")
        results.setdefault("tests", {})[cc] = {"exit": status, **counts, "missing": len(missing),
                                               "repeated": len(repeated), "unknown": len(unknown_ids)}
    report.append("")
    for cc, failed in failures.items():
        if failed:
            report.append(f"{cc} failures: {' '.join(failed[:40])}\n")

# the Tier 1 gate's input, shared by ts1 and specialized
gate_dir = work.parent / "tier1-gate"
cells_2026 = gate_dir / "cells_2026.csv"
micm_2026 = gate_dir / "micm_final_2026.csv"
if "ts1" in stages or "specialized" in stages:
    gate_dir.mkdir(exist_ok=True)
    if not cells_2026.exists():
        run([sys.executable, str(TIER2_ROOT / "reference" / "make_cells_csv.py"), "10000", "2026", "--output",
             str(cells_2026)], gate_dir, gate_dir / "make_cells.log")
    if not micm_2026.exists():
        driver = gate_dir / "micm_ts1_ref"
        if not driver.exists():
            run(["clang++", "-std=c++20", *OPTFLAGS.split(), "-I", str(TIER2_ROOT / "micm/include"),
                 "-I", str(TIER2_ROOT / "micm/test/benchmark"), str(TIER2_ROOT / "reference/micm_ts1_ref.cpp"),
                 "-o", str(driver)], gate_dir, gate_dir / "build_micm_ts1_ref.log")
        status, _, _ = run([str(driver), "10000", "100", "30", str(cells_2026)], gate_dir, gate_dir / "micm_ts1_ref.log")
        if status == 0:
            (gate_dir / "ts1_final.csv").rename(micm_2026)

# ts1
if "ts1" in stages:
    report.append("## Tier 1 TS1 gate: bin/ts1_bench_cpu (10,000 cells, seed 2026, 100 steps of 30 s)\n")
    report.append("| Compiler | Exit | elapsed_ms | compare_final.py (1e-10) |")
    report.append("|---|---|---:|---|")
    for cc in compilers:
        out = work / "runs" / cc / "tier1-ts1"
        out.mkdir(parents=True, exist_ok=True)
        (out / "ts1_final.csv").unlink(missing_ok=True)
        status, output, _ = run([str(tree(cc) / "bin/ts1_bench_cpu"), "10000", "100", "30", str(cells_2026)], out,
                                logs / f"{cc}-ts1_bench_cpu.log")
        elapsed, _ = elapsed_and_stats(output)
        verdict, largest, text = ("no result", "", "")
        if status == 0 and (out / "ts1_final.csv").exists():
            verdict, largest, text = compare(out / "ts1_final.csv", micm_2026, logs / f"{cc}-ts1-compare.log")
        report.append(f"| {cc} | {status} | {elapsed} | {verdict}: {text} |")
        results.setdefault("ts1", {})[cc] = {"exit": status, "elapsed_ms": elapsed, "gate": verdict, "largest": largest}
    report.append("")

# specialized
if "specialized" in stages:
    inputs = [("seed 2026", cells_2026, micm_2026)] + [(f"input {i + 1}", Path(c).resolve(), Path(m).resolve())
                                                       for i, (c, m) in enumerate(args.ts1_input)]
    report.append("## Specialized TS1 gate: bin/ts1_bench_specialized (10,000 cells, 100 steps of 30 s)\n")
    report.append("| Compiler | Input | Exit | elapsed_ms | Internal steps | compare_final.py (1e-10) |")
    report.append("|---|---|---|---:|---:|---|")
    for cc in compilers:
        for label, cells, reference in inputs:
            slug = re.sub(r"\W+", "-", label)
            out = work / "runs" / cc / f"specialized-{slug}"
            out.mkdir(parents=True, exist_ok=True)
            (out / "ts1_final.csv").unlink(missing_ok=True)
            status, output, _ = run([str(tree(cc) / "bin/ts1_bench_specialized"), "10000", "100", "30", str(cells)], out,
                                    logs / f"{cc}-ts1_bench_specialized-{slug}.log")
            elapsed, values = elapsed_and_stats(output)
            verdict, largest, text = ("no result", "", "")
            if status == 0 and (out / "ts1_final.csv").exists():
                verdict, largest, text = compare(out / "ts1_final.csv", reference, logs / f"{cc}-specialized-{slug}-compare.log")
            report.append(f"| {cc} | {label} | {status} | {elapsed} | {values.get('internal_steps', '')} | {verdict}: {text} |")
            results.setdefault("specialized", {}).setdefault(cc, {})[label] = {
                "exit": status, "elapsed_ms": elapsed, "internal_steps": values.get("internal_steps", ""),
                "gate": verdict, "largest": largest}
    report.append("")

# generator
if "generator" in stages:
    report.append("## The generator's tests: make test-specialized\n")
    report.append("| Compiler | Exit | Time | PASS lines | FAIL lines | Last lines |")
    report.append("|---|---|---:|---:|---:|---|")
    for cc in compilers:
        status, output, seconds = make(cc, "test-specialized", logs / f"{cc}-test-specialized.log")
        passes = len(re.findall(r"\bPASS(ED)?\b", output))
        fails = len(re.findall(r"\bFAIL(ED|URE)?\b", output))
        tail = " / ".join(line.strip() for line in output.strip().splitlines()[-3:])[:200]
        report.append(f"| {cc} | {status} | {seconds:.0f} s | {passes} | {fails} | {tail} |")
        results.setdefault("generator", {})[cc] = {"exit": status, "seconds": round(seconds, 1), "pass_lines": passes,
                                                   "fail_lines": fails}
    report.append("")

# daynight
if "daynight" in stages:
    scenarios = sorted(d for base in args.scenarios for d in base.resolve().iterdir() if (d / "scenario.txt").exists())
    report.append(f"## Day-night scenarios (--threads {args.threads}, --rtol 1e-6 --atol 1e-15, 32-cell groups)\n")
    report.append("Agreement is score.py against MICM's result at the same settings, with `--level 10`; "
                  "accuracy is score.py against the reference solution, for scale.\n")
    report.append("| Compiler | Scenario | Exit | elapsed_ms | Internal steps (MICM's) | Agreement with MICM | "
                  "Worst species | Gate | Level 10 | Accuracy (MICM's) |")
    report.append("|---|---|---|---:|---|---:|---|---|---|---|")
    for cc in compilers:
        variant, extra = "", ["--group-size", "32"]
        for scenario in scenarios:
            out = work / "runs" / cc / "daynight" / scenario.name
            out.mkdir(parents=True, exist_ok=True)
            for stale in ("ts1_final.csv", "ts1_snapshots.csv"):
                (out / stale).unlink(missing_ok=True)
            command = [str(tree(cc, variant) / "bin/ts1_scenario_specialized"), "--threads", str(args.threads), *DAYNIGHT,
                       *extra, str(scenario / "scenario.txt")]
            status, output, _ = run(command, out, out / "run.log")
            if status != 0 and extra and not variant:
                # the program may fix its group size at build time: make GROUP_SIZE=32 bench, as the prompt allows
                prepare(cc, "-g32")
                run(["make", "env"], tree(cc, "-g32"), logs / f"{cc}-g32-env.log")
                built = make(cc, "bench", logs / f"{cc}-g32-bench.log", "-g32", ["GROUP_SIZE=32"])
                report_note = f"{cc}: bin/ts1_scenario_specialized rejected --group-size 32 (exit {status}); " \
                              f"rebuilt with `make GROUP_SIZE=32 bench` (exit {built[0]})"
                results.setdefault("daynight_build", {})[cc] = report_note
                variant, extra = "-g32", []
                command = [str(tree(cc, variant) / "bin/ts1_scenario_specialized"), "--threads", str(args.threads),
                           *DAYNIGHT, str(scenario / "scenario.txt")]
                status, output, _ = run(command, out, out / "run.log")
            elapsed, values = elapsed_and_stats(output)
            micm_log = (scenario / "micm.log").read_text() if (scenario / "micm.log").exists() else ""
            micm_steps = re.search(r"^gate: .*?internal_steps=(\d+)", micm_log, re.M)
            micm_accuracy = re.search(r"^gate run against tighter: overall ([0-9.-]+)", micm_log, re.M)
            snapshots = (scenario / "micm_snapshots.csv").exists()
            result = out / ("ts1_snapshots.csv" if snapshots else "ts1_final.csv")
            micm = scenario / ("micm_snapshots.csv" if snapshots else "micm_final.csv")
            reference = scenario / ("reference_snapshots.csv" if snapshots else "reference_final.csv")
            row = {"exit": status, "elapsed_ms": elapsed, "internal_steps": values.get("internal_steps", ""),
                   "micm_internal_steps": micm_steps.group(1) if micm_steps else "",
                   "agreement": "", "worst": "", "gate": "FAIL", "level10": "FAIL", "accuracy": "",
                   "micm_accuracy": micm_accuracy.group(1) if micm_accuracy else ""}
            if status == 0 and result.exists():
                _, scored, _ = run([sys.executable, str(TIER2_ROOT / "scenarios/score.py"), str(result), str(micm),
                                    "--level", "10", "--json"], out, out / "score-micm.log")
                try:
                    summary = json.loads(scored.strip().splitlines()[-1])
                except (IndexError, ValueError):
                    summary = {"error": scored.strip()[-200:]}
                row.update(agreement=summary.get("overall_digits", ""), gate=summary.get("gate", "FAIL"),
                           level10=summary.get("level_result", "FAIL"),
                           worst=" ".join(f"{n} {v}" for n, v in summary.get("worst_species", [])[:1]))
                if reference.exists():
                    _, scored, _ = run([sys.executable, str(TIER2_ROOT / "scenarios/score.py"), str(result), str(reference),
                                        "--json"], out, out / "score-reference.log")
                    try:
                        row["accuracy"] = json.loads(scored.strip().splitlines()[-1]).get("overall_digits", "")
                    except (IndexError, ValueError):
                        pass
            results.setdefault("daynight", {}).setdefault(cc, {})[scenario.name] = row
            report.append(f"| {cc} | {scenario.name} | {status} | {elapsed} | {row['internal_steps']} "
                          f"({row['micm_internal_steps']}) | {row['agreement']} | {row['worst']} | {row['gate']} "
                          f"| {row['level10']} | {row['accuracy']} ({row['micm_accuracy']}) |")
            print(f"{cc} {scenario.name}: exit {status} agreement {row['agreement']} level10 {row['level10']}", flush=True)
    report.append("")
    for cc, note in results.get("daynight_build", {}).items():
        report.append(f"{note}\n")
    for cc in compilers:
        rows = results.get("daynight", {}).get(cc, {})
        if rows:
            ok = all(r["level10"] == "PASS" for r in rows.values())
            report.append(f"- {cc}: level 10 and the gate on every scenario: {'yes' if ok else 'NO'} ({len(rows)} scenarios)")
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
    sources = sorted(c_files(tree(cc)))
    banned = re.compile(r"__attribute__|__builtin_|\basm\b|__asm|vector_size|ext_vector_type|long\s+double|"
                        r"#\s*include\s*<(arm_|immintrin|xmmintrin|emmintrin|x86intrin|Accelerate|dispatch|mach/|dlfcn)|"
                        r"#\s*pragma\s+(clang|GCC)|omp\s+target|__APPLE__|__aarch64__|__x86_64__|__ARM_NEON|"
                        r"_Float16|__fp16|__bf16|\bsystem\s*\(|\bpopen\s*\(|\bexec[lv]p?e?\s*\(|posix_spawn|\bfork\s*\(|"
                        r"\bdlopen\s*\(")
    matches = []
    for source in sources:
        for number, line in enumerate(source.read_text(errors="replace").splitlines(), 1):
            if banned.search(line):
                matches.append(f"{source.relative_to(tree(cc))}:{number}: {line.strip()[:140]}")
    (logs / "banned.txt").write_text("\n".join(matches) + "\n")
    report.append(f"C sources scanned: {len(sources)} (tracked and generated, outside conda); lines matching a banned "
                  f"construct or starting a process: {len(matches)} (all in `logs/banned.txt`)\n")
    for match in matches[:30]:
        report.append(f"    {match}")
    report.append("")
    results["audit"] = {"flags_outside_rules": suspicious, "banned_matches": len(matches)}

results_path.write_text(json.dumps(results, indent=1, default=str) + "\n")
section = work / ("report.md" if set(stages) == set(STAGES) else f"report-{'-'.join(stages)}.md")
section.write_text("\n".join(report) + "\n")
print("\n".join(report))
