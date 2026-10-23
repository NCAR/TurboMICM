#!/usr/bin/env python3
"""Exercise C orchestration for sol-xhigh (the preferred route: the benchmark programs).

    run_mech.py MECH [MECH ...]

For each mechanism: put converted/MECH.json in place of the entry's TS1
configuration (mechanisms/ts1.json), time the generator (tools/specialize.py
with the Makefile's arguments, 3 runs, into gen/MECH/), run and time
`make bench CC=clang OPTFLAGS="-O3 -DNDEBUG"`, measure the generated C, run
bin/ts1_bench_cpu and bin/ts1_bench_specialized (--threads 1 4 60 60 CELLS)
each in its own empty directory under runs/MECH/, and compare the final
states with the experiment's compare_final.py. Results go to results/MECH.json
and logs to logs/. The original TS1 configuration is always restored at the
end (from backup/ts1.json, with a fresh modification time so the next build
regenerates TS1). Standard library only; writes only inside this directory.
"""
import json
import re
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

WORK = Path(__file__).resolve().parent
TREE = WORK / "sol-xhigh"
CONFIG = TREE / "mechanisms" / "ts1.json"
BACKUP = WORK / "backup" / "ts1.json"
CANON = WORK.parent.parent / "canonical"
COMPARE = Path.home() / "NSF-NCAR/ODE/benchmarks/tier2/reference/compare_final.py"
PYTHON = str(TREE / ".conda/bin/python")
# The Makefile's defaults (SPECIALIZE_* and GROUP_SIZE), as in the baseline build's command line.
SPECIALIZE_ARGS = ["--symbol", "ts1_specialized", "--batch", "8", "--split", "32", "--algebra-batch", "8",
                   "--width", "32", "--tile", "0", "--fuse", "0", "--joint-batch", "1", "--lu-order", "column"]
MAKE = ["make", "bench", "CC=clang", "OPTFLAGS=-O3 -DNDEBUG"]


def run(command, cwd, log=None):
    start = time.perf_counter()
    done = subprocess.run(command, cwd=cwd, capture_output=True, text=True)
    seconds = time.perf_counter() - start
    if log:
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(f"$ {' '.join(map(str, command))}\n# cwd {cwd}; exit {done.returncode}; {seconds:.2f} s\n"
                       f"--- stdout\n{done.stdout}--- stderr\n{done.stderr}")
    return done, seconds


def size(path):
    data = path.read_bytes()
    return {"lines": data.count(b"\n"), "bytes": len(data)}


def compare(a, b, log):
    done, _ = run([sys.executable, str(COMPARE), str(a), str(b)], WORK, log)
    text = (done.stdout + done.stderr).strip()
    match = re.search(r"largest relative difference (\S+)", text)
    largest = match.group(1) if match else "?"
    return {"verdict": f"{'PASS' if done.returncode == 0 else 'FAIL'} {largest}", "exit": done.returncode,
            "output": text}


def mechanism(name):
    logs, result = WORK / "logs" / name, {"mechanism": name}
    converted, cells = WORK / "converted" / f"{name}.json", CANON / name / "cells.csv"
    micm = CANON / name / "micm-run" / "ts1_final.csv"
    shutil.copyfile(converted, CONFIG)  # fresh mtime: make regenerates both descriptors
    result["converted"] = True

    gen_dir = WORK / "gen" / name
    gen_dir.mkdir(parents=True, exist_ok=True)
    times, errors = [], None
    for k in range(3):
        done, seconds = run([PYTHON, "tools/specialize.py", "mechanisms/ts1.json", str(gen_dir / f"sp_ts1_run{k}.c"),
                             *SPECIALIZE_ARGS], TREE, logs / f"generator-run{k}.log")
        if done.returncode:
            errors = done.stderr.strip()
            break
        times.append(seconds)
    done, seconds = run([PYTHON, "tools/configure.py", "mechanisms/ts1.json", str(gen_dir / "ts1.c"), "--symbol", "ts1"],
                        TREE, logs / "configure.log")
    result["configure_exit"], result["configure_s"] = done.returncode, round(seconds, 3)
    result["generator_ran"] = errors is None
    result["generator_times_s"] = [round(t, 3) for t in times]
    result["generator_s"] = round(statistics.median(times), 3) if times else None
    result["generator_error"] = errors

    done, seconds = run(MAKE, TREE, logs / "make-bench.log")
    result["build_exit"], result["build_s"] = done.returncode, round(seconds, 2)
    result["build_tail"] = "\n".join((done.stdout + done.stderr).strip().splitlines()[-15:])
    warnings = [line for line in (done.stdout + done.stderr).splitlines() if "warning:" in line]
    result["build_warnings"] = len(warnings)
    result["build_commands"] = [line.split(" -o ")[-1] if " -o " in line else line[:120]
                                for line in done.stdout.splitlines() if line.startswith(("clang", PYTHON))]

    generated = TREE / "generated"
    if (generated / "sp_ts1.c").exists():
        result["generated_specialized"] = size(generated / "sp_ts1.c")
        result["generated_c_lines"] = result["generated_specialized"]["lines"]
        result["generated_c_bytes"] = result["generated_specialized"]["bytes"]
        result["build_c_matches_timed_generator"] = (
            times != [] and (generated / "sp_ts1.c").read_bytes() == (gen_dir / "sp_ts1_run0.c").read_bytes())
        result["generator_deterministic"] = len({(gen_dir / f"sp_ts1_run{k}.c").read_bytes()
                                                 for k in range(len(times))}) == 1
    if (generated / "ts1.c").exists():
        result["generated_descriptors"] = size(generated / "ts1.c")

    finals = {}
    for label, program in (("generic", "ts1_bench_cpu"), ("specialized", "ts1_bench_specialized")):
        directory = WORK / "runs" / name / label
        if directory.exists():
            shutil.rmtree(directory)
        directory.mkdir(parents=True)
        binary = TREE / "bin" / program
        ok = result["build_exit"] == 0 and binary.exists()
        if ok:
            done, seconds = run([str(binary), "--threads", "1", "4", "60", "60", str(cells)], directory,
                                logs / f"{label}-run.log")
            ok = done.returncode == 0 and (directory / "ts1_final.csv").exists()
            result[f"{label}_output"] = (done.stdout + done.stderr).strip()
            result[f"{label}_exit"] = done.returncode
            result[f"{label}_wall_s"] = round(seconds, 3)
        result[f"{label}_ran"] = ok
        if ok:
            finals[label] = directory / "ts1_final.csv"

    for key, a, b in (("specialized_vs_generic", "specialized", "generic"),
                      ("generic_vs_micm", "generic", "micm"),
                      ("specialized_vs_micm", "specialized", "micm")):
        paths = {**finals, "micm": micm}
        if a in paths and b in paths:
            outcome = compare(paths[a], paths[b], logs / f"compare-{key}.log")
            result[key], result[key + "_output"] = outcome["verdict"], outcome["output"]
        else:
            result[key] = None
    (WORK / "results").mkdir(exist_ok=True)
    (WORK / "results" / f"{name}.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    try:
        for name in sys.argv[1:]:
            r = mechanism(name)
            print(f"{name}: generator {r['generator_s']} s, build exit {r['build_exit']} in {r['build_s']} s, "
                  f"C {r.get('generated_c_lines')} lines / {r.get('generated_c_bytes')} B; "
                  f"ran generic {r['generic_ran']} specialized {r['specialized_ran']}; "
                  f"S-G {r['specialized_vs_generic']}; G-M {r['generic_vs_micm']}; S-M {r['specialized_vs_micm']}",
                  flush=True)
    finally:
        shutil.copyfile(BACKUP, CONFIG)
        print("restored mechanisms/ts1.json from backup", flush=True)


if __name__ == "__main__":
    main()
