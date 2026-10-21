#!/usr/bin/env python3
"""One counting rule for every entrant, at its final commit (tracked files only).

- C benchmark: the TS1 driver.
- C tests: tests/ and examples/, and ports of MICM's tutorial programs wherever they live.
- C library: the other C under src/ and include/.
- Python tests: tests/, and tools/ scripts named test_*, run_tests*, check_*, validate_*.
- Python tooling: the other Python.
Function sizes and comment share are over the C library and benchmark.
"""
import json, subprocess, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from funcs import functions, comment_line_share

BASE = Path.home() / "NSF-NCAR/tier1-work/verify"
RUNS = ["sol-max", "astra-xhigh", "sol-xhigh", "astra-max", "sol-high", "astra-high", "luna-xhigh", "luna-max", "luna-high"]
BENCH = {"sol-high": "src/bench.c", "sol-xhigh": "src/bench.c", "sol-max": "src/bench.c",
         "astra-high": "src/ts1_bench.c", "astra-xhigh": "src/bench.c", "astra-max": "tools/ts1_bench.c",
         "luna-max": "src/benchmark.c", "luna-xhigh": "src/ts1_bench.c", "luna-high": "src/main.c"}
TUTORIAL_PORTS = {"src/programs.c", "src/tutorial.c", "tools/tutorials.c"}
GEN = ["generated/ts1.c", "build/generated/ts1.c", "build/ts1.c", "build/ts1_model.c"]

def cat(f, run):
    name = f.rsplit("/", 1)[-1]
    if f == BENCH[run]:
        return "c_bench"
    if f.endswith((".c", ".h")):
        if f.startswith(("tests/", "examples/")) or f in TUTORIAL_PORTS:
            return "c_tests"
        if f.startswith(("src/", "include/")):
            return "c_lib"
        return "c_other"
    if f.endswith(".py"):
        if f.startswith("tests/") or name.startswith(("test_", "run_tests", "check_", "validate_")):
            return "py_tests"
        return "py_tools"
    if f.endswith(".md"):
        return "docs"
    return "other"

out = {}
for run in RUNS:
    d = BASE / run / "clang" / run
    files = subprocess.run(["git", "-C", str(d), "ls-files"], capture_output=True, text=True, check=True).stdout.split()
    lines, nfiles, fs = {}, {}, []
    code = comm = 0
    for f in files:
        c = cat(f, run)
        text = (d / f).read_text(errors="replace")
        n = text.count("\n")
        lines[c] = lines.get(c, 0) + n
        nfiles[c] = nfiles.get(c, 0) + 1
        if c in ("c_lib", "c_bench"):
            if f.endswith(".c"):
                fs += [(l, f"{f}:{s}", nm) for nm, s, l in functions(text)]
            a, b = comment_line_share(text)
            comm += a
            code += b
    gen = {g: (d / g).read_text(errors="replace").count("\n") for g in GEN if (d / g).exists()}
    fs.sort(reverse=True)
    lens = sorted(l for l, _, _ in fs)
    out[run] = {"lines": lines, "files": nfiles, "generated": gen,
                "functions": len(fs), "median_function": lens[len(lens) // 2],
                "over_60": sum(l > 60 for l in lens), "over_100": sum(l > 100 for l in lens),
                "largest": [f"{nm} ({loc}, {l} lines)" for l, loc, nm in fs[:6]],
                "lines_with_comments": round(comm / code, 3)}
json.dump(out, open(Path(__file__).resolve().parent.parent / "metrics.json", "w"), indent=1)
for run, v in out.items():
    L = v["lines"]
    print(f"{run:12s} lib {L.get('c_lib',0):5d} ({v['files'].get('c_lib',0)} files) bench {L.get('c_bench',0):4d} | ctests {L.get('c_tests',0):5d} pytests {L.get('py_tests',0):5d} | pytools {L.get('py_tools',0):4d} | cother {L.get('c_other',0)} | gen {list(v['generated'].values())} | funcs {v['functions']} median {v['median_function']} >60 {v['over_60']} >100 {v['over_100']} comments {v['lines_with_comments']}")
