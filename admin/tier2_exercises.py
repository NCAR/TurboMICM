#!/usr/bin/env python3
"""Run two of the Tier 2 review's exercises on an entry: change TS1 and regenerate.

    admin/tier2_exercises.py RUN BRANCH WORK [--config PATH]

The prompt told the agents that reviewers would "add a reaction to TS1 and
regenerate, change a rate parameter, and run the generator on another
mechanism". Every entry builds both its generic TS1 benchmark and its
specialized one from one TS1 configuration file, so the first two exercises can
be the same for all of them. In a fresh worktree of BRANCH at WORK/RUN (with
links micm, reference, and scenarios beside it), with clang and OPTFLAGS
"-O3 -DNDEBUG":
  baseline  make env; make cpu bench; run bin/ts1_bench_cpu and
            bin/ts1_bench_specialized on reference/ts1_cells_100_seed1.csv
            (one thread, 100 steps of 30 s)
  add       append to the TS1 configuration a copy of its NO + O3 Arrhenius
            reaction, a second channel of the same rate; make bench; run both
  change    restore the configuration, then multiply that reaction's A by 1.5;
            make bench; run both
For each edit it records whether make regenerated and rebuilt (the build's
time and the specialized binary's change), and two checks: the specialized
result agrees with the generic one (reference/compare_final.py, 1e-10), and
both differ from the baseline (the same comparison fails). The configuration is
found from the Makefile rule that writes the specialized TS1 source (--config
names it instead). WORK/exercises.json and WORK/exercises.md hold the results.
Needs only the Python standard library.
"""

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
TIER2_ROOT = Path.home() / "NSF-NCAR/ODE/benchmarks/tier2"
OPTFLAGS = "-O3 -DNDEBUG"

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("run")
parser.add_argument("branch")
parser.add_argument("work", type=Path)
parser.add_argument("--config", help="the TS1 configuration, relative to the entry")
args = parser.parse_args()
work = args.work.resolve()
work.mkdir(parents=True, exist_ok=True)
tree = work / args.run
logs = work / "logs"
logs.mkdir(exist_ok=True)
for name in ("micm", "reference", "scenarios"):
    link = work / name
    if not link.is_symlink():
        os.symlink((TIER2_ROOT / name).resolve(), link)
if not tree.exists():
    subprocess.run(["git", "-C", str(REPO), "worktree", "add", "--detach", str(tree), args.branch],
                   check=True, capture_output=True)
cells = TIER2_ROOT / "reference" / "ts1_cells_100_seed1.csv"
compare = TIER2_ROOT / "reference" / "compare_final.py"


def run(command, log, cwd=tree):
    start = time.perf_counter()
    done = subprocess.run(command, cwd=cwd, capture_output=True, text=True, shell=isinstance(command, str))
    seconds = time.perf_counter() - start
    (logs / log).write_text(f"$ {command}\n# exit {done.returncode}; {seconds:.1f} s\n{done.stdout}{done.stderr}")
    return done.returncode, done.stdout + done.stderr, seconds


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else ""


def find_config():
    if args.config:
        return args.config
    makefile = (tree / "Makefile").read_text()
    for line in makefile.splitlines():
        target = re.match(r"^(\S*(specialized|sp_)\S*ts1\S*\.c|\S*ts1\S*specialized\S*\.c)\s*:(.*)$", line)
        if target:
            for word in target.group(3).split():
                if re.search(r"ts1\S*\.(json|ya?ml)$", word):
                    return word
    for candidate in ("configs/ts1.json", "mechanisms/ts1.json"):
        if (tree / candidate).exists():
            return candidate
    sys.exit("cannot find the TS1 configuration; give --config")


def reactions_of(document):
    """The list of reactions in a CAMP configuration document (camp-data MECHANISM records)."""
    data = document.get("camp-data", document) if isinstance(document, dict) else document
    for record in data:
        if isinstance(record, dict) and record.get("type") == "MECHANISM":
            return record["reactions"]
    raise SystemExit("no MECHANISM record in the TS1 configuration")


def names(side):
    if isinstance(side, dict):
        return set(side)
    return {item["species"] if isinstance(item, dict) else item for item in side}


def no_o3(reactions):
    for index, reaction in enumerate(reactions):
        if reaction.get("type") == "ARRHENIUS" and names(reaction.get("reactants", [])) == {"NO", "O3"}:
            return index
    raise SystemExit("no NO + O3 Arrhenius reaction in the TS1 configuration")


def bench_pair(label):
    results = {}
    for program in ("ts1_bench_cpu", "ts1_bench_specialized"):
        out = work / "runs" / label / program
        out.mkdir(parents=True, exist_ok=True)
        (out / "ts1_final.csv").unlink(missing_ok=True)
        status, output, _ = run(f"{tree}/bin/{program} --threads 1 100 100 30 {cells}", f"{label}-{program}.log", cwd=out)
        results[program] = out / "ts1_final.csv" if status == 0 and (out / "ts1_final.csv").exists() else None
    return results


def agree(a, b, log):
    if a is None or b is None:
        return None
    status, output, _ = run(f"{sys.executable} {compare} {a} {b}", log, cwd=work)
    found = re.search(r"largest relative difference (\S+)", output)
    return {"pass": status == 0, "largest": found.group(1) if found else ""}


config = find_config()
config_path = tree / config
original = config_path.read_text()
summary = {"run": args.run, "branch": args.branch, "config": config}

status_env, _, _ = run("make env", "env.log")
status_build, _, seconds = run(f"make cpu bench CC=clang OPTFLAGS='{OPTFLAGS}'", "build.log")
summary["baseline"] = {"env": status_env, "build": status_build, "build_s": round(seconds, 1)}
baseline = bench_pair("baseline")
summary["baseline"]["agree"] = agree(baseline["ts1_bench_specialized"], baseline["ts1_bench_cpu"], "baseline-agree.log")

edits = {}
for label in ("add", "change"):
    document = json.loads(original)
    reactions = reactions_of(document)
    index = no_o3(reactions)
    if label == "add":
        reactions.append(json.loads(json.dumps(reactions[index])))
        description = f"appended a copy of reaction {index} (NO + O3 -> NO2 + O2), a second channel"
    else:
        reactions[index]["A"] = reactions[index]["A"] * 1.5
        description = f"multiplied A of reaction {index} (NO + O3 -> NO2 + O2) by 1.5"
    config_path.write_text(json.dumps(document, indent=1) + "\n")
    before = digest(tree / "bin/ts1_bench_specialized")
    status, output, seconds = run(f"make bench CC=clang OPTFLAGS='{OPTFLAGS}'", f"{label}-build.log")
    rebuilt = digest(tree / "bin/ts1_bench_specialized") != before
    pair = bench_pair(label)
    edits[label] = {
        "edit": description,
        "make_bench": status,
        "build_s": round(seconds, 1),
        "specialized_rebuilt": rebuilt,
        "specialized_vs_generic": agree(pair["ts1_bench_specialized"], pair["ts1_bench_cpu"], f"{label}-agree.log"),
        "generic_changed": agree(pair["ts1_bench_cpu"], baseline["ts1_bench_cpu"], f"{label}-generic-vs-baseline.log"),
        "specialized_changed": agree(pair["ts1_bench_specialized"], baseline["ts1_bench_specialized"],
                                     f"{label}-specialized-vs-baseline.log"),
    }
    config_path.write_text(original)
summary["edits"] = edits
run(f"make bench CC=clang OPTFLAGS='{OPTFLAGS}'", "restore-build.log")


def verdict(edit):
    sv, gc, sc = edit["specialized_vs_generic"], edit["generic_changed"], edit["specialized_changed"]
    ok = (edit["make_bench"] == 0 and edit["specialized_rebuilt"] and sv and sv["pass"]
          and gc and not gc["pass"] and sc and not sc["pass"])
    return "PASS" if ok else "FAIL"


for edit in edits.values():
    edit["verdict"] = verdict(edit)
(work / "exercises.json").write_text(json.dumps(summary, indent=1, default=str) + "\n")
lines = [f"# TS1 exercises: {args.run} ({args.branch})\n", f"TS1 configuration: `{config}`\n",
         f"Baseline: make env exit {status_env}; make cpu bench exit {status_build} ({seconds:.0f} s); "
         f"specialized against generic: {summary['baseline']['agree']}\n",
         "| Exercise | Edit | make bench | Rebuilt | Specialized vs generic (1e-10) | Generic changed | Specialized changed | Verdict |",
         "|---|---|---|---|---|---|---|---|"]
for label, edit in edits.items():
    sv, gc, sc = edit["specialized_vs_generic"], edit["generic_changed"], edit["specialized_changed"]
    lines.append(f"| {label} | {edit['edit']} | exit {edit['make_bench']} ({edit['build_s']} s) | "
                 f"{'yes' if edit['specialized_rebuilt'] else 'no'} | "
                 f"{('PASS ' + sv['largest']) if sv and sv['pass'] else ('FAIL ' + sv['largest']) if sv else 'no result'} | "
                 f"{('yes, ' + gc['largest']) if gc and not gc['pass'] else 'no' if gc else 'no result'} | "
                 f"{('yes, ' + sc['largest']) if sc and not sc['pass'] else 'no' if sc else 'no result'} | {edit['verdict']} |")
(work / "exercises.md").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
