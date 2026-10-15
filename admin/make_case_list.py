#!/usr/bin/env python3
"""Write reference/micm_cases.txt, every MICM test case, from a built MICM.

    admin/make_case_list.py

ctest runs one executable per test. An executable linked with GoogleTest gives
one case per GoogleTest test, with id <ctest name>:<Suite.Test>; any other
executable is a plain program, with id <ctest name>.
"""

import argparse
import json
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--build", type=Path, default=ROOT / "micm" / "build")
parser.add_argument("--output", type=Path, default=ROOT / "reference" / "micm_cases.txt")
args = parser.parse_args()


def is_gtest(executable):
    symbols = subprocess.run(["nm", executable], capture_output=True, text=True, check=True).stdout
    return "testing8internal" in symbols  # mangled testing::internal


def gtest_cases(executable, cwd):
    listing = subprocess.run(
        [executable, "--gtest_list_tests"], capture_output=True, text=True, check=True, cwd=cwd
    ).stdout
    cases, suite = [], None
    for line in listing.splitlines():
        name = line.split("#")[0].strip()  # drop "# GetParam() = ..." and "# TypeParam = ..."
        if not name:
            continue
        if line.startswith(" "):
            cases.append(suite + name)
        else:
            suite = name  # "Suite." or "Prefix/Suite."
    return cases


tests = json.loads(
    subprocess.run(
        ["ctest", "--show-only=json-v1"], capture_output=True, text=True, check=True, cwd=args.build
    ).stdout
)["tests"]

lines, n_gtest, n_programs, n_executables = [], 0, 0, 0
for test in tests:
    executable = test["command"][0]
    cwd = next((p["value"] for p in test.get("properties", []) if p["name"] == "WORKING_DIRECTORY"), args.build)
    if is_gtest(executable):
        cases = gtest_cases(executable, cwd)
        lines += [f"{test['name']}:{case} gtest" for case in cases]
        n_gtest += len(cases)
        n_executables += 1
    else:
        lines.append(f"{test['name']} program")
        n_programs += 1

commit = subprocess.run(
    ["git", "-C", str(args.build.parent), "rev-parse", "--short", "HEAD"], capture_output=True, text=True, check=True
).stdout.strip()
header = [
    f"# MICM test cases at commit {commit}, written by admin/make_case_list.py.",
    f"# {n_gtest + n_programs} cases: {n_gtest} GoogleTest cases in {n_executables} executables"
    f" and {n_programs} plain programs.",
    "# Each line is <id> <kind>. The id is <ctest name>:<Suite.Test> for a GoogleTest case and",
    "# <ctest name> for a plain program; the kind is gtest or program.",
]
args.output.write_text("\n".join(header + lines) + "\n")
print(f"Wrote {args.output}: {n_gtest} GoogleTest cases in {n_executables} executables, {n_programs} plain programs")
