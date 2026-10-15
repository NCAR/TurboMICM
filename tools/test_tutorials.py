"""Compare complete standalone tutorial outputs with captured MICM outputs.

Only measured duration and concurrent thread-announcement interleaving vary
between runs. Validate their formats and contents before normalizing them.
All other bytes, including concentrations, statistics and spacing, must match.
"""
from pathlib import Path
import difflib
import math
import re
import subprocess
import sys


def normalize(name, output):
    if name == "solver_configuration":
        elapsed = re.findall(r"^Total solve time: (\S+) nanoseconds$", output, re.M)
        if len(elapsed) != 5 or any(not math.isfinite(float(t)) or float(t) < 0 for t in elapsed):
            raise ValueError("invalid elapsed-time report")
        output = re.sub(r"^Total solve time: \S+ nanoseconds$",
                        "Total solve time: <elapsed> nanoseconds", output, flags=re.M)
    if name == "openmp_tutorial":
        prefix, marker, rest = output.partition("Thread 1\n")
        phrase = "Running solver on thread "
        if not marker or prefix.count(phrase) != 3:
            raise ValueError("missing thread announcements")
        if sorted(prefix.replace(phrase, "").split()) != ["0", "1", "2"]:
            raise ValueError("wrong thread IDs")
        output = "<threads 0, 1, 2>\n" + marker + rest
    return output


def run_tutorials(root):
    results = {}
    for line in (root / "tests/cases.txt").read_text().splitlines():
        if not line.endswith(" program"):
            continue
        name = line.split()[0]
        try:
            reference = (root / "tests/data/tutorials" / (name + ".txt")).read_text()
            run = subprocess.run([root / "bin/tutorial_cpu", name], capture_output=True,
                                 text=True, timeout=120, cwd=root)
            if run.returncode or run.stderr:
                raise ValueError(f"exit={run.returncode}, stderr={run.stderr}")
            expected, actual = normalize(name, reference), normalize(name, run.stdout)
            if actual != expected:
                diff = "".join(difflib.unified_diff(expected.splitlines(True), actual.splitlines(True),
                                                    fromfile="MICM", tofile="C port"))
                sys.stderr.write(name + "\n" + diff[:8000])
                raise ValueError("stdout differs from MICM")
            results[name] = ("PASS", "reference stdout matched")
        except (OSError, ValueError, subprocess.TimeoutExpired) as error:
            sys.stderr.write(f"{name}: {error}\n")
            results[name] = ("FAIL", "tutorial output comparison failed")
    return results


if __name__ == "__main__":
    results = run_tutorials(Path(__file__).resolve().parents[1])
    for name, (status, reason) in results.items():
        print(name, status, reason)
    sys.exit(any(status != "PASS" for status, _ in results.values()))
