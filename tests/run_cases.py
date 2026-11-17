#!/usr/bin/env python3
"""Run ported cases and report every reference ID, including unported cases."""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
EXPECTED = ROOT / "tests" / "expected"
cases = [line.split()[0] for line in (ROOT / "tests" / "cases.txt").read_text().splitlines()
         if line.strip() and not line.startswith("#")]
results = {}
unexpected_failure = False

providers = [(name, [ROOT / "bin" / name]) for name in
             ("reference_rates", "reference_rate_store", "reference_solver", "reference_rosenbrock", "reference_constraints", "reference_constraint_set", "reference_equilibrium_integration", "reference_stiff_constraints", "reference_analytical", "reference_nonlinear_integration", "reference_state", "reference_matrix", "reference_matrix_functions", "reference_sparse", "reference_sparse_functions", "reference_grouped_views", "reference_jacobian", "reference_aerosol", "reference_external", "reference_version", "reference_initialization", "reference_lu")]
providers.append(("reference_ordering", [sys.executable, ROOT / "tests" / "reference_ordering.py"]))
for executable, command in providers:
    process = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    for line in process.stdout.splitlines():
        parts = line.split(" ", 2)
        if len(parts) != 3 or parts[0] not in cases or parts[0] in results or parts[1] not in ("PASS", "FAIL"):
            raise RuntimeError(f"invalid case report from {executable}: {line!r}")
        results[parts[0]] = (parts[1], parts[2])
    if process.returncode and not any(line.split(" ", 2)[1:2] == ["FAIL"] for line in process.stdout.splitlines()):
        unexpected_failure = True
        print(f"{executable}: {process.stderr}", file=sys.stderr)

for argument, case in (("no_user", "rate_constants_no_user_defined_example_by_hand"),
                       ("user", "rate_constants_user_defined_example_by_hand")):
    process = subprocess.run([ROOT / "bin" / "tutorial_program", argument], cwd=ROOT, capture_output=True, text=True)
    ok = process.returncode == 0 and process.stdout == (EXPECTED / f"{argument}.txt").read_text()
    results[case] = ("PASS", "printed reference output matched") if ok else ("FAIL", "printed output or program status differed")

for executable in ("test_core", "test_configurations", "test_constraints", "test_linear"):
    process = subprocess.run([ROOT / "bin" / executable], cwd=ROOT, capture_output=True, text=True)
    if process.returncode:
        unexpected_failure = True
        print(f"{executable}: {process.stderr}", file=sys.stderr)

for case in (
    "view_category_sfinae:ViewCategorySFINAE.VectorTypes",
    "view_category_sfinae:ViewCategorySFINAE.MatrixViews",
    "view_category_sfinae:ViewCategorySFINAE.SparseMatrixViews",
    "view_category_sfinae:ViewCategorySFINAE.MatrixTypesExcluded",
):
    results[case] = ("N/A", "C has no template concepts or SFINAE")

for case in cases:
    status, reason = results.get(case, ("FAIL", "case has not yet been ported"))
    print(case, status, reason)
    unexpected_failure |= status == "FAIL"
sys.exit(1 if unexpected_failure else 0)
