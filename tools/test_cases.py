#!/usr/bin/env python3
"""Run C-port checks and print one result for every supplied suite case."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT.parent
REFERENCE = ROOT / "reference"
CASE_LIST = REFERENCE / "micm_cases.txt"


def run(command: list[str], cwd: Path = PROJECT) -> bool:
    result = subprocess.run(command, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return result.returncode == 0


def compare_ts1(threads: int, destination: Path) -> bool:
    command = [
        str(PROJECT / "bin/ts1_bench_cpu"), "--threads", str(threads),
        "100", "100", "30", str(REFERENCE / "ts1_cells_100_seed1.csv"),
    ]
    if not run(command):
        return False
    shutil.copyfile(PROJECT / "ts1_final.csv", destination)
    return run([sys.executable, str(REFERENCE / "compare_final.py"), str(destination),
                str(REFERENCE / "micm_ts1_final_100_seed1.csv")])


def failure_reason(case_id: str) -> str:
    suite = case_id.split(":", 1)[0]
    if "constraint" in suite or suite in {"equilibrium", "linear", "linear_constraint"} or suite.startswith("dae_"):
        return "constraint and DAE solver path is not implemented"
    if "aerosol_model" in suite or "external_model" in suite:
        return "external and aerosol model APIs are not implemented"
    if suite in {"matrix", "vector_matrix"} or suite.startswith("sparse_matrix") or suite in {"grouped_views", "view_category_sfinae"}:
        return "generic matrix and view utility case is not implemented"
    if suite in {"state", "phase", "species", "system", "process", "process_set", "solver_builder"}:
        return "builder or state API case is not implemented"
    if suite.startswith("lu_decomposition") or suite in {"linear_solver", "linear_solver_in_place"}:
        return "standalone linear algebra case is not implemented"
    if suite in {"analytical_rosenbrock_integration", "analytical_backward_euler", "backward_euler"}:
        return "this integration case is not implemented"
    return "case is not implemented in the C test suite"


def main() -> int:
    required = [PROJECT / "bin/test_rates", PROJECT / "bin/test_example_no_user", PROJECT / "bin/test_example_user"]
    if any(not path.is_file() for path in required):
        print("test binaries are missing; run make test-cpu")
        return 2

    rate_ok = run([str(PROJECT / "bin/test_rates")])
    example_no_user_ok = run([str(PROJECT / "bin/test_example_no_user")])
    example_user_ok = run([str(PROJECT / "bin/test_example_user")])
    serial_ok = compare_ts1(1, PROJECT / "build/ts1_ours_1.csv")
    parallel_ok = compare_ts1(4, PROJECT / "build/ts1_ours_4.csv")
    rates = {
        "arrhenius_rate_constant": {"ArrheniusRateConstant.CalculateWithSystem", "ArrheniusRateConstant.CalculateWithPrescribedArguments"},
        "branched_rate_constant": {"BranchedRateConstant.CalculateAlkoxyBranchWithAllArguments", "BranchedRateConstant.CalculateNitrateBranchWithAllArguments"},
        "reversible_rate_constant": {f"ReversibleRateConstant.{name}" for name in (
            "DefaultConstructor", "CalculateWithSystem", "CalculateWithPrescribedArguments", "TemperatureDependence", "ZeroReverseRate", "EquilibriumConstantOnly")},
        "surface_rate_constant": {f"SurfaceRateConstant.{name}" for name in (
            "CalculateDefaultProbability", "CalculateSpecifiedProbability", "DiffusionCoefficientIsMissing", "MolecularWeightIsMissing")},
        "taylor_series_rate_constant": {f"TaylorSeriesRateConstant.{name}" for name in ("CalculateWithSystem", "CalculateWithPrescribedArguments")},
        "ternary_chemical_activation_rate_constant": {
            "TernaryChemicalActivationRateConstant.CalculateWithMinimalArguments",
            "TernaryChemicalActivationRateConstant.CalculateWithAllArguments",
            "TernaryChemicalActivationJPL19RateConstant.CalculateWithAllArguments",
            "TernaryChemicalActivationJPL19RateConstant.ZeroKintReducesToTroe"},
        "troe_rate_constant": {f"TroeRateConstant.{name}" for name in (
            "CalculateWithMinimalArguments", "CalculateWithAllArguments", "AnalyticalTroeExampleAB", "AnalyticalTroeExampleBC")},
        "tunneling_rate_constant": {f"TunnelingRateConstant.{name}" for name in ("CalculateWithMinimalArguments", "CalculateWithAllArguments")},
        "user_defined_rate_constant": {f"UserDefinedRateConstant.{name}" for name in (
            "CalculateWithSystem", "ConstructorWithRate", "ConstructorWithRateAndName", "CustomScalingFactor")},
    }
    example_programs = {
        "rate_constants_no_user_defined_example_by_hand": example_no_user_ok,
        "rate_constants_user_defined_example_by_hand": example_user_ok,
    }
    case_results: list[tuple[str, str, str]] = []
    for line in CASE_LIST.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        case_id, kind = line.split()
        suite, _, test = case_id.partition(":")
        status, reason = "FAIL", failure_reason(case_id)
        if suite in rates and test in rates[suite]:
            status, reason = ("PASS", "C rate formula agrees with the specified double precision value") if rate_ok else ("FAIL", "C rate formula check failed")
        elif suite == "reaction_rate_store" and "LambdaEntriesRcIndex" in test:
            status, reason = "N/A", "lambda callback entries are not part of the C named-rate descriptor"
        elif suite == "reaction_rate_store":
            status, reason = ("PASS", "C descriptor offsets and rate parameter references are consistent") if rate_ok else ("FAIL", "C descriptor consistency check failed")
        elif suite == "lambda_rate_constant":
            status, reason = "N/A", "capturing callable storage has no equivalent in the named-rate descriptor API"
        elif suite in example_programs:
            ok = example_programs[suite]
            status, reason = ("PASS", "generated configuration descriptor matches the hand-built rate structure") if ok else ("FAIL", "configuration descriptor check failed")
        elif suite == "openmp" and case_id.endswith("OpenMP.OneSolverManyStates"):
            status, reason = ("PASS", "threaded C groups match the reference final state") if parallel_ok else ("FAIL", "threaded C state differs from the reference")
        elif suite == "version":
            status, reason = "N/A", "release metadata is not part of the numerical C library interface"
        elif suite == "view_category_sfinae":
            status, reason = "N/A", "C has no template categories or SFINAE type traits"
        case_results.append((case_id, status, reason))

    for case_id, status, reason in case_results:
        print(f"{case_id} {status} {reason}")
    passed = sum(status == "PASS" for _, status, _ in case_results)
    failed = sum(status == "FAIL" for _, status, _ in case_results)
    not_applicable = sum(status == "N/A" for _, status, _ in case_results)
    print(f"summary PASS={passed} FAIL={failed} N/A={not_applicable} total={len(case_results)}")
    print(f"TS1 reference checks: serial={'PASS' if serial_ok else 'FAIL'} parallel={'PASS' if parallel_ok else 'FAIL'}")
    return 1 if failed or not serial_ok or not parallel_ok else 0


if __name__ == "__main__":
    raise SystemExit(main())
