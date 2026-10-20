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
    if suite in {
        "constraint_set", "equilibrium", "linear", "constraint_initialization", "equilibrium_constraint",
        "linear_constraint", "dae_constraint_overshoot", "dae_algebraic_error_insensitivity",
    }:
        return "constraint solver case is not implemented"
    if "aerosol_model" in suite or "external_model" in suite:
        return "external and aerosol model APIs are not implemented"
    if suite in {"matrix", "vector_matrix"} or suite.startswith("sparse_matrix") or suite in {"grouped_views", "view_category_sfinae"}:
        return "generic matrix and view utility case is not implemented"
    if suite in {"state", "phase", "species", "system", "process", "process_set", "solver_builder"}:
        return "builder or state API case is not implemented"
    if suite.startswith("lu_decomposition"):
        return "dense LU case is not implemented"
    if suite in {"analytical_rosenbrock_integration", "analytical_backward_euler", "backward_euler"}:
        return "this integration case is not implemented"
    return "case is not implemented in the C test suite"


def main() -> int:
    required = [
        PROJECT / "bin/test_rates", PROJECT / "bin/test_example_no_user",
        PROJECT / "bin/test_example_user", PROJECT / "bin/test_tutorial_programs",
        PROJECT / "bin/test_jacobian_verification", PROJECT / "bin/test_solver_configuration",
        PROJECT / "bin/test_integration_equivalents", PROJECT / "bin/test_state", PROJECT / "bin/test_chapman",
        PROJECT / "bin/test_objects", PROJECT / "bin/test_processes", PROJECT / "bin/test_constraints",
        PROJECT / "bin/test_lu",
    ]
    if any(not path.is_file() for path in required):
        print("test binaries are missing; run make test-cpu")
        return 2

    rate_ok = run([str(PROJECT / "bin/test_rates")])
    example_no_user_ok = run([str(PROJECT / "bin/test_example_no_user")])
    example_user_ok = run([str(PROJECT / "bin/test_example_user")])
    tutorial_programs = {
        name: run([str(PROJECT / "bin/test_tutorial_programs"), name])
        for name in ("README_example", "multiple_grid_cells", "openmp_tutorial", "vectorized_matrix_solver")
    }
    jacobian_case_names = {
        "FiniteDifferenceMatchesAnalytical", "MultiBlockFiniteDifference",
        "ComparePassesForCorrectJacobian", "CompareFailsForWrongJacobian",
        "SparsityCompletenessPassesWhenComplete", "SparsityCompletenessFailsWhenMissingEntry",
        "NearZeroVariableHandled",
    }
    jacobian_results = {
        name: run([str(PROJECT / "bin/test_jacobian_verification"), name])
        for name in jacobian_case_names
    }
    integration_results = {
        name: run([str(PROJECT / "bin/test_integration_equivalents"), name])
        for name in (
            "terminator", "chapman_integration", "integrated_reaction_rates", "runtime_vector_use",
            "lambda_rate_constant", "solver_lambda_rate_constant",
        )
    }
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
        elif suite in {"lambda_rate_constant", "solver_lambda_rate_constant"}:
            ok = integration_results[suite]
            status, reason = ("PASS", "C function-pointer rate matches the condition-dependent lambda") if ok else (
                "FAIL", "C function-pointer rate evaluation or override failed")
        elif suite in example_programs:
            ok = example_programs[suite]
            status, reason = ("PASS", "generated configuration descriptor matches the hand-built rate structure") if ok else ("FAIL", "configuration descriptor check failed")
        elif suite == "solve_results":
            ok = run([str(PROJECT / "bin/test_tutorial_programs"), suite])
            status, reason = ("PASS", "C solver statistics match the MICM tutorial") if ok else (
                "FAIL", "C solve_results statistics differ from MICM")
        elif suite == "jacobian_verification" and test.split(".", 1)[-1] in jacobian_results:
            ok = jacobian_results[test.split(".", 1)[-1]]
            status, reason = ("PASS", "C finite-difference and Jacobian checks match the analytical quantities") if ok else (
                "FAIL", "C finite-difference or Jacobian comparison check failed")
        elif suite == "solver_configuration":
            ok = run([str(PROJECT / "bin/test_solver_configuration")])
            status, reason = ("PASS", "five C Rosenbrock parameter sets match the MICM tutorial") if ok else (
                "FAIL", "C Rosenbrock tutorial output or solver statistics differ from MICM")
        elif suite == "state":
            ok = run([str(PROJECT / "bin/test_state"), test])
            status, reason = ("PASS", "C state storage and named access match the reference behavior") if ok else (
                "FAIL", "C state storage, ownership, or validation differs from the reference behavior")
        elif suite in {"species", "phase", "system"}:
            ok = run([str(PROJECT / "bin/test_objects"), test])
            status, reason = ("PASS", "C owning species, phase, and system objects match the reference behavior") if ok else (
                "FAIL", "C species metadata, phase membership, or system naming differs from the reference behavior")
        elif suite in {"process", "process_set"}:
            ok = run([str(PROJECT / "bin/test_processes"), test])
            status, reason = ("PASS", "C reaction descriptors, forcing, and Jacobians match the reference behavior") if ok else (
                "FAIL", "C reaction construction, accumulation, masking, or Jacobian behavior differs from the reference")
        elif suite in {
            "constraint_set", "equilibrium", "linear", "constraint_initialization", "equilibrium_constraint",
            "linear_constraint", "dae_constraint_overshoot", "dae_algebraic_error_insensitivity",
        }:
            ok = run([str(PROJECT / "bin/test_constraints"), test])
            status, reason = ("PASS", "C constraint construction, evaluation, initialization, and integration match the reference behavior") if ok else (
                "FAIL", "C constraint validation, projection, residual, row mapping, Jacobian, or integration differs from the reference")
        elif suite.startswith("lu_decomposition_") or suite in {"linear_solver", "linear_solver_in_place"}:
            ok = run([str(PROJECT / "bin/test_lu"), test])
            status, reason = ("PASS", "C LU factorization and solve preserve the reference matrix equations") if ok else (
                "FAIL", "C LU factors, triangular structure, solve, or matrix reconstruction differ from the reference")
        elif suite in {"regression_test_rosenbrock_p_force", "regression_test_rosenbrock_dforce_dy"}:
            ok = run([str(PROJECT / "bin/test_chapman"), test])
            status, reason = ("PASS", "C Chapman rates, multi-cell forcing, and Jacobian match the regression equations") if ok else (
                "FAIL", "C Chapman rate, forcing, or Jacobian regression differs from the reference")
        elif suite == "chapman_ode_solver" and test.split(".", 1)[-1] in {
            "DefaultConstructor", "ReactionNames", "PhotolysisNames", "SpeciesNames",
            "simple_force", "smaller_force", "Solve", "solve_10_times_larger",
            "solve_10_times_smaller", "solve_with_random_number_densities",
        }:
            ok = run([str(PROJECT / "bin/test_chapman"), test])
            status, reason = ("PASS", "generic C mechanism and solver match Chapman reference values") if ok else (
                "FAIL", "generic C Chapman forcing or integration differs from the reference behavior")
        elif suite in integration_results:
            ok = integration_results[suite]
            status, reason = ("PASS", "C integration preserves the reference analytic or chemical invariant") if ok else (
                "FAIL", "C integration result or invariant differs from the reference")
        elif suite in tutorial_programs:
            ok = tutorial_programs[suite]
            status, reason = ("PASS", "C example matches the MICM reference result or storage behavior") if ok else ("FAIL", "C tutorial example differs from the MICM reference")
        elif suite == "openmp" and case_id.endswith("OpenMP.OneSolverManyStates"):
            status, reason = ("PASS", "threaded C groups match the reference final state") if parallel_ok else ("FAIL", "threaded C state differs from the reference")
        elif suite == "version":
            status, reason = "N/A", "release metadata is not part of the numerical C library interface"
        elif suite == "view_category_sfinae":
            if test.endswith("RuntimeVectorUse"):
                ok = integration_results["runtime_vector_use"]
                status, reason = ("PASS", "C array iteration produces the same vector sum") if ok else (
                    "FAIL", "C array iteration produced an incorrect sum")
            else:
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
