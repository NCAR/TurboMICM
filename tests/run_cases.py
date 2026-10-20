#!/usr/bin/env python3
"""Run the C-port cases and validate the generated configuration pathways."""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REFERENCE = ROOT.parent / "reference"


def load_generator():
    spec = importlib.util.spec_from_file_location("chem_generate_c", ROOT / "tools" / "generate_c.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def verify_configurations() -> None:
    generator = load_generator()
    cases = {
        "rate_constants_no_user_defined": (7, 0),
        "rate_constants_user_defined": (10, 3),
        "robertson": (3, 0),
    }
    for directory, (reaction_count, custom_minimum) in cases.items():
        extension = "yaml" if directory == "rate_constants_no_user_defined" else "json"
        config = ROOT / "examples" / directory / f"config.{extension}"
        model = generator.read_config(config)
        c_source = generator.compile_model(model)
        if len(model["reactions"]) != reaction_count:
            raise RuntimeError(f"{directory}: expected {reaction_count} parsed reactions")
        custom_labels = {reaction.get("label") for reaction in model["reactions"] if reaction.get("label")}
        if len(custom_labels) < custom_minimum:
            raise RuntimeError(f"{directory}: expected user-defined rate parameters")
        if "const chem_mechanism chem_model" not in c_source:
            raise RuntimeError(f"{directory}: no C mechanism table was emitted")


def read_case_ids():
    result = []
    for line in (REFERENCE / "micm_cases.txt").read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            result.append(line.split()[0])
    return result


def verify_ts1_gate() -> tuple[bool, str]:
    binary = ROOT / "bin" / "ts1_bench_cpu"
    run = subprocess.run(
        [str(binary), "--threads", "1", "100", "100", "30", str(REFERENCE / "ts1_cells_100_seed1.csv")],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if run.returncode:
        return False, "TS1 benchmark failed: " + run.stderr.strip()
    compare = subprocess.run(
        [sys.executable, str(REFERENCE / "compare_final.py"), str(ROOT / "ts1_final.csv"),
         str(REFERENCE / "micm_ts1_final_100_seed1.csv")],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return compare.returncode == 0, compare.stdout.strip() or compare.stderr.strip()


def verify_example_runtime() -> dict[str, tuple[bool, str]]:
    programs = {
        "rate_constants_no_user_defined_example_by_hand": ROOT / "bin" / "test_example_no_user",
        "rate_constants_user_defined_example_by_hand": ROOT / "bin" / "test_example_user",
    }
    results = {}
    for case_id, binary in programs.items():
        run = subprocess.run(
            [str(binary)], cwd=ROOT, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        message = run.stdout.strip() or run.stderr.strip()
        results[case_id] = (run.returncode == 0, message)
    return results


def main() -> int:
    verify_configurations()
    run = subprocess.run(
        [str(ROOT / "bin" / "test_cpu")], cwd=ROOT, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    if run.returncode not in (0, 1):
        print(run.stderr, file=sys.stderr)
        return run.returncode
    rows = []
    for line in run.stdout.splitlines():
        fields = line.split(maxsplit=2)
        if len(fields) != 3:
            print(f"invalid case output: {line}", file=sys.stderr)
            return 2
        rows.append([fields[0], fields[1], fields[2]])
    expected = read_case_ids()
    if [row[0] for row in rows] != expected:
        print("C test runner output does not match the MICM case list", file=sys.stderr)
        return 2

    gate_ok, gate_message = verify_ts1_gate()
    example_results = verify_example_runtime()
    regression_id = "regression:RosenbrockRegression.FlowTube"
    for row in rows:
        if row[0] == regression_id:
            row[1] = "PASS" if gate_ok else "FAIL"
            row[2] = "TS1 final state matched MICM" if gate_ok else gate_message
        if row[0] in example_results:
            matched, message = example_results[row[0]]
            row[1] = "PASS" if matched else "FAIL"
            row[2] = message
    for row in rows:
        print(" ".join(row))
    return 0 if gate_ok and run.returncode == 0 and all(ok for ok, _ in example_results.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
