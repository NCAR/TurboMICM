"""Compare standalone tutorial output to the copied reference executables."""
from pathlib import Path
import re
import subprocess
import sys

root = Path(__file__).resolve().parents[1]
programs = (
    ("rate_constants_no_user_defined_example_by_hand", "tutorial_no_custom.expected"),
    ("rate_constants_user_defined_example_by_hand", "tutorial_custom.expected"),
    ("README_example", "readme.expected"),
    ("multiple_grid_cells", "multiple_grid_cells.expected"),
    ("solve_results", "solve_results.expected"),
    ("solver_configuration", "solver_configuration.expected"),
)
for program, expected_name in programs:
    result = subprocess.run([root / "build" / program], capture_output=True, check=True)
    expected = (root / "tests" / "reference" / expected_name).read_bytes()
    actual = result.stdout
    if program == "solver_configuration":
        pattern = rb"Total solve time: ([0-9.eE+-]+) nanoseconds"
        actual_times = re.findall(pattern, actual)
        expected_times = re.findall(pattern, expected)
        if len(actual_times) != 5 or len(expected_times) != 5 or any(float(t) < 0 for t in actual_times):
            print(f"{program}: invalid timing output", file=sys.stderr)
            sys.exit(1)
        actual = re.sub(pattern, b"Total solve time: <time> nanoseconds", actual)
        expected = re.sub(pattern, b"Total solve time: <time> nanoseconds", expected)
    if actual != expected:
        print(f"{program}: output differs from the reference", file=sys.stderr)
        sys.exit(1)
print("Tutorial programs match the reference output")
