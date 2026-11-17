"""Compare standalone tutorial output to the copied reference executables."""
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parents[1]
programs = (
    ("rate_constants_no_user_defined_example_by_hand", "tutorial_no_custom.expected"),
    ("rate_constants_user_defined_example_by_hand", "tutorial_custom.expected"),
    ("README_example", "readme.expected"),
    ("multiple_grid_cells", "multiple_grid_cells.expected"),
    ("solve_results", "solve_results.expected"),
)
for program, expected_name in programs:
    result = subprocess.run([root / "build" / program], capture_output=True, check=True)
    expected = (root / "tests" / "reference" / expected_name).read_bytes()
    if result.stdout != expected:
        print(f"{program}: output differs from the reference", file=sys.stderr)
        sys.exit(1)
print("Tutorial programs match the reference output")
