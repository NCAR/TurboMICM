# Reference fixtures

The eight `.txt` files are stdout captured from the corresponding supplied
MICM tutorial executables at commit `e15898ce` (v3.7.0-259), in a copy of the
provided build directory inside this repository. The original build directory
was not run or modified. These files are data, not implementations.

`cases_tutorials.py` compares every deterministic line, including printed
concentrations, dimensions, and solver statistics. Wall-clock nanoseconds are
validated as numeric timing lines; their values cannot be identical between
implementations. OpenMP startup messages are compared as a set of thread IDs
because their order is nondeterministic. All other line content must match.

The two wall-loss CSV files are copied from
`micm/test/regression/expected_results/`. The C port checks the full 3,601-row
trajectory at the original absolute tolerance of 1e-10, including the ten rows
checked by MICM's regression helper.

MICM source and fixtures: Copyright (C) 2023-2026 University Corporation for
Atmospheric Research, Apache-2.0. See the repository's `LICENSE`.
