# Design

## Mechanisms

Python readers accept normalized JSON and CAMP configurations. The CAMP reader
resolves referenced files relative to the input file, then normalizes phases,
species, reactions, and rate parameters. `tools/import_ts1.py` translates the
permitted TS1 builder input. `tools/generate_c.py` emits immutable C tables from
that normalized representation. Build-time generation keeps Python, YAML, and
configuration libraries out of the executable.

The C API separates mechanism data from runtime computation. Species and
reaction tables carry metadata and index-based stoichiometry. Third bodies and
general density-parameterized species fold their density factors into rate
constants while remaining distinct from kinetic state dependencies. Evaluation
uses a rate-kind enum. `chem_forcing` and `chem_jacobian` apply
mass action and the mechanism's generic reaction terms. The same routines feed
the generic three-stage Rosenbrock integration path exercised by TS1.

## Solver and benchmark

`src/solver.c` implements the three-stage method, adaptive error control,
symbolic fill construction, natural-order in-place LU, and solver counters.
Its defaults are relative tolerance `1e-6`, absolute tolerance `1e-3`, initial
step `1e-6 * DT`, minimum step `1e-15 * DT`, maximum step `DT`, and a 1000-step
limit. The inert third-body species is excluded from the error norm.

The public `chem_solver_workspace` owns reusable matrix, stage, symbolic fill,
and group-output buffers for one mechanism. The single-cell API uses one cell
per workspace; the benchmark initializes one workspace per cell group. `src/main.c`
loads CSV values by species or parameter label, initializes omitted values
from generated defaults, runs the untimed warm-up and timed steps, and writes
`ts1_final.csv`. `--threads N` creates N contiguous cell groups. Each group
shares one adaptive step controller and error norm, matching MICM's benchmark
grouping. The current C implementation still factors each cell's matrix
sequentially inside the group.

## Extension points

The solver method coefficients and controller are currently localized in
`src/solver.c`; they should be moved into explicit method and controller
interfaces before supporting alternate methods. The factorization is similarly
isolated from rate evaluation and can be replaced behind a linear-solver
interface. Add a rate type by extending the mechanism enum, its evaluator and
Jacobian handling, and the configuration-to-C generator. Configuration cases
should exercise the generated path, not a handwritten mechanism-specific path.

## Tuning and portability

The current runtime setting is thread count, defaulting to available OpenMP
workers; it also sets the number of cell groups. Integration tolerances and
step limits are solver defaults. The C
sources use ISO C11, libm, and CPU OpenMP only. The Makefile selects OpenMP flags
for the chosen compiler and takes optimization flags solely from `OPTFLAGS`.
