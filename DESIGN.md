# Design

## Mechanisms

Python readers accept normalized JSON and CAMP configurations. The CAMP reader
resolves referenced files relative to the input file, then normalizes phases,
species, reactions, and rate parameters. `tools/import_ts1.py` translates the
permitted TS1 builder input. `tools/generate_c.py` emits immutable C tables from
that normalized representation. Build-time generation keeps Python, YAML, and
configuration libraries out of the executable.

The C API separates mechanism data from runtime computation. `chem_state_view`
provides non-owning, cell-major access to concentration and custom-parameter
arrays, with label-based setters and status-code validation. Batch updates validate
all labels and lengths before writing. Species and
reaction tables carry metadata and index-based stoichiometry. Third bodies and
general density-parameterized species fold their density factors into rate
constants while remaining distinct from kinetic state dependencies. Evaluation
uses a rate-kind enum. `chem_forcing` and `chem_jacobian` apply
mass action and the mechanism's generic reaction terms. The same routines feed
all solver schemes and TS1. Rate constants are evaluated before integration
and cached per workspace until temperature, pressure, density, or custom rates
change.

`chem_matrix_view` is a non-owning rectangular view with explicit row and column
strides. It provides checked element access, fill, conversion to and from a
row-major vector, shape-checked AXPY, and min/max/sum reductions. It supports
strided storage without imposing an allocator or matrix ownership policy.
Element traversal accepts a C callback and caller context.
Row callbacks can read two source matrices with different column counts and
write into a destination matrix, while validating that all row counts agree.
`chem_const_matrix_view` provides checked read-only access and row conversion
when the caller holds const storage.
`chem_vector_matrix_view` packs logical rows into fixed-width lanes, with the
same row and column strides as MICM's grouped vector matrices.
Its row callback accepts multiple source matrices with distinct column counts
and lane widths.

`chem_phase` and `chem_system` are non-owning compositions of caller-owned
species records. Phase state sizes omit density-parameterized species, and
phase/system name helpers emit qualified or plain names into caller buffers.

`chem_constraint_set` describes caller-owned linear conservation or
equilibrium constraints using species indexes. Its utilities validate records,
mark algebraic rows, construct the constraint portion of the Jacobian pattern,
and add residual and negative-Jacobian terms to cell-major forcing and dense
Jacobian arrays. Equilibrium constraints use a van't Hoff temperature correction.
The masked process helpers add reaction forcing and negative-Jacobian values to
caller arrays while skipping algebraic rows.
These are evaluation primitives; the solver does not yet use them for DAE
stage equations or algebraic-row replacement. The optional constrained solve
projects algebraic variables by Newton iteration at the start and after each
accepted ODE step, using the constraint Jacobian restricted to algebraic columns.
The initial Newton iteration count is reported in solver statistics.
Projection stops when the algebraic correction is at most ten percent of the
configured absolute-plus-relative state error scale. A backtracking line search
keeps trial algebraic concentrations nonnegative and requires the normalized
constraint residual to decrease.
The constrained group entry point currently runs cells sequentially with an
independent scalar controller for each cell.

## Solver and benchmark

`src/solver.c` implements MICM's two-, three-, four-, four-stage DAE-, and
six-stage DAE Rosenbrock parameter sets, adaptive error control, symbolic fill
construction, natural-order in-place LU, and solver counters. The public
`chem_rosenbrock_parameters` stores stage count, packed lower-triangular A and C
coefficients, M and E vectors, gamma, function-evaluation reuse flags, and
error-estimator order. Add a method by filling this structure and passing it to
`chem_rosenbrock_solve_with_parameters`; the new set uses the generic solver
path.

The default controller uses relative tolerance `1e-6`, absolute tolerance
`1e-3`, initial step `1e-6 * DT`, minimum step `1e-15 * DT`, maximum step `DT`,
and a 1000-step limit. The inert third-body species is excluded from the error
norm. Ordinary state rows are clamped to nonnegative values after each solve
call, as MICM does in its `PostSolveClamp` step.

The public `chem_solver_workspace` owns reusable matrix, stage, symbolic fill,
and group-output buffers for one mechanism. The single-cell API uses one cell
per workspace; the benchmark initializes one workspace per cell group. `src/main.c`
loads CSV values by species or parameter label, initializes omitted values
from generated defaults, runs the untimed warm-up and timed steps, and writes
`ts1_final.csv`. `chem_version*` exposes the C library version as `0.1.0+0`.
`--threads N` creates N contiguous cell groups. Each group
shares one adaptive step controller and error norm, matching MICM's benchmark
grouping. The current C implementation still factors each cell's matrix
sequentially inside the group.

## Extension points

The controller tolerances live in `chem_solver_options` and its adaptation
formula is a separate block in `src/solver.c`, so it can be changed without
changing method coefficients. The factorization is isolated from rate
evaluation and can be replaced behind a linear-solver interface. Add a rate
type by extending the mechanism enum, its evaluator and Jacobian handling, and
the configuration-to-C generator. Configuration cases should exercise the
generated path, not a handwritten mechanism-specific path.

## Tuning and portability

The current runtime setting is thread count, defaulting to available OpenMP
workers; it also sets the number of cell groups. Integration tolerances and
step limits are solver defaults. The C
sources use ISO C11, libm, and CPU OpenMP only. The Makefile selects OpenMP flags
for the chosen compiler and takes optimization flags solely from `OPTFLAGS`.
