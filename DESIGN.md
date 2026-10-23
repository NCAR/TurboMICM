# Design

## Tier 2 mechanism specialization

`tools/generate_c.py` is the authoritative code generator. It orders reactions
and species at generation time and emits the same immutable mechanism tables
used by Tier 1. With `CHEM_SPECIALIZED` enabled, that source also provides fixed
TS1 operations:

- forcing is emitted as straight-line per-reaction expressions with fixed
  species indexes, reactant powers, yields, and coefficients;
- negative-Jacobian updates are emitted in species-column, reaction, and
  reactant-term order, preserving accumulation order and the generic library's
  diagonal coefficient handling;
- the generator constructs the same symbolic fill pattern as `make_pattern`,
  then emits natural-order LU multipliers, updates, and forward/backward solve
  terms with constant indexes.

The Rosenbrock controller, stage coefficients, error estimate, tolerances,
step-size adaptation, rate-constant evaluation, and rate-constant cache remain
the Tier 1 implementation. Temperature, pressure, air density, concentrations,
custom rate parameters, duration, and tolerances remain runtime data. The
specialized callbacks consume cached constants from the generic evaluator and
do not alter rate formulas. A generic mechanism has a null specialization table
and follows the table-driven operations.

The generated LU uses the same no-pivot natural ordering and symbolic fill as
the generic factorization. It checks every pivot for zero and non-finite values
and returns the same solver status. The generated triangular solve visits lower
and upper entries in the same ascending column order as the generic solve.
Callbacks are guarded by `CHEM_SPECIALIZED`, so the generic binary does not
compile the expanded schedules.

`make bench` builds the generic Tier 1 benchmark, the TS1 specialized benchmark,
and the specialized scenario runner. Scenario runs use explicit cell groups,
with one adaptive controller per group. `--group-size` and `--threads` are
runtime parameters; when omitted, one group is assigned to each worker.
Compiler optimization flags remain controlled by `OPTFLAGS`.

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
change. Solver workspaces also precompute a mechanism-generic per-species list
of reaction reactant terms for the negative-Jacobian accumulation. Each column
still visits those dependencies in reaction order, preserving the reference
summation sequence without scanning unrelated reactions on every step.

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
Typed species properties are also caller-owned; getters return distinct status
codes for missing names and mismatched types.

`chem_constraint_set` describes caller-owned linear conservation or
equilibrium constraints using species indexes. Its utilities validate records,
mark algebraic rows, construct the constraint portion of the Jacobian pattern,
and add residual and negative-Jacobian terms to cell-major forcing and dense
Jacobian arrays. Equilibrium constraints use a van't Hoff temperature correction.
The masked process helpers add reaction forcing and negative-Jacobian values to
caller arrays while skipping algebraic rows.
The constrained Rosenbrock solver adds constraint residuals and Jacobians to
the stage equations, omits the alpha shift and interstage mass coupling on
algebraic rows, and includes algebraic step changes in the normalized error
estimate. It also projects the initial state and each accepted step by Newton
iteration, solving only for the algebraic variables. Initial projection uses
the configured state tolerance and reports its iteration count; accepted-step
projection uses a strict residual tolerance to preserve the constraint
manifold. A backtracking line search keeps trial algebraic concentrations
nonnegative and requires the normalized constraint residual to decrease.
The constrained group entry point currently runs cells sequentially with an
independent scalar controller for each cell.

## Solver and benchmark

`src/solver.c` implements MICM's two-, three-, four-, four-stage DAE-, and
six-stage DAE Rosenbrock parameter sets, adaptive error control, symbolic fill
construction, natural-order in-place LU, and solver counters. It also provides
a fully implicit backward Euler path using Newton corrections, nonnegative
state updates, the configurable sequence of failed-step reductions, and the
state residual convergence test. Backward Euler uses the same generic mechanism
forcing and Jacobian evaluators and can be configured with one absolute
tolerance per species. Its constraints remain unsupported, matching MICM's
backward Euler policy. The public
`chem_rosenbrock_parameters` stores stage count, packed lower-triangular A and C
coefficients, M and E vectors, gamma, function-evaluation reuse flags, and
error-estimator order. Add a method by filling this structure and passing it to
`chem_rosenbrock_solve_with_parameters`; the new set uses the generic solver
path.

The default controller uses relative tolerance `1e-6`, absolute tolerance
`1e-3`, initial step `1e-6 * DT`, minimum step `1e-15 * DT`, maximum step `DT`,
and a 1000-step limit. Set `initial_step` to override the starting step, and
use `chem_mechanism_absolute_tolerances` to build the per-species absolute
tolerance vector from `absolute tolerance` species properties (default `1e-3`).
The inert third-body species is excluded from the error norm. Ordinary state
rows are clamped to nonnegative values after each solve call, as MICM does in
its `PostSolveClamp` step.

The public `chem_solver_workspace` owns reusable matrix, stage, symbolic fill,
and group-output buffers for one mechanism. The single-cell API uses one cell
per workspace; the benchmark initializes one workspace per cell group. `src/main.c`
loads CSV values by species or parameter label, initializes omitted values
from generated defaults, runs the untimed warm-up and timed steps, and writes
`ts1_final.csv`. `chem_version*` exposes the C library version as `0.1.0+0`.
`--threads N` creates N contiguous cell groups. Each group
shares one adaptive step controller and error norm, matching MICM's benchmark
grouping. The current C implementation still factors each cell's matrix
sequentially inside the group. The generic path's Jacobian zeroing, matrix
copying, factorization, and solve use its symbolic sparse pattern. The
specialized path keeps the same per-cell execution model but replaces the
factor and solve loops with generated constant-index schedules.

## Extension points

The Rosenbrock controller tolerances live in `chem_solver_options`; backward
Euler parameters and its convergence predicate use a separate public options
structure. Both methods share the workspace factorization buffers. Their
controller implementations are separate from the method coefficient tables,
so either can be changed independently. The factorization is isolated from rate
evaluation and can be replaced behind a linear-solver interface. Add a rate
type by extending the mechanism enum, its evaluator and Jacobian handling, and
the configuration-to-C generator. Configuration cases should exercise the
generated path, not a handwritten mechanism-specific path.

## Tuning and portability

Runtime settings are thread count (default: available OpenMP workers), scenario
group size (default: one group per worker), and solver tolerances (MICM defaults
unless passed to the scenario runner). `OPTFLAGS` defaults to `-O3 -DNDEBUG`.
The C
sources use ISO C11, libm, and CPU OpenMP only. The Makefile selects OpenMP flags
for the chosen compiler and takes optimization flags solely from `OPTFLAGS`.
