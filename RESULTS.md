# Results

## Build and configuration checks

- `make cpu CC=clang`: passed on the provided macOS machine.
- `make cpu CC=gcc-15`: passed on the provided macOS machine.
- `tests/test_core.c`: passed with Clang and GCC 15, including rate formula
  checks across Arrhenius, branched, lambda, reversible, surface, Taylor, Troe,
  ternary activation, tunneling, and user-defined rates; also forcing and small
  analytic integration checks, C reaction construction and value-copy checks,
  density-parameterized reaction multipliers, indexed lambda rates, version
  getters, state setter behavior, and process rate constants, forcing, and
  Jacobians across scalar and grouped cells.
- `tests/test_codegen.py`: passed with both compilers through `make test-cpu`,
  checking rate-type ordering, parameter defaults and indexes, derived branched
  parameters, surface parameters, density multipliers, and reaction counts in
  generated C.
- The C core test ports the README tutorial's Foo/Bar/Baz mechanism and checks
  all ten printed states at the tutorial's three-significant-digit precision.
- The C core test ports `test_solver_results` with three grouped cells and
  checks the printed solver state and counters. Its 18 accepted steps, 36
  function calls, 18 Jacobian updates/decompositions, and 54 solves match MICM.
- The C core test also runs the grouped chemistry for ten 200-second intervals
  and checks all 90 printed concentrations from `test_multiple_grid_cells` at
  the tutorial's three-significant-digit precision.
- The C core test ports `test_openmp` with three worker-owned states and solver
  workspaces, checking each thread's final concentrations against MICM output.
- The C core test reproduces the first ten one-second FlowTube regression states
  for both the two-stage Rosenbrock solver and backward Euler, matching the
  checked-in MICM CSV values within the reference tolerance.
- Stub aerosol coverage uses generic C reactions with the same qualified state
  names. It checks single and multi-cell state updates, forcing, Jacobians, and
  single and three-cell Rosenbrock integrations against the analytical
  partitioning results. Backward Euler tests check single and three-cell solves
  with one aerosol process at 2e-4 tolerance and both processes at 1e-1, matching
  the reference case tolerances.
- The C constrained solver integrates simple equilibrium DAE cases with a
  temperature dependent equilibrium constant, a reordered species array, two
  coupled algebraic species, and non-unit stoichiometry. It checks the
  equilibrium residual after every 0.001 second solve interval.
- A four-stage DAE Rosenbrock run combines the constrained Robertson system
  with the unconstrained chlorine terminator system over twelve decade-spaced
  intervals, checking mass conservation, chlorine conservation, nonnegative
  concentrations, and final-state evolution.
- The C core test ports the vectorized-matrix tutorial's two cell layouts and
  both sparse element orders. The tutorial's solve loop is skipped by its own
  initial converged-state condition, so the layout checks cover its output.
- The C solver accepts configurable Rosenbrock coefficient sets. The core test
  checks all 50 printed concentration rows and accumulated work counters for
  MICM's two-, three-, four-, four-stage DAE-, and six-stage DAE methods.
- The no-user and user-defined rate tutorials now check all ten printed states
  against MICM. YAML/JSON reader tests confirm that the two configurations share
  the same seven base reactions and convert Troe/ternary parameters to the
  values used by the hand-built mechanisms.
- Scalar and grouped solves clamp ordinary state rows to nonnegative values
  after each solver call, matching MICM's `PostSolveClamp` behavior.
- Repeated-reactant Jacobians now match for both separate unit terms (`B + B`)
  and one stoichiometric term (`2 B`), covering a coefficient case that the
  previous negative-Jacobian implementation mishandled.
- Solver workspaces precompute the mechanism-generic Jacobian reactant
  dependencies per species column. The schedule retains reaction accumulation
  order and avoids scanning unrelated reactions on every step; Rosenbrock
  Jacobian clearing and matrix copying also visit only the symbolic sparse
  pattern.
- Linear conservation and equilibrium constraint utilities validate
  caller-owned indexed records, construct algebraic masks and Jacobian patterns,
  and add residual and negative-Jacobian contributions for scalar and grouped
  cells. The checks include finite-difference derivatives and invalid records;
  the constrained Rosenbrock entry point projects linear and equilibrium
  algebraic variables from inconsistent initial data and leaves consistent data
  unchanged, reports initial projection iterations, and restores the caller's
  state when projection fails. Row-scaled constraints and repeated solve calls
  also preserve the corrected manifold. Initial-condition projection follows
  the configured tolerance, while accepted constrained steps solve algebraic
  rows to a strict residual tolerance. Constrained Rosenbrock stages include
  constraint forcing and Jacobians and use a zero mass-matrix diagonal for
  algebraic rows; the embedded step-change error estimate includes those rows.
  Cascade-equilibrium regressions confirm that tightening the algebraic balance
  tolerance increases accepted steps and prevents deep negative overshoot in
  the super-stiff case. A grouped constrained call projects all cells and
  aggregates statistics while using an independent scalar controller per cell.
- Masked reaction forcing and negative-Jacobian accumulation preserve selected
  algebraic rows while updating differential rows across multiple cells.
- The dense matrix view test covers strided indexing, conversion to and from a
  flat vector or individual rows and columns, AXPY, fill, min/max/sum reductions,
  and invalid dimensions. AXPY rejects row and column shape mismatches. Read-only
  views support indexed reads and row conversion. Callback traversal works for
  both nonempty and empty mutable views, can be reused for different row counts,
  and can use vector data through its callback context. Multi-matrix row
  callbacks handle different column counts and reject mismatched row counts.
- The grouped vector-matrix view test checks lane-packed storage, grouped sizes,
  row and column strides, mutable and read-only indexing, vector conversion,
  AXPY, fill, reductions, callbacks, and loops across widths one, two, and four.
  Callback reuse across widths and row counts, including captured vector data,
  is also checked. Row and column vector copies enforce exact dimensions and
  accept empty views with empty vectors. Grouped views also support logical
  matrix copying and multi-matrix row callbacks across lane widths.
- Phase and system metadata checks cover optional diffusion coefficients,
  parameterized-species state-size filtering, and phase-qualified or plain names.
- Species property tests cover temperature/density parameterization and typed
  double, integer, boolean, and string values, including missing-name and
  type-mismatch status codes.
- JSON and YAML tutorial mechanisms normalize to 7 species; the input reactions
  expand to 10 and 7 generated reactions, respectively (JSON and YAML).
- TS1 import contains 210 phase species, 547 reactions, and the expected initial
  defaults. The generated mechanism compiles with both compilers.
- `make test-cpu CC=clang` and `make test-cpu CC=gcc-15`: both emit all 680 case
  statuses. 662 cases pass and 18 C++-only cases are marked N/A; neither target
  reports a failed case. Both compiler runs were forced full rebuilds.

## TS1 final-state comparison

Seed-1 and seed-2026 final states passed the supplied MICM comparison utility.
The 10,000-cell run exercised the full required scale with both compilers:

| Run | Compiler | C elapsed time | MICM elapsed time | Accepted steps | Largest relative difference |
| --- | --- | ---: | ---: | ---: | ---: |
| First cell, seed 1 | Clang | 38.237 ms | 8.722 ms | 900 | `0` |
| 100 cells, seed 1, 4 threads | Clang | 787.578 ms | 374.132 ms | 3,600 | `1.32e-14` |
| 10,000 cells, seed 2026, 10 threads | Clang | 48,528.301 ms | 27,298.520 ms | 9,000 | `5.34e-14` |
| 10,000 cells, seed 2026, 10 threads | GCC 15 | 55,820.845 ms | 27,298.520 ms | 9,000 | `1.43e-13` |

Both full runs passed all 2,110,000 values (Clang's largest relative difference
was for BCARYO2VBS in cell 9034; GCC's was for BCARYO2VBS in cell 4054). The
Clang implementation is 1.8 times slower than MICM on the 10,000-cell run. The
final inventory has 662 passing cases and eighteen C++-only cases marked N/A.
