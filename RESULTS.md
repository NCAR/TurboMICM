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
- Linear conservation and equilibrium constraint utilities validate
  caller-owned indexed records, construct algebraic masks and Jacobian patterns,
  and add residual and negative-Jacobian contributions for scalar and grouped
  cells. The checks include finite-difference derivatives and invalid records;
  the constrained Rosenbrock entry point projects linear and equilibrium
  algebraic variables from inconsistent initial data and leaves consistent data
  unchanged, reports initial projection iterations, and restores the caller's
  state when projection fails. Row-scaled constraints and repeated solve calls
  also preserve the corrected manifold. The projection tolerance follows ten
  percent of the configured absolute-plus-relative state error scale, and
  backtracking handles nonlinear cold starts. Full DAE stage equations remain
  incomplete. A grouped constrained call projects all
  cells and aggregates statistics while using an independent scalar controller
  per cell.
- Masked reaction forcing and negative-Jacobian accumulation preserve selected
  algebraic rows while updating differential rows across multiple cells.
- The dense matrix view test covers strided indexing, conversion to and from a
  flat vector or individual rows and columns, AXPY, fill, min/max/sum reductions,
  and invalid dimensions. AXPY rejects row and column shape mismatches. Read-only
  views support indexed reads and row conversion. Callback traversal works for
  both nonempty and empty mutable views.
- JSON and YAML tutorial mechanisms normalize to 7 species; the input reactions
  expand to 10 and 7 generated reactions, respectively (JSON and YAML).
- TS1 import contains 210 phase species, 547 reactions, and the expected initial
  defaults. The generated mechanism compiles with both compilers.
- `make test-cpu CC=clang` and `make test-cpu CC=gcc-15`: both emit all 680 case
  statuses. One hundred eighty-three cases pass, nine are marked N/A, and 488 cases
  still report incomplete ports, so both targets exit nonzero.

## TS1 final-state comparison

The seed-1 first cell, 100 steps at `DT = 30 s`, passed the supplied MICM
comparison utility with both compilers. The 100-cell seed-1 run also passed:

| Run | Compiler | C elapsed time | MICM elapsed time | Solver steps | Largest relative difference |
| --- | --- | ---: | ---: | ---: | ---: |
| First cell | Clang | 168.976 ms | 10.767 ms | 900 | `0` |
| First cell | GCC 15 | 156.311 ms | 10.767 ms | 900 | `5.69e-15` |
| 100 cells, 4 threads | Clang | 4,398.141 ms | 354.051 ms | 3,600 | `1.32e-14` |
| 100 cells, 4 threads | GCC 15 | 4,068.024 ms | 354.051 ms | 3,600 | `3.27e-14` |

The 10,000-cell seed-2026 run, 100 steps, and 10 threads completed in
298,417.976 ms; the reference completed in 26,129.420 ms. The comparison passed
all 2,110,000 values; the largest relative difference was `5.34e-14` for
BCARYO2VBS in cell 9034. The C implementation is about 11.4 times slower than
MICM at this scale. The 100-cell and first-cell runs above also pass for both
compilers. Overall, 488 of the 680 inventory cases remain incomplete; nine
C++-only cases are marked N/A.
