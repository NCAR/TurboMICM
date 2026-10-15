# Completion ledger

The full task is active. This file records checkpoints, not a reduced definition
of completion. `tests/cases.txt` is the authoritative inventory of all 680 cases.
`tools/run_tests.py` emits FAIL for each missing case implementation.

## Implemented and checked

- Local, pinned conda environment; portable Clang/GCC build and compiler-setting
  invalidation; only C11/libm/OpenMP in the library and benchmark.
- All analytic rate types, callback rates, parameterized reactants, owned
  mechanism metadata, state allocation/copy/setters, forcing and dense Jacobian.
- YAML/JSON tutorial reader, generic C builder generation, retained TS1 importer
  and configuration (210 species, 547 reactions, 209 solution variables).
- All three tutorial configurations compared with hand-built mechanisms at 1
  and 7 cells for rates and ten integration steps.
- Sparse symbolic fill; in-place Mozart and Doolittle LU; triangular solves;
  Markowitz permutation; five original Rosenbrock parameter sets; ODE
  Rosenbrock steps, error calculation and controller. TS1 uses this same solver.
- Benchmark CLI, warm-up, static OpenMP groups, statistics, defaults, CSV input
  and final-state output. Group size is configurable at runtime.
- 31 rate cases, 26 in-place LU/linear-solver cases and 2 backward Euler cases
  ported with original or
  tighter numerical tolerances. No N/A declarations yet.
- 100-cell, 100-step, 30-second TS1 gate passes on Clang and GCC. GCC group size
  7 with 4 threads also passes. Small Clang reference was rebuilt with the same
  flags and compared successfully. Full 10,000-cell seed-2026 gate also passes
  all 2,110,000 values with Clang (identical) and GCC (maximum 1.43e-13 relative).
- Backward Euler Newton iteration, convergence test, step reductions and
  status handling are implemented in `src/backward_euler.c`.
- Public solver post-interval concentration clamping is implemented; a retained
  TS1 fixture exercises the negative-intermediate-state case.

## Remaining work (required)

- Port all remaining 621 cases, individually checking original test quantities.
  Read each source/policy; do not assign a PASS based on a related smoke test.
- Constraints and their consistent initial-condition Newton/line-search
  algorithm, algebraic error control, external model hooks and aerosol models.
- Backward Euler analytical/regression coverage, reaction-rate tracking, solver status
  behavior, solver-builder validations and full species/property/phase APIs.
- Dense/sparse matrix utilities, grouped views and operations, generic C
  equivalents of C++ ownership and dimension checks. Only genuinely
  compile-time-only cases may receive documented N/A reasons.
- Out-of-place LU paths and their tests, solver regression/analytical tests,
  all eight complete tutorial programs with reference stdout comparisons.
- Validate all five Rosenbrock parameter sets and rejected-step behavior. The
  TS1 gate exercises only ROS3 with no rejections.
- Validate the configuration reader's malformed inputs and additional rate
  types beyond the three supplied examples; compare TS1 default fields by name.
- Performance work using generic kernels; broader seeds/group sizes and both
  compilers; full final portability/source audit and complete RESULTS.md.
- Final commit and requirement-by-requirement completion audit. Mark goal
  complete only after all requirements are proven.

## Implementation notes for continuing work

`src/internal.h` holds internal mechanism/process structs. `chem_mechanism_reorder`
mutates index maps and must be called before states/solvers exist. The mechanism
pattern excludes forced diagonals during Markowitz selection, matching MICM;
the sparse constructor adds the diagonal for numerical factorization.

`src/solver_internal.h` shares the workspace between integrators.
`src/rosenbrock.c` builds derivative records ordered by independent variable then
reaction and evaluates the negative sparse Jacobian. Dense `chem_jacobian`
returns the positive derivative. The solver owns scratch matrices; states own
cell workspaces and must not be used concurrently. Independent states can share
a mechanism. Linear-solver and step-controller callbacks are public interfaces.

Existing generated objects and validation files are under ignored `build/`.
All Python invocations must use `.conda/bin/python`. Read-only input trees are
`../micm` and `../reference`; do not inspect other sibling directories. Run MICM
tests only from a copy of its build tree inside this repository. The explicit
`make micm-ref` reference executable is the only linked MICM benchmark allowed.
