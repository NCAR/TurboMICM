# Design (work in progress)

The generic linear workspace can also be constructed directly from a matrix
pattern. It exposes the same symbolic fill, Mozart and Doolittle numeric
factorizations, and triangular solves used by the integration workspace.

`tools/mechanism.py` reads and validates the configuration, then
`tools/generate.py` emits immutable C tables. `chem.h` also exposes the table
format for mechanisms constructed directly in C. The tables list repeated
reactants explicitly, preserving the product rule even at zero concentration.
The generator chooses a diagonal minimum-degree species ordering from the
Jacobian graph to reduce factorization fill. Names in the output preserve the
mapping back to the input.

`matrix.h` supplies contiguous and grouped matrix storage with runtime group
width, row functions, column operations, and reductions. `sparse.h` provides
compressed row and column layouts, grouped blocks, and block functions. Function objects retain
expected column shapes and validate row and vector sizes on each invocation.
Dense and sparse group views expose contiguous active entries; a partial group's
count excludes padded entries. Mutable and const views have distinct C types.
`state.h` stores per-cell conditions, variables, custom parameters, and rate
constants; names are borrowed while data and tolerances are owned. Bulk
setters validate all labels and dimensions before changing the state.
Rate updates use each cell's conditions and custom parameters. Process-set
operations also expose additive forcing, negative Jacobians, symbolic patterns,
and optional algebraic-row masks.

Rate constants, forcing, and Jacobian evaluation are separate from the solver.
External models can add forcing and Jacobian entries through callbacks and
declare their symbolic pattern; they may also update per-cell parameters and
provide algebraic residuals and Jacobians in named rows. External process
entries remain represented in the pattern even when constraints replace the
corresponding numerical rows.
The group Rosenbrock interface shares a time step and error norm across cells,
while retaining a workspace and external context for each cell.
The solver builder validates configuration and calls model finalizers. Solver
options can be replaced for a solve; grouped Rosenbrock states share control,
and backward Euler advances cells independently.
`ChemRosenbrockParameters` contains stage coefficients and control parameters;
two-, three-, and four-stage ODE parameter sets are provided.
`chem_rosenbrock_solve` implements the stage loop and error controller. The
workspace derives the symbolic Jacobian and fill once, then uses in-place
Mozart or Doolittle elimination selected in the solver parameters. A replacement linear solver can use the same Jacobian and
stage interfaces, though the interface should be made pluggable before swapping
it. A new rate type requires an enum member, parameter interpretation in the
rate calculation, generator mapping, and tests. C callers may also install
named rate callbacks and replace them between rate updates.

Linear and equilibrium constraints expose validation, residuals, and negative
Jacobian rows independently or over a set of cells. A constrained workspace borrows the constraint
and term arrays, adds their dependencies to the sparse pattern, and replaces
algebraic rows during Rosenbrock stages. Initialization uses weighted Newton
corrections with a backtracking line search and restores the input on failure.
Four- and six-stage DAE parameter sets are available. This path has preliminary
coverage and does not yet satisfy the full constraint test inventory.

Backward Euler uses Newton iterations, nonnegative iterates, and the reference
step-reduction sequence. Both solvers share the chemistry and linear algebra.

The benchmark uses independent solver state per cell and OpenMP static
scheduling. `--threads N` controls the thread count; absent the option OpenMP's
maximum is used. Each thread owns a reusable workspace. Current Jacobians use
dense storage with a symbolic sparse pattern and precomputed row and column
indices for factorization and substitution; group size is one cell. The
initial step is `1e-6 * DT`, the minimum is `1e-15 * DT`, relative tolerance is
`1e-6`, and absolute tolerance is `1e-3`, matching the reference defaults.
These solver values can be set in `ChemRosenbrockParameters`, including an
optional per-species absolute tolerance array.
