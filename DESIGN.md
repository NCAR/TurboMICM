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

Rate constants, forcing, and Jacobian evaluation are separate from the solver.
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
