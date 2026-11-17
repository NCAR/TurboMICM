# CPU chemistry port (in progress)

Create the local tooling environment with `make env`. `make cpu` builds the C11
mass-action library. `make check CC=clang` or `make check CC=gcc-15` builds and
runs the current focused tests; the tests read tutorial configuration input from
the adjacent `micm` checkout. `OPTFLAGS` can be supplied to make.

A preliminary single-cell Rosenbrock solver with five parameter sets is available through
`include/rosenbrock.h`. The required full test suite, remaining rate types, additional
solver features, and TS1 benchmark remain to be implemented. Run `make clean`
before switching compilers.

Regenerate the TS1 configuration from the reference source with
`.conda/bin/python tools/import_ts1.py ../micm/test/benchmark/ts1_mechanism.hpp configs/ts1`.
The current checks compare the translation's initial values with the reference
defaults CSV and exercise generated TS1 rate and forcing evaluation.

Build the preliminary TS1 driver with `make bench`, then run
`bin/ts1_bench_cpu [--threads N] CELLS STEPS DT CELLS_CSV`. It writes
`ts1_final.csv`. `make check-ts1` compares 100 cells and 100 steps with the
supplied reference, and `make micm-ref` builds the independent timing reference.
The current solver uses independent per-cell integration and sparse no-pivot
factorization; further optimization and broader verification are still needed.

`make test-cpu` reports every ID from the reference case list. Currently
588 cases have been ported, so the target intentionally fails and labels the
other 92 as unported. The focused `make check` tests cover additional behavior
but do not claim full equivalence with those remaining cases.

The single-cell backward Euler API is declared in `include/backward_euler.h`.
It shares forcing, Jacobian, and optional sparse linear-solver interfaces with
the Rosenbrock solver.

`include/constraint.h` provides indexed linear and equilibrium constraints,
validation, residuals, Jacobians, sparsity, and algebraic error estimates.
Rosenbrock options accept a constraint set and per-cell equilibrium parameters
or temperature for index-one DAE integration. Sparse plans for constrained
systems must include both kinetic and constraint Jacobian structure.

`include/state.h` provides an optional owned, multi-cell state with named
concentration and custom parameter updates. Its arrays can be passed to the
existing rate and integration APIs. Initialize states to zero before use and
release them with `chem_state_destroy`. External models can contribute names
and tolerances through `chem_state_init_composed` and resolve their indices
through the state lookup functions.

`include/species.h` provides owned species, typed properties, phase metadata,
parameterizations, diffusion coefficients, and gas-system name queries. These
metadata objects are separate from the indexed reaction tables used by the
numerical kernels.

Custom DAE constraints can provide residual and negative-Jacobian callbacks
and declare their state dependencies through `ChemCustomConstraint`.

`include/matrix.h` provides owned row-major or interleaved group matrices,
with checked row and column operations and configurable group width.

`include/sparse_matrix.h` provides compressed row and column sparse blocks
with optional interleaved cell groups.
