# CPU chemistry port (in progress)

Create the local tooling environment with `make env`. `make cpu` builds the C11
mass-action library and TS1 benchmark. `make check CC=clang` or `make check CC=gcc-15` builds and
runs the current focused tests; the tests read tutorial configuration input from
the adjacent `micm` checkout. `OPTFLAGS` can be supplied to make. Set `CONDA`
if Conda is not on `PATH`; on macOS with clang, `LIBOMP_PREFIX` can override
the Homebrew libomp location.

Rosenbrock and backward Euler solvers are available through
`include/rosenbrock.h` and `include/backward_euler.h`. Run `make clean`
before switching compilers.

Regenerate the TS1 configuration from the reference source with
`.conda/bin/python tools/import_ts1.py ../micm/test/benchmark/ts1_mechanism.hpp configs/ts1`.
The current checks compare the translation's initial values with the reference
defaults CSV and exercise generated TS1 rate and forcing evaluation.

Build the TS1 driver with `make bench`, then run
`bin/ts1_bench_cpu [--threads N] CELLS STEPS DT CELLS_CSV`. It writes
`ts1_final.csv`. `OMP_SCHEDULE` selects the OpenMP schedule and chunk size;
without it, the benchmark uses static scheduling. `make check-ts1` compares 100 cells and 100 steps with the
supplied reference, and `make micm-ref` builds the independent timing reference.
The current solver uses independent per-cell integration and sparse no-pivot
factorization. Its default generated ordering follows the reference process
and state order. To try the alternate diagonal Markowitz ordering, run
`make clean` and rebuild with `TS1_ORDER_FLAGS=--reorder`.

`make test-cpu` reports every ID from the reference case list. Currently
675 cases pass and five C++ template-dispatch cases are N/A. The focused
`make check` tests also cover additional behavior.

The backward Euler API is declared in `include/backward_euler.h`. It supports
single-cell and shared-step grouped solves and shares forcing, Jacobian, and
optional sparse linear-solver interfaces with the Rosenbrock solver.

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
through the state lookup functions. `chem_state_init_constrained` registers
constraint names and algebraic rows; its temperature update evaluates the
equilibrium parameters for each cell.

`include/model.h` composes external process callbacks, state fragments,
constraints, and sparse patterns. Model contexts and constraint data remain
caller-owned; the resulting state and sparse structure are owned by the model
system.

`include/species.h` provides owned species, typed properties, phase metadata,
parameterizations, diffusion coefficients, and gas-system name queries. These
metadata objects can be compiled with `include/process.h` into indexed
reaction tables used by the numerical kernels.

Custom DAE constraints can provide residual and negative-Jacobian callbacks
and declare their state dependencies through `ChemCustomConstraint`.

`include/matrix.h` provides owned row-major or interleaved group matrices,
with checked row and column operations and configurable group width.
`include/typed_matrix.h` adds arbitrary trivially copyable element types and
const input and mutable output row views for callbacks.

`include/sparse_matrix.h` provides compressed row and column sparse blocks
with optional interleaved cell groups. `include/typed_sparse_matrix.h` provides
the same layouts for trivially copyable values; pointers stored there remain
caller-owned.

`include/mechanism_metadata.h` stores owned reaction, photolysis, and species
names for any mechanism. `include/reorder.h` derives a runtime Markowitz
permutation from kinetic and constraint sparsity and remaps mechanism and
constraint indices. Custom constraint callbacks require a caller-provided
adapter when the state is reordered.

`include/cpu_solver.h` compiles a system and reactions into an owned solver,
creates compatible states, updates rate parameters, and solves with either
Rosenbrock or backward Euler. A positive group width limits cells per solve;
zero permits arbitrary cell counts. Solve-time option overrides persist.
