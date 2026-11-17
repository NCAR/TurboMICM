# Design

## Data flow

`tools/config.py` resolves `camp-files`, validates species and reaction names,
and normalizes the examples' JSON and YAML forms. `tools/generate.py` converts
the normalized mechanism into immutable C tables. It converts molecule-based
rate coefficients to molar units, expands branched reactions, and can compute
a symmetric minimum-fill species order from the reaction graph. The generated
tables use the same `ChemMechanism` representation that callers can construct
in C. The tests compare generated and hand-built mechanisms and also compare a
reordered mechanism with its original form. Optional configuration constraints
produce a separate `ChemConstraint` array and count, with all row and term
indices remapped together when ordering is enabled.

The C library stores each state variable, rate, and sparse-matrix entry across
the cells of one group contiguously. Forcing and Jacobian evaluation use the
reaction tables; symbolic sparse elimination is prepared once and shared by
states created with `chem_state_create_like`. Each state owns its numerical
workspaces. OpenMP SIMD loops operate across cells. The benchmark partitions
cells into groups with independent adaptive solver state and parallelizes
groups with OpenMP. There is no mechanism-specific runtime code.

Rate updates are explicit: set concentrations, conditions, custom parameters,
and any density overrides, then call `chem_update_rates` before solving. The
mechanism and constraint arrays and callback contexts are borrowed and must
outlive their states. Mutable arrays returned by the state API are per state.
The returned species and parameter arrays contain one value per cell; the
absolute-tolerance array contains one value per species.
`chem_state_clone` copies mutable values and rate data while sharing the
immutable symbolic plan; callback contexts and mechanism data remain borrowed.
Named bulk setters check every label and vector length before changing a state.
Parameterized reactants may supply linear functions of temperature, pressure,
and air density; their values multiply the associated reaction rate during a
rate update. A stable optional rate index groups reaction IDs by rate kind.
`chem_system.h` supplies owned, typed species metadata and borrowed phase and
system views for hand-built mechanisms. Phase export creates the runtime
species array, excluding parameterized species. Checked state construction
can reject species unused by native reactions; extended mechanisms may instead
use their model declarations when validating dependencies.
Native-reaction diagnostics can expose the unfilled Jacobian pattern or add
forcing and negative-Jacobian contributions to existing arrays, with optional
algebraic-row masking.

`chem_matrix.h` provides typed-by-element-size dense matrices with row-major or
grouped storage, checked row and column operations, reusable matrix functions,
and strict or padded reductions. These utilities have independent reference
tests; the chemistry kernel keeps its own specialized species-major workspace.
`chem_sparse.h` adds compressed-row and compressed-column structures whose
blocks use the same row-major or grouped storage. Sparse coordinates are
validated and deduplicated when a structure is created.
Dense columns and sparse block entries expose explicit grouped views with
checked fill and copy operations, including partial final groups.

## Numerical components

The Rosenbrock implementation has the two-, three-, four-, and differential-
algebraic four- and six-stage parameter sets. `ChemRosenbrock` holds the
coefficients, tolerances, controller bounds, and projection controls. Its
optional step-controller callback supplies the next step-size factor; a null
callback uses the default formula. Backward Euler has its own parameters and
entry point. Both methods use the same forcing and linear-solver interfaces.
The Jacobian diagnostics expose finite differences, numerical comparisons, and
sparsity checks for independently supplied forcing functions.

The default linear solver performs an in-place sparse Mozart factorization and
triangular substitutions. `chem_doolittle_linear_solver` supplies an
alternative Doolittle implementation. A caller can install a `ChemLinearSolver`
pair of factor/solve callbacks; the public sparse view exposes the row offsets,
column indices, lookup table, values, and algebraic-row mask. A factor callback
receives the diagonal shift and must apply it only to differential rows.
`chem_linear_factor_copy` provides a separate-factor path that preserves the
original matrix; the default state solver uses precomputed in-place traversal.

Linear and equilibrium constraints replace designated species rows. The
Rosenbrock path projects inconsistent initial conditions using weighted Newton
corrections and backtracking, restores the caller's state on projection
failure, applies the mass-matrix diagonal in stage coupling, and leaves
algebraic variables unclamped after solving. Equilibrium constants are updated
from temperature alongside rate constants. The current backward-Euler entry
point rejects constrained states.

External models can add species, parameters, process forcing, Jacobian
entries, and algebraic constraints. `chem_mechanism_compose` appends declared
species and parameters to a base mechanism; the resulting names can then be
resolved to indices for the models' sparse patterns. `chem_state_create_extended`
combines those patterns with the built-in reaction and constraint patterns,
then calls each model's finalizer with the resulting sparse view. Its update
callback runs with rate updates, and its forcing and negative-Jacobian callbacks
run during either integrator. Constraint callbacks supply a residual and its
gradient for a declared algebraic row. Process entries in algebraic rows remain
in the sparse pattern for model finalizers, but their numeric values are cleared
before the constraints are applied. Projection uses the same callbacks.

The composition, model descriptors, callbacks, and contexts are borrowed by
states and must remain alive while those states are used. A model finalizer may
cache sparse entry indices; a context that caches indices must not be shared
between different symbolic plans. States created with `chem_state_create_like`
share their original symbolic plan, while owning separate numeric workspaces.

To add a Rosenbrock method, populate a new `ChemRosenbrock` parameter set. To
change the controller or linear solver, install the respective callback. To add
a rate type, add a `ChemRateKind`, implement its calculation and validation in
`src/rates.c`, and teach the normalizer and generator to emit its data. A
condition-dependent C function may instead be supplied as a callback rate.

## Tuning controls

- Group size: `--group-size G` or `CHEM_GROUP_SIZE`, default 64. A final group
  may contain fewer cells; every positive group size is supported.
- Threads: `--threads N` or the OpenMP `OMP_NUM_THREADS` setting; the default
  is the OpenMP runtime's maximum available thread count.
- Scheduling: `OMP_SCHEDULE`, default `static`; OpenMP runtime schedules such
  as `dynamic,1` are accepted.
- SIMD width: chosen by the compiler for `omp simd`; there is no fixed or
  machine-specific vector length.

The case runner does not classify remaining reference cases as passing until
their checks exist.
