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

## Numerical components

The Rosenbrock implementation has the two-, three-, four-, and differential-
algebraic four- and six-stage parameter sets. `ChemRosenbrock` holds the
coefficients, tolerances, controller bounds, and projection controls. Its
optional step-controller callback supplies the next step-size factor; a null
callback uses the default formula. Backward Euler has its own parameters and
entry point. Both methods use the same forcing and linear-solver interfaces.

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

The matrix, vector, external-model, and remaining reference-case interfaces are
still being ported; the case runner does not classify them as passing until
their checks exist.
