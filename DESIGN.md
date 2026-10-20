# Design

## Build path

`tools/generate.py` reads YAML or JSON with PyYAML, normalizes the configuration,
and emits a C source file with immutable species, reaction, stoichiometry, and
parameter arrays. Make rules compile every generated C file with the selected
`$(CC)` and `$(OPTFLAGS)`. There is no parser or Python dependency at runtime.
The checked-in `config/ts1.json` comes from `tools/translate_ts1.py`, whose sole
mechanism input is the provided TS1 builder header.

The library API in `include/chem.h` takes a `chem_mechanism` descriptor. It can
also evaluate a single rate constant through `chem_evaluate_rate`, or integrate
reaction forcing from precomputed rates through `chem_evaluate_forcing`, or
integrate cell-major state arrays in a group through `chem_integrate_group` and
`chem_integrate_group_steps`. Callers can select a Rosenbrock coefficient set
with `chem_integrate_group_with_parameters`; the public struct stores up to six
stages, row-wise `a` and `c` coefficients, weights, and step controls. The
constructors include two-, three-, four-, four-stage DA, and six-stage DA sets.
The parameter set also selects Doolittle LU or the in-place Mozart LU used by
the TS1 benchmark.
The public `jacobian_verification.h` helpers use a C forcing callback to compute
finite-difference Jacobians and compare them with a sparse pattern; they support
central differences and forward differences at the nonnegative boundary.
Code-built mechanisms may also use `CHEM_RATE_CALLBACK` with a caller-owned
context pointer; function pointers are not representable in the YAML/JSON
configuration format.

For process evaluation, `chem_evaluate_forcing` clears and fills a cell-major
forcing array, while `chem_accumulate_forcing` adds into an existing array and
can skip algebraic species rows. `chem_accumulate_negative_jacobian` adds the
dense `-df/dy` for each cell and accepts the same row mask. Masked species
values still contribute to rates. These additive APIs let C callers preserve
MICM's ProcessSet-style accumulation without an expression-template layer.

`chem_constraints.h` defines borrowed linear and equilibrium constraint
descriptors indexed against the mechanism's integrated species order. A
constraint set validates algebraic rows, state-parameter indices, and species
dependencies, then exposes its algebraic mask and dense sparsity pattern. Its
forcing operation overwrites algebraic rows with residuals; its Jacobian
operation adds `-dG/dy` to a dense cell-major matrix. Equilibrium residuals use
the Van't Hoff constant supplied in each cell's state-parameter row and clamp
negative concentrations to zero. These APIs support independent constraint
construction and numerical checks. The initialization operation projects only
algebraic species with damped Newton steps, measures correction against the
state tolerances, supports multiple cells, and restores the original state if
projection fails. `chem_integrate_group_with_constraints` assembles these
residuals and Jacobians into the Rosenbrock DAE system: algebraic rows have zero
mass diagonal, stage coupling applies only to differential species, and the
error norm uses algebraic step changes. The caller supplies per-cell state
parameters such as temperature-adjusted equilibrium constants. The DAE path
uses dense cell matrices with the selected Mozart or Doolittle factorization;
the unconstrained mechanism path retains its sparse symbolic factorization.
The public dense LU routines in `chem.h` expose the same no-pivot algorithms in
packed in-place form or as separate unit-diagonal `L` and upper-triangular `U`
matrices, with a reusable forward/back substitution function.
The coordinate sparse LU path accepts a caller-supplied fill-complete pattern,
uses no pivoting, and stores reciprocal diagonal pivots for direct back solve.
The sparse matrix-vector helper and LU routines also cover the packed Chapman
reference operations without adding a mechanism-specific implementation.

`chem_state` is an owning C state container for clients that want state
management as well as the lower-level solver arrays. Initialization copies
species and custom-parameter names and allocates cell-major concentrations,
custom parameters, and rate constants plus per-species tolerances. Copy is
deep; move transfers ownership and clears the source. Named setters support
single-cell scalar values, per-cell vectors, bulk named values, and a checked
raw custom-parameter matrix path. C has no overloaded bracket proxy, so the
equivalent operations use named get/set functions and explicit cell indices.
Errors are returned as status codes instead of C++ exceptions.

`chem_species_info`, `chem_phase`, and `chem_system` provide owning code-built
metadata. Species properties use a tagged value type for strings, doubles,
integers, and booleans. Parameterized species store four linear condition
coefficients and are excluded from phase/system state size and names. Phase
species optionally carry diffusion and density values. Copies own their names,
properties, and nested phase entries; initialize destination structs and
`chem_string_list` values with `{0}`, then release them with the matching
destroy function. String property reads return a borrowed pointer valid until
the property is changed or the species is destroyed. Property and allocation
errors use distinct return codes.

## Solver path

`src/chem.c` builds the reaction dependency graph and sparse Jacobian pattern
once per group call. The benchmark's in-place Mozart LU applies a Markowitz
species permutation to reduce fill; the Doolittle path preserves mechanism
order. Both reuse the symbolic pattern and group-local work arrays. Rate
constants are cached for the duration of the call because temperature,
pressure, and rate parameters remain fixed during an ODE solve. The forcing and
Jacobian loops traverse cells contiguously; the cell loops use `#pragma omp simd`
where applicable.

The default integrator uses MICM's three-stage Rosenbrock coefficients, error
estimate, step-size controller, and mechanism tolerances. Callers may pass
another `chem_rosenbrock_parameters` set through the same integration API.
Independent groups
run concurrently with a static OpenMP schedule. The benchmark creates
`min(threads, cells)` contiguous balanced groups, so each active thread owns one
adaptive solver history. The default thread count is OpenMP's maximum; use
`--threads N` to choose another value. There is no independently tunable group
size or fixed SIMD width in this implementation. Allocation and symbolic
construction are per group call, while arrays within the call are group-local.

At each interval boundary, the implementation clips negative values that are
within the species' configured absolute tolerance to zero. This prevents tiny
accepted roundoff undershoots from becoming artificial sources in very fast
reactions at the next interval. Larger negative values are retained.

## Extension points

- **Add a Rosenbrock parameter set:** fill a `chem_rosenbrock_parameters` value
  with the stage count, triangular `a` and `c` arrays, stage flags, output/error
  weights, and `gamma0`. Pass it to `chem_integrate_group_with_parameters`.
- **Replace the step-size controller:** replace `chem_step_size_factor` and
  adjust the h bounds in `chem_integrate_group_steps_with_parameters`; the
  method stages and error norm remain separate.
- **Swap the linear solver:** select `CHEM_LU_DOOLITTLE` or
  `CHEM_LU_MOZART_IN_PLACE` in `chem_rosenbrock_parameters`; a new method can be
  added at the factor/solve boundary while reusing the cell-batched storage.
- **Add a rate type:** extend `chem_rate_kind`, normalization and descriptor
  generation in `tools/generate.py`, then implement the formula in
  `chem_rate_value` and add direct formula and generated-configuration checks.

## Tuning and portability

The supported run-time performance parameter is thread count (`--threads`),
which defaults to the OpenMP runtime maximum. Group count follows thread count
and is capped by the cell count. OpenMP uses static scheduling; SIMD width is
chosen by the compiler. No architecture-specific flags, intrinsics, libraries,
or APIs are used. The Makefile adds only language, warning, include, OpenMP,
and link flags around `OPTFLAGS`.
