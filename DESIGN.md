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
cell-major state arrays in a group through `chem_integrate_group` and
`chem_integrate_group_steps`.

## Solver path

`src/chem.c` builds the reaction dependency graph and sparse Jacobian pattern
once per group call. It applies a Markowitz species permutation to reduce fill,
constructs the in-place Mozart sparse LU structure, then stores mutable state and
work arrays in the group context. Rate constants are cached for the duration of
the call because temperature, pressure, and rate parameters remain fixed during
an ODE solve. The forcing and Jacobian loops traverse cells contiguously; the
cell loops use `#pragma omp simd` where applicable.

The integrator uses MICM's default three-stage Rosenbrock coefficients, error
estimate, step-size controller, and mechanism tolerances. Independent groups
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

- **Add a Rosenbrock parameter set:** add its coefficients and selection to the
  integrator parameter interface, keeping the error estimate and solver stages
  together as a named method.
- **Replace the step-size controller:** isolate the error-to-next-step policy
  from the stage evaluation in `chem_integrate_group_steps`; preserve the
  current tolerance norm or make its replacement explicit in the API.
- **Swap the linear solver:** replace the symbolic pattern/factor/solve helpers
  while retaining the `chem_factor`/`chem_solve` boundary and the cell-batched
  matrix storage.
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
