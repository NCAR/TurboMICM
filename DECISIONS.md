# Implementation decisions

- Runtime code is ISO C11 with libm and OpenMP. Python is confined to the local
  prefix environment and build/test tooling. The public API uses explicit error
  codes, ownership, and callbacks in place of C++ exceptions and policies.
- Keep the reference reaction/species ordering, numerical expressions, Rosenbrock
  parameter sets, norm, step controller, and unpivoted sparse LU algorithms.
  Numerical alternatives are selected through explicit solver interfaces.
- Configuration extension: `parameterization: air-density` represents a species
  whose concentration is supplied by environmental conditions. `initial state`
  carries default temperature, pressure, concentrations, and custom parameters.
  These extensions preserve TS1's parameterized third body and InitState.
- Generation emits immutable C descriptors for arbitrary parsed mechanisms.
  All forcing, Jacobian and solver kernels are shared, mechanism-generic C.
- Tests must have concrete assertions for their named behavior. Unimplemented
  cases remain failures; N/A is reserved for C++-only language mechanics.
- The default CAMP concentration unit is molecule cm^-3, matching the tutorials;
  the TS1 configuration declares `concentration units: mol m-3` because its C++
  parameters have already been converted. Conversion depends on reaction order.
- Reactants and products additionally accept ordered lists. This preserves repeated
  species and the order of individual yield updates without combining terms,
  which avoids unnecessary differences in floating-point roundoff.
- Use the reference's default diagonal Markowitz species permutation, including
  its stable tie-breaking, before generating sparse LU schedules. The input
  configuration keeps its original species ordering. Default absolute tolerance
  is 1e-3, as set by MICM's builder, and relative tolerance is 1e-6.
- Public solves apply MICM's post-solve nonnegative clamp only to differential
  variables. `ck_solve_raw` exposes the underlying integrator for direct tests.
- Reference random linear-algebra tests use C++ implementation-specific random
  distributions. Their C equivalents use a portable seeded generator with the
  same lognormal distributions and test the same reconstruction quantities and
  bounds. Random numbers are test inputs, never expected solver results.
- The vectorized tutorial's original loop condition never enters integration.
  Preserve that program's printed behavior; numerical integration is exercised
  by the solver tests and the other tutorials.
- Tutorial fixture comparisons omit elapsed timing and OpenMP startup-log
  interleaving, both nondeterministic. All printed numerical states/statistics
  remain checked. Fixtures were collected from the copied reference build.
- State setters require exactly one value per cell, matching the reference's
  scalar/single-cell and vector/multiple-cell behavior; no implicit broadcast.

- C++ matrix-policy instantiations use one runtime C implementation with explicit
  element size, row/group width, CSR/CSC ordering, and writable-view metadata.
  Policy variants that select the same solver layout use the same C numerical
  path; all five Rosenbrock methods and both LU algorithms are exercised.
- Five cases are N/A: four `ViewCategorySFINAE` compile-time category/concept
  cases and `GroupedView.CategoryTagsAreDistinct`. Runtime access, mutability,
  shape errors, iteration and grouped-view semantics have C checks.
- The legacy hard-coded Chapman solver is represented by the generic Chapman
  descriptor with its original three-stage method, tolerances and initial step.
  Regression comparisons retain the original bounds. The solve comparison uses
  the reference engine's frozen mathematical random stream, expressed in C, to
  retain its particular input samples; otherwise loose legacy integration error
  can exceed its comparison bound for different random samples. Forcing and
  Jacobian expectations use independently written explicit equations.
- Normalized error includes every algebraic variable. Its error component is
  the step change, matching the reference's DAE fix; algebraic variables are
  never hidden by the public post-solve clamp.
- The benchmark uses independent groups with static scheduling by default.
  `OMP_SCHEDULE` can select a different OpenMP runtime schedule, and
  `OMP_NUM_THREADS` or `--threads` controls the thread limit. Positive group
  sizes are configurable at build time or runtime; partial groups are unpadded.

## Tier 3

- Preserve the Tier 1 solver and its default path exactly; experiments use a
  separate `ck_integrator` interface and runtime method names. Refactor only CSV
  helpers shared by the two drivers, with a direct Tier 1 regression comparison.
- Begin with Rosenbrock order/stage sweeps and linearly implicit Euler
  extrapolation. This tests the literature's favorable low-accuracy Rosenbrock
  results against another generic stiff family without importing solver code.
- Use a per-cell maximum of RMS local error by default, an explicit absolute
  tolerance far below the scoring threshold, and negativity-aware rejection.
  These are generic choices, independent of TS1 species names or reference data.
- Cache Jacobians only within an external chemistry interval. Photolysis changes
  at its boundary invalidate both Jacobian and LU; warm starts retain only the
  proposed step. Reuse of ordinary Rosenbrock Jacobians is measured as an
  approximation rather than described as a formally high-order W-method.
- Use exact equality of h for LU reuse. A step-size deadband allows holding h,
  making reuse mathematically consistent with the stored Jacobian and step.
- Separate reading, snapshot storage and CSV output from measured step intervals.
  Count all attempts, including failed initial attempts, in the new engine.
  The legacy statistics retain MICM's convention, which omits some initial
  rejections; comparisons must acknowledge this distinction.
- Prototype in standard-library Python. Share the coefficient tables as data,
  but independently implement nonlinear chemistry, Jacobian, dense factor/solve,
  stage recurrences, error control and extrapolation. Test invariant conservation
  independently in C. No speed claims come from Python.
- Initial scan orchestration accidentally used a shared output directory for
  two concurrent scans. Those initial rows are preserved as diagnostic history;
  affected results are superseded by isolated runs. The scan tool now assigns
  a directory per output CSV. Development timings may overlap builds; final
  operating-point timings are repeated with no other task processes running.
