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
- Generated kernels are derived from arbitrary parsed mechanisms. There are no
  hand-written mechanism-specific runtime kernels.
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
