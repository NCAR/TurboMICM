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
