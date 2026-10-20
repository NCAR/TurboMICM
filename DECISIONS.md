# Decisions

- The library uses an immutable mechanism descriptor and a separate mutable
  solver context per cell group. The same descriptor and integration path serve
  the TS1 benchmark, direct rate checks, and configuration examples.
- YAML and JSON are normalized at build time. Generated C uses ordinary
  standard-C arrays, so the C library and benchmark do not parse configuration
  files and do not link Python.
- C-built mechanisms can carry a rate callback and caller-owned context to
  represent condition-dependent functions and mutable captures. Configuration
  files cannot serialize function pointers or context objects.
- TS1's supplied C++ builder is the translation source because the prompt names
  it as the source of truth. The checked-in translation script emits the
  configuration consumed by the public generator. TS1-specific generated data
  is allowed; hand-written TS1 behavior is not.
- The format supports named user-defined rate parameters and a parameterized
  third-body species (`air_density`) that contributes to pressure-dependent
  rates without being an integrated state species. Surface rates also carry a
  surface species reference and use that species' molecular weight and
  diffusion coefficient. These extensions are represented in the normalized
  descriptor schema.
- CAMP tutorial falloff parameters are normalized from their documented
  concentration units into the internal SI rate units. The tutorial YAML and
  JSON descriptors are checked against equivalent hand-built C mechanisms.
- Each benchmark thread receives one balanced contiguous group, capped at one
  group per cell. This keeps adaptive solver histories independent per thread;
  group size is not independently configurable yet.
- The sparse LU uses MICM's in-place Mozart factor/solve approach. A symbolic
  Markowitz permutation is applied to reduce fill while preserving the original
  mechanism order in the public state layout.
- The solver parameter API also exposes Doolittle LU in natural species order
  for MICM CPU tutorial configurations that use the default Doolittle solver.
  TS1 keeps the in-place Mozart default used by its benchmark reference.
- Mutable per-cell data has a public owning `chem_state` option in addition to
  the solver's low-level arrays. Names are copied into the state, copy is deep,
  and move clears the source. C++ bracket proxies are represented by named
  getter/setter functions with explicit cell indices; invalid names and
  dimensions return error codes.
- Code-built species, phase, and system objects use explicit owning C structs.
  Species property values are tagged as string, double, integer, or boolean;
  phase/system state size omits parameterized species, and their names follow
  MICM's qualified phase and unqualified system conventions. C initialization
  and copy operations require zero-initialized destinations, with destroy
  functions for each owner; string property access returns a borrowed value.
- Process forcing and the dense negative Jacobian are additive C operations.
  Their optional species mask skips output rows while leaving masked species
  in the rate law, matching algebraic-variable handling in process sets.
- Linear and equilibrium constraints are borrowed, immutable C descriptors
  indexed by the mechanism's state order. Each equilibrium constraint names a
  state-parameter column for per-cell K_eq, and set operations provide a dense
  sparsity mask, algebraic rows, residual replacement, and negative Jacobian
  accumulation. Set validation rejects duplicate algebraic rows, matching
  MICM's construction check. Algebraic initialization is a separate operation
  that uses damped Newton projection and restores the input state on failure;
  this keeps projection tests independent of full DAE time integration.
- Tiny negative accepted concentrations are set to zero at interval boundaries
  only when their magnitude is within the configured absolute tolerance. This
  addresses roundoff undershoots that otherwise propagate through fast
  reactions; the solver's integration error control remains unchanged.
