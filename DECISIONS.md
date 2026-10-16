# Decisions

- The public C interface uses owned model and state handles, explicit error codes,
  and callback tables in place of C++ builders, exceptions, and policy templates.
  All numerical storage is double precision. Mechanisms built in C and mechanisms
  read from configuration use the same model finalization and solver paths.
- A group stores each species contiguously across cells. OpenMP distributes groups;
  `omp simd` expresses independent cell arithmetic. Each group has independent
  adaptive step control and private work arrays.
- Configuration parsing and optional C generation happen only at build time, in
  the local conda environment. The native library has no Python dependency.
- The original prompt stays an untracked input because it contains experiment
  identities. Generated and authored deliverables contain no participant identity.

## Numerical and execution choices

- The port retains the five Rosenbrock tables, backward Euler, weighted error
  control, constraint initialization, Markowitz ordering, and Mozart/Doolittle
  algorithms. Configuration generation emits model construction calls; it does
  not select a different numerical implementation for TS1. All numerical work
  uses double precision and ordinary C/OpenMP.
- Each state is an independent group, following the grouping permitted by the
  benchmark. The group defaults to 32 cells, measured against 16, 64, 128, and
  256. Group size, threads, and scheduling are exposed and documented in
  DESIGN.md. There is no dependence on a seed or a particular cell count.
- The benchmark uses dynamic group scheduling by default, with the standard
  `OMP_SCHEDULE` override. A group's entire history stays together; scheduling
  does not change its arithmetic. `ch_solve_many`, the small public batch API,
  retains static scheduling. Both call the same `ch_solve` implementation.
- Rates are evaluated once before TS1 warm-up because temperature, pressure,
  density, and custom parameters do not change during the run. The supplied
  MICM reference also updates parameters once before warm-up. Every subsequent
  forcing and Jacobian evaluation uses the current concentrations.
- LU values use a conventional unit-lower/actual-upper representation. The
  integration workspace caches the same pivot reciprocals that factorization
  already computes, avoiding repeated divisions for stage solves. The public
  standalone LU interface retains the conventional representation. Reconstructed
  factors, both algorithms, alternative callbacks, analytical solves, and TS1
  are checked rather than depending on a particular internal representation.
- API failures use error codes. Backend callback contexts stay caller-owned;
  structural arrays are copied. The documentation requires a new builder after
  failed finalization, which avoids promising rollback of partial symbolic work.
- The native code contains no operating-system-specific API. Platform-specific
  SDK and OpenMP paths are confined to the Makefile. The `.so` library filename
  is shared across platforms to keep the ctypes harness simple; macOS links it
  with `-dynamiclib`.

## Configuration schema and TS1 extensions

The reader implements the supplied tutorial format rather than linking a legacy
parser. `camp-files` are relative to their containing file, and `camp-data`
accepts `CHEM_SPEC` and `MECHANISM` records. JSON and YAML have the same semantics.
The following explicit extensions permit a lossless TS1 translation and make
the schema usable for the other rate types:

- A mechanism's `rate units` can be `mol m-3`; the default is `molecule cm-3`
  as in the tutorials. TS1 C++ already contains native mol/m3 constants, so
  applying the tutorial unit conversion to it would convert twice.
- Reactants can be an ordered list of species names, with duplicates; products
  can be an ordered list of `{species, yield}` mappings. The existing object
  form remains supported, including reactant `qty`. Lists preserve repeated
  terms and summation order from TS1 rather than merging duplicate products.
- A species can specify `parameterization: air density`. It remains part of the
  mechanism's species metadata, but its value is supplied by the conditions
  rather than integrated. This represents TS1's third-body species `M`, giving
  210 declared species and 209 evolving variables.
- A top-level `initial state` mapping contains `temperature`, `pressure`,
  `concentrations`, and `parameters`. Name/value mappings supply defaults and
  permit scalar broadcasting to all cells. CSV values override these defaults;
  air density is then recomputed. The translator retains all source defaults,
  including the six unphysical photolysis values identified by the prompt; the
  supplied input CSVs override them with zero.
- `USER_DEFINED` explicitly represents a labeled scaled custom rate. Its fields
  are `MUSICA name` and optional `scaling factor`, matching the tutorial label
  convention. Surface labels expand to adjacent `.effective radius [m]` and
  `.particle number concentration [# m-3]` parameters. Species molecular weight
  and diffusion properties provide the surface formula's other inputs.
- Ternary activation, JPL19, Taylor, and reversible formulas have explicit type
  names and fields in `tools/configure.py`. Taylor adds `coefficients` (ascending
  powers of temperature). A branch can specify `branch: 0` or `1`; otherwise
  the two tutorial product lists expand to both channels. Runtime C function
  pointers cannot be represented in a data file and remain a native API feature.

`tools/translate_ts1.py` is a deliberately narrow parser of the supplied C++ data
declarations, not a C++ interpreter. It checks the expected declaration counts
and fails on unsupported rate classes. It extracts all numerical input values;
there is no handwritten TS1 numerical kernel. Its output passes through the
ordinary configuration reader, which is also exercised by all three independent
tutorial configuration comparisons.

## Mapping the CPU cases to C

Every supplied ID is registered explicitly in `tests/cases.txt`; missing
registrations fail. No cases are marked N/A. Several policy instantiations call
a shared C fixture with the corresponding storage/order/LU parameters. The
port checks numerical values, structure, ownership, bounds, and status codes
through the public APIs rather than trying to reproduce template syntax.

| C++ concept | C equivalent exercised |
| --- | --- |
| Copy/move constructors and assignment | Explicit clone, copy, move, owned buffers, independence after copying |
| Exceptions and compile-time misuse checks | Error codes, NULL constructors, dimensions/type/category validation |
| Policy and vector-width templates | Runtime ordering/group-width selection and short final groups |
| Const views, mutable views, temporary views | Read-only flags, bounded views, strict padded operations |
| Capturing lambdas and reusable functions | Function pointer plus context; reusable validated callback descriptors |
| Typed species properties and phases | Owned typed property maps and qualified phase names |
| Device copy calls in CPU policies | Direct host storage and copy/no-op-equivalent visibility checks |
| String matrices in the printing tutorial | Integer intern IDs with the original strings in a constant table |

Randomized tests use fixed NumPy seeds and the same matrix sizes, distributions,
and algebraic checks. C++ standard-library RNG sequences are not portable.
Deterministic numerical fixtures and analytical error bounds are retained.
External aerosol and terminator process fixtures are expressed through the
generic callback interface. The Chapman regression oracle independently uses
explicit equations and dense NumPy solves; it does not call the C implementation
to obtain its expected solution. The full wall-loss history is checked at the
original bound, extending the reference helper's ten-row comparison.

The eight plain-program cases execute native C ports and compare their stdout
with captured reference output. Timing values and concurrent thread-start order
are inherently nondeterministic: timing lines must have the expected format and
count, and startup thread IDs must match. All numerical and other deterministic
output must match exactly. `vectorized_matrix_solver` initializes its reference
solver result as already converged, so its integration loop does not run; the
C tutorial preserves that observable output. Actual grouped integration is
covered by the numerical tests, OpenMP tests, and TS1 gate.

OpenMP checks run in standalone C processes. This avoids loading the conda
numerical stack's OpenMP runtime together with the native compiler's runtime
inside one Python process. There is no duplicate-runtime suppression setting.
The native solver, controller replacement, and linear-solver replacement checks
run with the same compiler and optimization flags as the benchmark.
