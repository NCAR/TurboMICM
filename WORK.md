# Completion record

This is a work log, not a claim of completion. The original scope is the full
CPU port and every case in the supplied case list.

## Requirements still to implement and verify

- Generic C mechanism/species/phase/state APIs, all rate types, external models,
  algebraic constraints and initialization.
- Dense/vector/sparse layouts and matrix operations, ordering, both LU algorithms
  in separate and in-place forms, and linear solves.
- Rosenbrock formulas and every parameter set, backward Euler, extensible solver
  and controller APIs, statistics and result states.
- JSON/YAML configuration reader, all tutorial configurations compared with
  hand-built mechanisms, retained TS1 translation script and configuration.
- 680 explicit case mappings with real assertions and matching tolerances;
  eight plain programs with reference output comparisons.
- Generic automatic kernel generation exercised in tests and benchmark.
- TS1 C benchmark: CSV defaults/overrides, warm-up, timing, statistics, final CSV,
  all-core default and thread override.
- Clang and GCC builds/tests; TS1 per-species gate; tuning-parameter correctness;
  performance comparison using reference built with identical optimization flags.
- README.md, DESIGN.md, DECISIONS.md, RESULTS.md, environment.yml, final commit.

## Current evidence

- Inspected the supplied prompt, reference README and case inventory: 680 cases,
  including 672 GoogleTest cases and eight programs. Fresh standalone repository.
- The initial inventory has been replaced by the measured checkpoints below.

## Numerical core checkpoint

Implemented the tagged rate API and all formula types; immutable mechanism
model, custom slots, third-body multiplication, SoA cell states; symbolic
Markowitz and both in-place LU plans; all five Rosenbrock parameter sets and the
ODE integration loop; backward Euler; JSON/YAML reader, C mechanism/kernel
emitter, and TS1 transcription script/config; standalone C TS1 benchmark.

Clang passed 32 mapped rate/lambda tests and the supplied TS1 100-cell gate for
threads/groups (1,32), (4,7), and generic evaluator (1,32). Largest relative error
9.7e-13. The test runner explicitly FAILs all not-yet-ported cases. No N/A cases
have been declared. Full requirements above remain binding.

At this checkpoint, constraints/external models, rich species/phase properties,
dense/vector/sparse utility APIs, separate LU storage, most case mappings,
configuration-vs-hand equivalence, plain tutorial programs, and final validation
are still outstanding. Source reads also found that the reference builder sorts
reactions by rate type; the current model retains declaration order, within the
TS1 gate, but tests must account for reference ordering when comparing rate slots.

GCC 15 also passed all 32 mapped cases and the 100-cell TS1 gate (4 threads,
groups of 7; largest relative difference 9.69e-13). Eight plain-program reference
outputs are retained in tests/reference, produced by tools/capture_reference.py
from build/reference-tests, a complete private copy of the supplied build.

## Constraint checkpoint

Added linear, equilibrium, and callback constraints; algebraic-row masks and
sparsity; finite-difference verification; damped Newton initialization with
weighted corrections and rollback; and the DAE mass-matrix Rosenbrock path.
All 54 constraint unit/initialization cases are mapped, for 86 total implemented
cases. Clang and GCC pass these cases. Integration and external-model cases
remain outstanding. `bin/test_cpu --implemented` is a development-only filter;
the default test target continues to report every unported case as a failure.

## Matrix checkpoint

Added native typed dense/vector arrays, CSR/CSC sparse matrices and builders,
raw and bound grouped views, reusable callbacks, padding-aware reductions,
copy/fill/arithmetic, and printing. All 117 dense/vector matrix cases, 190 sparse
matrix cases, 22 grouped-view cases, and five view-category cases are ported.
Clang and GCC pass all 420 implemented cases; GCC reports no diagnostics for
this checkpoint. The full default suite still reports the 260 unported cases as
failures. Linear-algebra policy cases, system/state/process APIs, solver and
integration regressions, configuration equivalence, and tutorial programs remain.

## Linear solver checkpoint

Added all 52 LU/linear-solver policy cases, eight Rosenbrock core cases, two
backward-Euler core cases, and five version cases: 487 implemented cases total.
The array adapters support separate and combined factors with both LU algorithms.
A per-state backend hook is used by both integrators and constraint initialization;
the rejected-step regression verifies separate L/U against the default backend.
The post-change TS1 gate passes all 21,100 values, worst relative difference
9.7e-13. System/state/process, integration, configuration equivalence, and tutorial
coverage remain outstanding; no N/A cases have been declared.

## Species and system checkpoint

Added owned typed species properties, affine parameterizations, phase species
with optional diffusion/density, named phases, and gas systems with permutation
support. Fourteen reference cases exercise deep copy/destruction, property type
errors, third bodies, state sizes, and names. Clang and GCC pass all 501 mapped
cases; 179 cases remain unported.

## State and process checkpoint

All 22 state cases now cover copying, moving, scalar/vector and batch assignment,
const and mutable column access, raw custom-parameter import, and dimension/name
errors. Added owned reaction construction and a system-to-model builder with
stable rate-type sorting and species-name resolution. All four process cases,
11 reaction-rate-store cases, five process-set cases, and seven Jacobian utility
cases pass. The process-set tests include the three original random-system sizes,
kinetic sparsity, exact forcing/Jacobian entries, and algebraic-row masking.
Clang and GCC pass all 550 mapped cases, leaving 130 unported; no N/A cases.

The integration checkpoint passes 606 mapped cases with both Clang and GCC.
Added eight owned-solver builder cases, all 32 analytical cases and their original
method/layout combinations, twelve built-in DAE integration cases, two Terminator
cases, integrated reaction rates, and the eight-thread shared-solver case. The
normal owned-solver entry point reads immutable solver options; explicit parameter
overrides retain their documented mutation semantics. No N/A cases; 74 remain.
