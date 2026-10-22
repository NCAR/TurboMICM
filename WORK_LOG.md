# Completion ledger

The task is complete. This file retains the implementation checkpoints and
final acceptance evidence. `tests/cases.txt` is the authoritative inventory of all 680 cases.
`tools/run_tests.py` emits FAIL for each missing case implementation.

## Implemented and checked

- Local, pinned conda environment; portable Clang/GCC build and compiler-setting
  invalidation; only C11/libm/OpenMP in the library and benchmark.
- All analytic rate types, callback rates, parameterized reactants, owned
  mechanism metadata, state allocation/copy/setters, forcing and dense Jacobian.
- YAML/JSON tutorial reader, generic C builder generation, retained TS1 importer
  and configuration (210 species, 547 reactions, 209 solution variables).
- All three tutorial configurations compared with hand-built mechanisms at 1
  and 7 cells for rates and ten integration steps.
- Sparse symbolic fill; in-place Mozart and Doolittle LU; triangular solves;
  Markowitz permutation; five original Rosenbrock parameter sets; ODE
  Rosenbrock steps, error calculation and controller. TS1 uses this same solver.
- Benchmark CLI, warm-up, static OpenMP groups, statistics, defaults, CSV input
  and final-state output. Group size is configurable at runtime.
- 31 rate cases, 26 in-place LU/linear-solver cases and 2 backward Euler cases
  ported with original or
  tighter numerical tolerances. No N/A declarations yet.
- 100-cell, 100-step, 30-second TS1 gate passes on Clang and GCC. GCC group size
  7 with 4 threads also passes. Small Clang reference was rebuilt with the same
  flags and compared successfully. Full 10,000-cell seed-2026 gate also passes
  all 2,110,000 values with Clang (identical) and GCC (maximum 1.43e-13 relative).
- Backward Euler Newton iteration, convergence test, step reductions and
  status handling are implemented in `src/backward_euler.c`.
- Public solver post-interval concentration clamping is implemented; a retained
  TS1 fixture exercises the negative-intermediate-state case.

- Linear and Van't Hoff equilibrium constraints, algebraic sparsity/row
  replacement, mass-matrix coupling, algebraic step-change error estimates,
  weighted Newton initialization and damped line search with rollback.
- 76 additional cases (135 total) pass with Clang and GCC: 42 constraint
  type/set cases, 12 initialization cases, 8 DAE integration/overshoot cases,
  7 Rosenbrock utility cases and 7 finite-difference Jacobian utility cases.
  Configuration and TS1 interval-clamp regressions still pass on both.

- Out-of-place factorization preserves its input and uses compact combined L/U
  storage with triangular accessors. Both original algorithms are exercised;
  the separate-storage linear-solver tests use MICM's default Doolittle.
- 36 further cases pass with both compilers (171 total): 26 separate-factor
  LU/linear cases, rejected-step alpha regression, two DAE integrations and
  seven process/rate accumulation cases. GCC rebuild has no warnings.
- Full 10,000-cell TS1 comparison rerun after the DAE/process changes: Clang
  identical to MICM at 8,624.21 ms; GCC max 1.43e-13 relative at 10,293.17 ms.
  These are diagnostic single runs, not ranking medians.

- All 11 reaction-rate store tests and all 22 state tests pass on both compilers
  (204 total). State assignment is transactional; moves transfer ownership;
  named column access has mutable and const variants; bulk setters validate
  dimensions. Named setters now require one value per cell, matching MICM;
  explicit fill functions broadcast, including generated configuration inputs.

- All 59 matrix and 58 vector-matrix cases pass on both compilers (321 total).
  The utility API supports real/integer storage, row-major and grouped layouts,
  deep copies, checked row assignment, scalar operations, element callbacks,
  reusable shape-checked functions, borrowed const/mutable views, row temporaries,
  vector arguments, and strict/non-strict reductions. Group sizes 1–4 retain
  original raw offsets, padding, dimensions and numerical assertions.

- Species typed-property maps and affine parameterizations, phase metadata,
  system names/reordering and owning named reactions now pass all 16 associated
  tests. A system builds the existing runtime mechanism with an owned metadata
  snapshot, so parameter contexts survive input destruction. Supplemental checks
  exercise this lifetime and rate evaluation.
- Five version tests, the callback-rate override integration and the eight-thread
  shared-mechanism/independent-workspace test pass on Clang and GCC (344 total).
  Both builds are warning-free; configuration and TS1 regressions pass.

- All 190 sparse utility cases pass on both compilers (534 total), including
  standard/grouped CSR and CSC layouts, builder deduplication, original raw
  offsets, diagonal operations, exact print output, reusable functions, mixed
  dense/sparse/vector operations, all group-width errors, const views, temporaries,
  vector-index captures, fills and copies. No fill entries are inserted into
  these utility matrices; the solver's symbolic LU structure remains separate.
  Group widths 1–4 and the required width-10 combinations are exercised.
  Both compiler builds have no warnings; configuration and TS1 regressions pass.

- All 22 grouped-view and 5 view-category cases pass on both compilers (561
  total). Group views precompute base pointers and preserve read-only access.
  Grouped fill/copy and non-strict iteration include allocated tail padding;
  strict iteration excludes it. Source template-category tests use C enum tags
  and constructor classification, so none require N/A. Dense/sparse mixed group
  paths, original widths 1/2/4/8, partial tails and cross-view copies are covered.
  All dense and sparse utility inventory cases are now ported. Both builds are
  warning-free and configuration/TS1 regressions continue to pass.

- Added generic external model descriptors with owned contexts, added variables,
  state parameters, resolved Jacobian entries and separate process/constraint
  callbacks. All 20 external constraint cases and 18 aerosol cases pass on both
  compilers (599 total). Tests cover conservative/ratio-only systems, composed
  equilibria, temperature parameters, model order, repeated builds, reordering,
  process entries in algebraic rows and finite-difference Jacobians. Aerosol
  tests cover both integrators, both LU algorithms, factor ownership and ROS3/
  backward-Euler cell batches. Both builds remain warning-free.
- Reran the full 10,000-cell seed-2026 TS1 gate after the external model changes:
  Clang 8,669.65 ms and GCC 10,173.56 ms, unchanged step/call counts. Full final
  state comparisons retain exact Clang agreement and GCC maximum relative error
  1.43e-13. The three configuration comparisons and clamp regression also pass.

- Added the two algebraic-error sensitivity regressions (cascade equilibrium
  tolerance sensitivity and ultra-stiff overshoot), both Terminator cases across
  one through four cells and all twenty longitudes, and integrated reaction-rate
  monotonicity across 100 intervals with a switched-off rate after interval 50.
  The original analytical bounds and solver settings are retained. All five pass
  on both compilers, bringing coverage to 604 PASS / 76 pending FAIL / 0 N/A.
  Both builds have no warnings and configuration/TS1 regressions pass. These
  additions change tests only; the full-size numerical gate above remains valid.

- Ported all 13 Chapman unit cases, the Chapman integration case, all 12
  Chapman forcing/Jacobian/solver regression cases, and both flow-tube cases.
  Chapman tests use generic runtime mechanisms and sparse operations, with a
  separate test-only fixed Chapman oracle for the original comparisons. The
  flow-tube data are retained unchanged and checked after all 3600 intervals.
  All 28 additions pass on both compilers: 632 PASS / 48 pending FAIL / 0 N/A.
  No warnings; configuration equivalence and TS1 clamp regression still pass.
  The only production addition is generic unfactored sparse matrix/vector
  multiplication; the existing benchmark kernels are unchanged.

- Ported all 32 analytical integration cases, including the six rate-law chains,
  their stiff reversible variants, surface uptake, Robertson, Oregonator, HIRES
  and E5. All five Rosenbrock sets, backward Euler where present in the source,
  both LU algorithms and both factor ownership modes retain original tolerances.
  Clang and GCC: 664 PASS / 16 pending FAIL / 0 N/A, no compiler warnings.
  Configuration equivalence and TS1 clamp regression pass. Production kernels
  are unchanged by this test-only addition.

- Added the owning solver-plan builder, retained solve options, explicit cell
  limits, unused-species validation and external model finalization callbacks.
  All eight builder cases retain their original counts and parameter checks.
- Ported all eight standalone tutorials; captured reference stdout from the
  copied MICM build and compare every printed numerical result and counter.
  Duration measurements and concurrent announcement order are validated then
  normalized. Both compiler suites now report 680 PASS / 0 FAIL / 0 N/A.
  All supplementary configuration and TS1 clamp checks pass without warnings.

## Final acceptance audit

- Both `make cpu test-cpu bench CC=clang` and the GCC build pass all 680 cases,
  with zero FAIL/N/A entries and no compiler warnings. All eight program outputs
  match the reference, with only documented nondeterministic fields normalized.
- Supplemental configuration validation passes 49 checks. Generated JPL19,
  Taylor and reversible rates match independent formulas. All 379 TS1 initial
  fields match the supplied reference exactly by name. The importer reproduces
  the committed 210-species/547-reaction mechanism exactly.
- Fresh full-size runs pass all 2,110,000 values: Clang is identical to MICM;
  GCC's largest relative difference is 1.43e-13. Seven additional combinations
  per compiler cover seeds 17/73, groups 1/3/7/32/64/129/1000, threads 1/3/4/10
  and static/dynamic/guided schedules. All pass the original gate.
- Final single-run timings are 9,625.25 ms (Clang C), 10,602.78 ms (GCC C),
  and 27,421.51 ms (paired Clang++ MICM). No optimization flags were tuned.
- Source audit verifies C11/stdlib/libm/OpenMP only, no prohibited extensions,
  OS APIs or target flags; generated C is compiled with the same Makefile flags.
  Only the allowed reference executable and copied MICM tests use MICM code.
- README, DESIGN, DECISIONS, RESULTS and the pinned environment specification
  are present. The implementation and its verification are committed in this
  repository. No required work remains.

## Maintenance notes

`src/internal.h` holds internal mechanism/process structs. `chem_mechanism_reorder`
mutates index maps and must be called before states/solvers exist. The mechanism
pattern excludes forced diagonals during Markowitz selection, matching MICM;
the sparse constructor adds the diagonal for numerical factorization.

`src/solver_internal.h` shares the workspace between integrators.
`src/rosenbrock.c` builds derivative records ordered by independent variable then
reaction and evaluates the negative sparse Jacobian. Dense `chem_jacobian`
returns the positive derivative. The solver owns scratch matrices; states own
cell workspaces and must not be used concurrently. Independent states can share
a mechanism. Linear-solver and step-controller callbacks are public interfaces.

Existing generated objects and validation files are under ignored `build/`.
All Python invocations must use `.conda/bin/python`. Read-only input trees are
`../micm`, `../reference`, and (for Tier 2) `../scenarios`; do not inspect other sibling directories. Run MICM
tests only from a copy of its build tree inside this repository. Only the explicit `make micm-ref` and Tier 2 `make micm-scenario-ref`
reference executables link MICM; the C library and production drivers do not.


## Tier 2 completion audit

- The shared CAMP reader feeds the new Python symbolic planner and C emitters.
  Fixed chemistry indices/yields, sparse fill, LU/substitution operations, rate
  formulas and three-stage traversal are generated; runtime inputs remain data.
  Re-running the generator reproduced the built TS1 C byte for byte. Make
  dependency files name all transitive CAMP includes.
- Generic library and benchmark remain operational. Clean `cpu bench` builds
  pass with Clang 21 and GCC 15.3 using only the required optimization flags.
  Both `test-cpu` runs report 680 PASS, zero failures/skips, plus supplemental
  configuration, tutorial and TS1 checks.
- `test-specialized` passes on both compilers with a strict 1e-10 relative bound
  and no absolute allowance. It covers all five Tier 1 file mechanisms and a
  generated mutation, five cell counts, two tolerances, rate cache invalidation,
  factorization, substitution, singular pivots, concentrations and statistics.
  An alternate variable-stride/chunk/fusion build also passes.
- Both C drivers pass the full 10,000-cell/100-step gate at 18 and 1 threads on
  both compilers. Supplied seed-1 data also passes with groups of 7 and 32,
  static/dynamic scheduling and four threads, for both drivers/compilers.
- Every disclosed scenario passes score.py level 10 for snapshots and final
  states with both compilers, correct grouping and required tolerances. An extra
  48-hour case verifies a different season/latitude/start/NOx, interpolation
  between table points and a one-cell tail. All nonfinite/negative gates pass.
- Source/generated-source portability audit found no forbidden constructs.
  C remains C11 plus libm/OpenMP. Python is build/test orchestration only.
- SPECIALIZATION, README, DESIGN, DECISIONS, RESULTS and environment.yml cover
  the generator, all tuning settings, regeneration, experiments, full timings,
  clean-build times and numerical results. Historical Tier 1 data is retained
  and identified separately. Default Clang executables are restored locally.

Authoritative run evidence is retained under ignored `validation/clang/`,
`validation/gcc-15/`, `validation/seasonal/` and `validation/tuning/`.
