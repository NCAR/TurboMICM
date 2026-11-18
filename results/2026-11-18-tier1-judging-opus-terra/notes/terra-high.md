# terra-high (a8a2ea0): Judge 1's reading notes, 2026-11-18

Root: ~/NSF-NCAR/tier1-work/verify/terra-high/clang/terra-high/ (read-only).

## Counts (tools/metrics_terra.py)
- C library and benchmark 4,057 (lib 3,842 in 47 files + bench/ts1_bench.c 215); 226 functions,
  median 8; over 60: 7, over 100: 1 (chem_rosenbrock_solve 134, src/rosenbrock.c:139); comments
  4.3 %, the highest of the twelve (astra-high 3.9 % the highest of the nine); tests 5,467 C (61
  files, tests/ and examples/) + 539 Python (7 files); tooling 385; generated build/ts1.c 23 long
  lines (89,564 bytes).
- 70 library and benchmark lines over 120 characters, longest 192 (src/cpu_solver.c 10,
  src/state.c 7, src/sparse_matrix.c 6, src/rosenbrock.c 6, src/process.c 6).

## Documents
- README.md 84, DESIGN.md 205, DECISIONS.md 78, RESULTS.md 55.
- DESIGN.md per change: a parameter set is a ChemRosenbrockMethod the integrator consumes
  (:25-27); a rate family: enum and formula, then parsing and unit conversion in the generator,
  "p holds parameters in the order used by src/rates.c" (:69-74); the linear solver: "implement
  the factor and solve operations for that plan, add an algorithm choice to ChemLuAlgorithm, and
  connect it at the factor/solve calls in both integrators" (:78-80); the controller: the block in
  chem_rosenbrock_solve, "shared by single-cell and grouped solves" (:186-191). All match the code.
- DESIGN.md:30 "It allocates scratch space per solve" (src/rosenbrock.c:167-175).
- DESIGN.md:123-165 describe matrix, sparse-matrix and view utilities that "the current solver"
  does not use (:127-128,143-144).

## API and integrator
- include/rosenbrock.h: ChemRosenbrockMethod (stages, a/c, m/e, gamma0, error_order,
  new_function; :10-14); ChemRosenbrockOptions (tolerances, h limits, factor bounds, safety,
  rejection factor, max steps, sparse_plan or NULL for dense, lu_algorithm, constraints,
  factor_out_of_place; :15-29); ChemSolveStatus with NAN and INFINITY apart (:35-40).
- chem_rosenbrock_solve (src/rosenbrock.c:139-272): forcing and Jacobian are callbacks
  (ChemOdeFunction); validation (:147-163); calloc of algebraic mask and n*n+10n work per call
  (:167-175); Jacobian callback each attempt, counted once per step (:193-198); eight factor calls
  chosen inline (:208-216), two solve calls (:234-235); general stage loop with reuse copied
  forward before c-terms (:218-237); NaN and Inf statuses (:248-249); controller block (:250-262);
  clamp at the end (:268).
- Grouped Rosenbrock is a wrapper building a block system for chem_rosenbrock_solve
  (src/rosenbrock_group.c:58-132; :107-126): the controller is written once.
- Backward Euler has its own factor/solve dispatch (src/backward_euler.c:64-70).
- Sparse plan (src/sparse.c:39-110): CSR and CSC lists with diagonal positions; factor returns 1
  on success, 0 on a zero or non-finite pivot (:90-110), the reverse of the 0-is-OK status codes.

## Rates, forcing, Jacobian
- include/rates.h:5-15: ChemRateType, ChemRate {type, p[21], custom_index}; src/rates.c:11-67 one
  evaluator (`return NAN` for an unhandled type), wrapped by batch functions (:70-113).
- src/chem.c: chem_rhs / chem_accumulate_jacobian (:4-49), used by the benchmark's callbacks
  (bench/ts1_bench.c:111-118), and a second, masked and additive pair chem_add_forcing /
  chem_subtract_jacobian (:51-91).
- src/rate_store.c:8-24 maps each type to one of CHEM_RATE_STORE_GROUPS = 11 groups
  (include/rate_store.h:17).

## Generator and reader
- tools/config.py:17-35: reads camp-files; exactly one MECHANISM entry (:34-35).
- tools/generate.py: three orderings — markowitz_order/markowitz_pattern (--reorder, :10-38),
  reference_markowitz_order (--reference-order, MICM's unsigned cost and stable type grouping,
  :41-77), and the C chem_permutation_init (src/reorder.c:11-35) used by the runtime and the tests.
  TS1_ORDER_FLAGS ?= --reference-order (Makefile:4).
- Unit rule per reaction from `rate coefficient units` (SI or molecules cm-3) and the reactant
  count (:143-194); kinds by name `CHEM_RATE_%s` (:233); fixed symbols (:211-234).

## Benchmark (bench/ts1_bench.c, 215 lines)
- CSV columns by name, unknown and duplicate rejected (:58-83); rhs/jac callbacks (:111-118);
  one sparse plan shared (:166-169); rates in the warm-up loop (:178); timed loop (:181-190).

## Tests
- tools/test_cpu.py: CASES table id -> (program, reason) (:9-293); a program's exit status for
  every id mapped to it (:317-323); N/A five (:295-301); unmapped ids FAIL (:313-314). The comment
  at :8 says "Each listed executable must check every assertion in each mapped case."
- Robertson and the two tutorial configurations are generated from MICM's own tree
  (Makefile:73-74,147-150).
- The audit of 2026-11-17 (results/2026-11-17-terra-high-verification): 36 vector-ordered ids rest
  on checks that build no vector-ordered layout; 205 ids not compared with MICM's tests.
- Sample: 8 of 12 faithful; case 3 (all 32 LU ids) a stand-in.
