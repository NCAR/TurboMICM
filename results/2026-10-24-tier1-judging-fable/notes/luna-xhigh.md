# luna-xhigh: Judge 2's reading notes

Final commit ada55d5. Worktree ~/NSF-NCAR/tier1-work/verify/luna-xhigh/clang/luna-xhigh.

## Shape
- `src/chem.c` (1,744 lines) holds the rate formulas, the five parameter sets, Markowitz, the symbolic
  build, forcing, Jacobian, both LU factorizations, the solves, the error norms, the ODE integrator, the
  dense DAE integrator and the finite-difference utilities; `include/chem.h` is 502 lines.
- `chem_rate.p[24]` positional with an irregular layout: Troe in `p[5..12]`, JPL19 D's at `p[18], p[19]`,
  kint at `p[20..23]`, branched at `p[13..16]`, k_r at `p[17]`, user scaling at `p[1]`
  (`src/chem.c:52-151`); the only key is Python's `PARAM_INDEX` (`tools/generate.py:26-34`).
- Integrator `chem_integrate_group_steps_with_parameters` (`src/chem.c:1061-1267`, 207 lines): per call it
  builds the symbolic structure (Markowitz, fill, Jacobian jobs, LU update lists: `chem_symbolic_build`,
  `:444-620`), allocates every work array (`:1086-1098`), transposes y into species-major storage
  (`:1099-1101`), computes the rates once (`:1113`), runs `interval_count` intervals, transposes back and
  frees. The benchmark passes all 100 timed steps in one call per group (`src/ts1_bench.c:364-365`), so the
  build and allocation run once per group inside the timer, and once more for the warm-up (`:343-344`).
- Species-major `y[s*cells + c]` with `omp simd` over cells in forcing, Jacobian, LU and solves
  (`:665-674,865-877,904-915,989-1007`); the stage axpy loops are plain loops over `state_size` (`:1175-1199`).
- Controller: a fixed function `chem_step_size_factor` (`:1048-1051`), no hook; LU: `lu_method` enum,
  static `chem_factor`/`chem_factor_doolittle`/`chem_solve` called directly (`:1152-1154,1189`).
- Interval-boundary clipping: negatives within the species' absolute tolerance are set to zero, larger
  negatives kept (`:1236-1246`), a documented departure from MICM's clamp (DESIGN.md).
- Error norm floors with `fmax(result, 1e-10)` (`:1045`), so a NaN step passes `!isfinite` (`:1206`) and is
  accepted.
- Benchmark: `--threads` only; groups = threads (`src/ts1_bench.c:328`); `omp_set_dynamic(0)`.

## Readability
- Long names and one statement per line; 2-space indent; common count 2.0 % comments, median 10 lines,
  12 functions over 60 and 7 over 100 (the most of the nine).
- The library prints diagnostics to stderr on failure (`src/chem.c:1127,1132,1157,1207`).
- Documents: README 74, DESIGN 154 (prose; the extension list says to replace `chem_step_size_factor` and
  edit the integrator for a new linear solver), DECISIONS 63, RESULTS 62.

## Maintainability
- Tests: `tools/test_cases.py` maps suites to 23 test programs; many suites pass the test name to the
  program (`test_state`, `test_matrix`, `test_lu`, `test_constraints`, `test_regression`...), while the
  rate suites, the reaction-rate store, the example programs and `solver_configuration` take one program's
  exit status for all their ids (`tools/test_cases.py:74,145-150,166-169`); unported ids print FAIL with a
  reason (`:35-52`); 18 N/A (version suite, SFINAE, some builder and external-model cases, `:148,200-204,
  237-238,287-288,310-318`). The TS1 gate runs at 1 and 4 threads (`:111-112`).
- The generator runs on the two rate-constant tutorials in the build (`Makefile` NO_USER_GEN, USER_GEN);
  no Robertson; the CAMP reader sets every initial concentration to 0 (`tools/generate.py:71`) and converts
  units with fixed exponents for TROE and TERNARY only (`:131-136`).
- Makefile: a compile and a link rule per test program, 24 times (about 200 lines of repetition).
- Error handling: statuses are -1; allocation per call on the solve path (by design); `chem_accumulate_forcing`
  mallocs per call (`src/chem.c:715`).

## Exercises
- E1: public `chem_rosenbrock_parameters` with constructors (`include/chem.h:299-334`), validated
  (`src/chem.c:323-353`), generic loop to 6 stages (`:1166-1191`); DESIGN says how. Ease 5.
- E2: no hook; `chem_step_size_factor` is a static function (`:1048-1051`) and the acceptance logic is in
  the 207-line integrator; DESIGN says to edit both. The DAE integrator repeats the logic (`:1476-1502`).
  Ease 2.
- E3: `lu_method` enum chooses between two static functions inside the integrator (`:1152-1154`); no
  function pointer; the plan struct `chem_symbolic` is private. Ease 2.
- E4: C: enum (`include/chem.h:12-24`), a `chem_rate_value` case with new slots (`src/chem.c:65-151`),
  `chem_evaluate_rate` checks; Python: `RATE_ENUM`, `PARAM_INDEX`, `type_map`, the `fields` dict and the
  unit branch (`tools/generate.py:14-34,74-88,115-136`); a test. Ease 3.
- E5: `tools/generate.py CONFIG out.c`; the single mechanism is reached through `chem_generated_mechanism()`
  (`tools/generate.py:298-301`), so one mechanism per program; `make cpu CONFIG=...` regenerates the
  benchmark's mechanism; the C side reorders at run time, so the input order does not matter; the fixed
  unit exponents and the zero initial concentrations are hazards; tests generate two tutorials. Ease 3.
