# sol-max: notes for the judge (readability and maintainability)

Facts only. Paths are the entry's at `001fb64` unless they start with `include/micm/` (MICM).
Function counts are from `judging/metrics.json`.

## Tracing the TS1 path
- Timed path: `main` (`src/bench.c:141-149`) -> `chem_solve` (`src/rosenbrock.c:114-192`) ->
  `chem_forcing`/`chem_jacobian` (`src/process.c:49-79`) -> generated `ts1_force`/`ts1_jacobian`
  -> 18 + 30 chunk functions; `factor_shift`, `chem_alpha_minus_jacobian` (`src/rosenbrock.c:28-38`),
  `chem_lu_mozart`, `chem_linear_solve` (`src/linear.c:16-39,63-85`), `chem_axpy`
  (`src/state.c:116-118`), error norm and controller (`src/rosenbrock.c:12-27`), plus `chem_solver`,
  `chem_rosenbrock_parameters` (`src/bench.c:143`), `chem_initialize_constraints` (returns at once,
  `src/rosenbrock.c:66`), `add_stats` (`src/bench.c:80-84`): 14 functions in 6 files plus the
  generated one, depth at most 4. Five calls are indirect (`src/rosenbrock.c:37,152,160`;
  `src/process.c:52,63`), the kernels chosen in the loader (`generated/ts1.c:21401-21402`).

## Size, naming, style
- 113 C functions, median 9 lines; over 100 lines: `chem_rosenbrock_parameters` (202, emitted once by
  `tools/translate_parameters.py:7-39`), `chem_model_finalize` (141: unused-species check, sort,
  derived rate fields, labels, pattern, jobs, fill, LU lists, callbacks; `src/model.c:168-308`),
  `main` (102, `src/bench.c:85-186`). Comment lines 1.7 %; many statements per line
  (`src/model.c:45-47`; `src/bench.c:115`); 35 C lines over 120 characters (`src/bench.c:107`: 150).
- Single letters change meaning by file: `s` is the stride in kernels and LU (`src/process.c:5`;
  `src/linear.c:16`) but the state in the solver (`src/rosenbrock.c:114`); `p` is the parameter set
  (`src/rosenbrock.c:118`), a species (`src/state.c:24`), a rate (`src/model.c:193`).
- `ChemRate` is one flat struct for eleven rate types (`include/chem.h:25-42`) with reused fields:
  `B` is k_r for reversible rates (`src/rates.c:68`), `A` the user-defined scaling (`src/rates.c:70`).
- All structs are public; the benchmark, generated code, and tests write fields directly
  (`src/bench.c:127-130`; `generated/ts1.c:21402`; `tests/solver_cases.py:98,105`).

## Logic kept in two places
- Markowitz order and fill in C (`src/reorder.c:6-25`; `src/model.c:262-263`) and in Python
  (`tools/configure.py:118-133,152-157`); the loader compares only the nonzero count (`tools/configure.py:230`).
- The Jacobian pattern is built twice in C (`src/reorder.c:35-52`; `src/model.c:239-261`); only the
  second skips algebraic rows (`src/model.c:247,251`).
- Rate-type order in Python (`tools/configure.py:10-12,101`) and in C (`include/chem.h:21-24`; `src/model.c:186-191`).
- `tests/api.py:17-87,98-99` mirrors every public struct field by field.
- Statistics summed in three places (`src/bench.c:80-84`; `src/advance.c:19-22`; `tests/tutorials.c:139-142`);
  branched and surface constants in finalize and again in `chem_rate_value` (`src/model.c:194-200`; `src/rates.c:49-53,73`).

## Unused code
- `external_force`, `external_jacobian`, `external_context` (`include/chem.h:113-114`) are read
  (`src/process.c:53,64`; `src/model.c:170,261`; `src/reorder.c:52`) but assigned nowhere (`git grep`).

## Error handling
- Builders and setters return `ChemStatus` (`include/chem.h:10-14`; `src/model.c:76-166`;
  `src/state.c:55-84`); the benchmark checks arguments, CSV values, columns, allocations, and output
  (`src/bench.c:61-64,73-79,88-98,105-110,115,121,165,176-177`).
- Void routines turn allocation failure into NaN: `chem_lu_mozart` (`src/linear.c:18`), so the solve
  returns `CHEM_NAN`, not `CHEM_MEMORY` (`src/rosenbrock.c:162`); `chem_calculate_rates` (`src/state.c:91-92`).

## Documents against the code
- `DECISIONS.md:15` says "The benchmark clamps"; the clamp is in the library's `chem_solve`
  (`src/rosenbrock.c:180-189`), for every caller.
- `DESIGN.md:21` says generated indices "are verified against the finalized model"; the code checks
  the nonzero count only (`tools/configure.py:230`); generic and generated kernels are compared in
  tests (`tests/tutorials.c:261-292`), and for TS1 only inside the `README_example` case
  (`tests/program_cases.py:30-32`).
- `DECISIONS.md:5` is still in the future tense from the first commit ("will exercise").

## Build hygiene
- Per-compiler object directories and flag stamps (`Makefile:38,75-86`); generated sources depend on
  generator and flags (`Makefile:63-72`); `.SILENT:` hides commands (`Makefile:2`); fixed Homebrew
  libomp path for Apple clang (`Makefile:23-26`); `-fPIC` everywhere (`Makefile:56-62`).

## Possible defects
- NULL dereference on allocation failure: if `calloc` of `t->names` fails, `table_read` returns 0 with
  `t->columns >= 1` and `t->names == NULL` (`src/bench.c:41-44`), and `main` calls `table_free`
  (`src/bench.c:102`), which indexes `table->names` (`src/bench.c:33`). Certain by reading; reached
  only when `calloc` fails; not run.
- Not MICM's control flow: after a rejection `chem_solve` stops with `CHEM_SMALL_STEP` when
  `time + 0.1*h == time` (`src/rosenbrock.c:176`), a test not in MICM's rejection branch
  (`include/micm/solver/rosenbrock.inl:266-291`). Read in both; TS1 never rejects a step.
- `chem_lu_mozart` mallocs and frees its reciprocal buffer in each factorization (`src/linear.c:17,38`),
  281,700 times in the TS1 run; the rest of the solve path allocates nothing.

## Wiki
- Tier-1-Implementations' at-a-glance gives Python tests as 3,993 lines; the counts given for this page
  are 3,609 of tests and 384 of tooling (3,993 in all). Its sol-max citations that I checked matched.
