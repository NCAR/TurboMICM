# astra-high: Judge 2's reading notes

Final commit 4055d7c. Worktree ~/NSF-NCAR/tier1-work/verify/astra-high/clang/astra-high.

## Shape
- Opaque `ChemMechanism` (`src/internal.h:27-40`), public `ChemState` (`include/chem.h:99-106`), opaque
  `ChemSolver` per (mechanism, cells) that owns the sparse plan, the derivative plan and the work arrays
  (`src/solver_internal.h:10-19`, `src/rosenbrock.c:26-81`): the benchmark creates one solver per group
  (`src/ts1_bench.c:133-134`), so each of the 313 groups builds its own copy of the Jacobian plan.
- `chem_mechanism_create` copies species and reactions and sorts the reactions by rate type with a counting
  sort over `CHEM_RATE_TYPE_COUNT` (`src/mechanism.c:55-67`); the enum carries that sentinel
  (`include/chem.h:22-26`), and every range check uses it (`src/mechanism.c:59`, `src/rates.c:47`).
- `ChemRate` is a tagged union of per-type structs with named fields (`include/chem.h:28-58`); the
  formulas read as MICM's (`src/rates.c:64-72`).
- Reordering: `chem_mechanism_reorder` (Markowitz in C, `src/sparse.c:179-200`, `src/mechanism.c:181-208`)
  is called by the benchmark after `ts1_create` (`src/ts1_bench.c:129`); the generated TS1 is 6,459 lines
  of `ChemSpecies`/`ChemReaction` initializers passed to `chem_mechanism_create` (`build/ts1.c`), no kernels.
- Forcing: generic per-reaction loop with a scratch copy (`src/state.c:107-134`); negative sparse Jacobian from
  the solver's derivative plan (`src/rosenbrock.c:82-114`); a dense Jacobian for tests (`src/state.c:181-191`).
- Integrator `chem_solve` (`src/rosenbrock.c:146-225`, 80 lines): generic stage loop to 6 stages,
  `s->controller(error, p, context)` returns the factor (`:199`), `s->linear.factor/solve` with context
  (`:170,189`); the clamp is inside `chem_solve` (`:217-223`).
- Benchmark: `chem_update_rates` is called inside the timed loop before every solve
  (`src/ts1_bench.c:154`), unlike MICM's driver; group size only by `CHEM_GROUP_SIZE` env (`:121-122`),
  `--threads` is the only flag (`:109-113`); `schedule(runtime)` set to static unless `OMP_SCHEDULE` (`:138`).

## Readability
- Compact but structured: 4-space indent, several short statements per line in places
  (`src/rosenbrock.c:9-16,203-207`), `goto fail` cleanup; common count 3.9 % of lines with a comment
  (the highest of the nine), median function 7 lines, 5 functions over 60 and 1 over 100.
- Names match MICM (`new_function`, `reject_last`, `chem_add_mass_diagonal`, `local_order`); the rate
  union makes `k0_A`, `kinf_C`, `Fc`, `N` explicit.
- The Jacobian plan lives in the solver, the forcing loop in `state.c`, the dense Jacobian in `state.c`,
  the sparse one in `rosenbrock.c`: the TS1 path crosses four files for forcing and Jacobian.
- Error norm floors with a conditional (NaN propagates, `src/rosenbrock.c:139`), `!isfinite(error)` returns
  CHEM_NAN or CHEM_INF (`:201`).
- Documents: README 74, DESIGN 193 (file map, a paragraph per exercise, a parameter table), DECISIONS 310,
  RESULTS 104, WORK_LOG 181.

## Maintainability
- Tests in C: `bin/test_cpu` runs 24 `test_*()` groups (`tests/main.c`) that print `<id> PASS|FAIL`
  lines; `tools/run_tests.py` collects them, rejects duplicates and unknown ids, adds the eight tutorial
  comparisons, prints `FAIL case port pending` for any id without a report (`tools/run_tests.py:18-34`),
  and also runs the config-equivalence programs, the reader's unit tests and `tools/test_ts1.py`.
  No N/A: all 680 pass with both compilers.
- The generator runs in the build on the three tutorials and `mechanisms/config_extra.json`
  (`Makefile` build/no_user.c, with_user.c, robertson.c, config_extra.c) and the results are compared
  with hand-built mechanisms (`bin/test_config`, `bin/test_config_extra`).
- Ordering logic once, in C. No generated kernels, so nothing to keep in step with the plan.
- Error handling: statuses with 19 values (`include/chem.h:12-16`), `chem_alloc` overflow guard, `goto fail`
  with complete frees; `chem_lu_factor` returns CHEM_SINGULAR on a zero or non-finite pivot (`src/sparse.c:91`).
- Build: per-source objects into `build/libchem.a`, a compiler-settings stamp, `-Wall -Wextra -Wpedantic`;
  the test binary is one 25-file compile line repeated twice in the Makefile (prerequisites and recipe).
- Validation tooling: `tools/validate_benchmarks.py` (three seeds, group sizes, schedules, the 10,000-cell gate).

## Exercises
- E1: `ChemRosenbrockMethod` enum + `chem_rosenbrock_parameters` switch (`include/chem_solver.h:29-39`,
  `src/rosenbrock_parameters.c`); the public struct is passed to `chem_solve`; generic loop; DESIGN says how.
  Ease 5.
- E2: `chem_solver_set_controller(solver, fn, ctx)` per solver, i.e. per group (`src/rosenbrock.c:13-14`);
  the callback gets error and parameters, not h (`include/chem_solver.h:46`); DESIGN says acceptance and
  bounds stay in the solver. A PI controller needs the signature extended in one module; per-group state
  fits the per-solver context. Ease 4.
- E3: `ChemLinearSolver {factor, solve, context}` per solver through `chem_solver_set_linear_solver`
  (`include/chem_solver.h:48-57`), called at `src/rosenbrock.c:170,189` and in backward Euler and
  constraint initialization; the callback receives the `ChemSparse` plan (order, elements, index) and
  the packed matrix. Ease 5.
- E4: C: a union member and an enum value before the sentinel (`include/chem.h:22-58`), defaults and the
  evaluation branch (`src/rates.c:22-44,73-119`); Python: one `TYPES` entry `(ctype, member, defaults)`
  (`tools/configure.py:16-27`) and a unit branch (`:137-153`); the emitter writes `.rate.data.<member>.<key>`
  generically (`:224`). Tests in `tests/test_rates.c` and a configuration case. No hidden range table.
  `CHEM_CALLBACK` with function and context exists. Ease 4.
- E5: `tools/configure.py config.yaml build/x.c --symbol x`, then `x_create`/`x_initialize`
  (`README.md`), a Makefile rule; the benchmark hard-codes `ts1_create`/`ts1_initialize`
  (`src/ts1_bench.c:12-13`); no kernels; ordering in C only; tests generate four configurations. Ease 4.
