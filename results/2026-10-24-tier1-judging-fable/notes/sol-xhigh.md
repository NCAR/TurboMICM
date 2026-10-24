# sol-xhigh: Judge 2's reading notes

Final commit 171f54a. Worktree ~/NSF-NCAR/tier1-work/verify/sol-xhigh/clang/sol-xhigh.

## Shape
- `ck_` API (`include/chem.h`, 319 lines, one statement per line). One generic path: `ck_model_create`
  copies the description, maps species to variables (parameterized species excluded, `src/model.c:297-312`),
  reorders variables by diagonal Markowitz at run time (`src/model.c:313-380`, `ck_markowitz_order`
  `src/model.c:598-643`), stably reorders reactions by rate-type priority (`src/model.c:239-243,381-392`),
  builds pattern, fill, maps, divide/update schedules (`symbolic`, `src/model.c:88-188`) and derivative jobs
  with precomputed targets (`derivatives`, `src/model.c:189-238`).
- TS1 path: `src/bench.c:293-294` `ck_rosenbrock_parameters_init(3, &parameters)` (return unchecked) ->
  `ck_solve` (`src/solver.c:189-201`) -> `ck_solve_raw` (`src/solver.c:44-176`, 133 lines) ->
  `ck_forcing`/`ck_jacobian` (`src/process.c:41-54,93-107`), `ck_dispatch_factor`/`ck_dispatch_solve`
  (`src/process.c:285-300`) -> `ck_factor` Mozart or Doolittle by enum (`src/process.c:108-180`),
  `ck_linear_solve` (`src/process.c:181-209`).
- Layout species-major `y[variable * cells + cell]` (`include/chem.h:208-209`); group 32 default
  (`src/bench.c:11-13`), `CK_GROUP_SIZE` env, `GROUP_SIZE` make variable.
- Benchmark: warm-up, timed loop, `schedule(runtime)` set to static unless `OMP_SCHEDULE` (`src/bench.c:252-253,295-319`).

## Readability
- Allman braces, one statement per line; consistent `ck_` names; rate parameters are named fields
  (`A, B, C, D, E, k0_A, ... Fc, N`, `include/chem.h:54-67`), so formulas read as in MICM (`src/rates.c:62-65,76-80`).
- Comments: 1.0 % of lines (common count), the lowest of the Sol/Astra six; the solver has two comments
  in 201 lines (`src/solver.c:1-3,28`).
- `process.c` holds forcing, Jacobian, both LU factorizations, the triangular solves, the separate-LU
  variants, the dispatch, and the matrix-vector product (309 lines): the file name does not say so, and
  DESIGN.md names no file for the LU.
- `ck_model_create_with_order` runs 297 lines (`src/model.c:248-544`) with the same three-line cleanup
  repeated four times (`src/model.c:320-324,333-336,346-349,360-363`). Common count: 11 functions over 60
  lines and 6 over 100, the most of the Sol/Astra six.
- Documents: README (102), DESIGN (172) with a paragraph per exercise, DECISIONS (62), AUDIT (31) mapping the
  prompt's requirements to evidence, WORKLOG (12).

## Maintainability
- Tests in C: `tests/runner.c` reads `tests/cases.txt`, splits each id into suite and test, and asks 23
  `test_*` dispatchers in turn, each returning 0 pass, 1 fail, 2 not implemented, 3 N/A
  (`tests/runner.c:15-89`); an id no dispatcher knows prints `FAIL port not implemented`. Programs go to
  `test_tutorial(id)`. So every PASS line comes from a check selected by that id's name.
- `tools/test_configuration.py` (10 unittest cases on the reader) runs before the cases (`Makefile` test-cpu).
- The three tutorial configurations are generated in the build (`Makefile` generated/tutorial_*.c) and
  compared with hand-built descriptors (README: rates, forcing, Jacobians, trajectories).
- Error handling: statuses everywhere; `ck_allocate` guards overflow (`src/model.c:4-9`); `goto fail`
  with a complete `ck_model_free` (`src/model.c:33-86`); the error norm propagates NaN
  (`src/solver.c:38`) and the integrator returns CK_NAN/CK_INF (`src/solver.c:134-143`).
- Ordering logic lives once, in C at run time; Python keeps the input order (DECISIONS.md).
- Build: static library from per-file objects, `build/configuration.txt` stamp, `-Wall -Wextra -Wpedantic`.
- `ck_rate_prepare` rejects `(unsigned)r->type > CK_LAMBDA` (`src/rates.c:47`): a type appended after
  `CK_LAMBDA` must move this check; `rate_priority` falls back to 11 for unknown types (`src/model.c:242`).

## Exercises
- E1: public `ck_rosenbrock_parameters` (`include/chem.h:275-288`), generic loop to 6 stages
  (`src/solver.c:46,90`), sets as cases of `ck_rosenbrock_parameters_init` (`src/parameters.c:6-197`),
  DESIGN.md says how and names the capacity. Ease 5.
- E2: `step_factor(parameters, error, context)` (`include/chem.h:281`, `src/solver.c:131-132`) returns the
  factor; no h, no acceptance, context shared through the parameters struct across threads
  (`src/bench.c:293-299,314`). DESIGN.md says the accept/reject rules stay in `ck_solve_raw`. Ease 3.
- E3: `ck_linear_operations {factor, solve, context}` through `parameters.linear`
  (`include/chem.h:268-273,280`), dispatched at `src/solver.c:88,116` and in backward Euler; receives
  (model, cells, packed matrix, algorithm, context); `ck_jacobian_index` maps storage. Ease 5.
- E4: C: enum + fields (`include/chem.h:39-67`), defaults (`src/rates.c:24-44`), prepare check (`:47`),
  calculate (`:66-122`); Python: `TYPES` (`tools/configure.py:20-24`), `DEFAULTS` (`:16-19`, field names must
  match the C struct, `:209`), unit scaling (`:140-150`); tests in `tests/test_rates.c` and a config check.
  `CK_LAMBDA` + `ck_model_set_lambda` for C callbacks. Ease 3.
- E5: `tools/configure.py CONFIG generated/x.c --symbol x`, a Makefile rule, a program; the benchmark
  hard-codes `extern const ck_mechanism ts1` (`src/bench.c:14`); the C model reorders at run time and names
  come from the model (`ck_species_name`), so the input order does not matter; generic kernels; tests
  generate the three tutorials. Ease 4.
