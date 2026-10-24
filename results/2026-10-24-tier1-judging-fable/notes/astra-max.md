# astra-max: Judge 2's reading notes

Final commit 0dc7f67. Worktree ~/NSF-NCAR/tier1-work/verify/astra-max/clang/astra-max.

## Shape
- Builder API with opaque `ch_model` and `ch_state` (`include/chem.h:72-104,135-152`; structs in
  `src/internal.h:51-74`): `ch_model_add_species/parameter/reaction`, then `ch_model_finalize(m, reorder,
  ignore_unused)` (`src/model.c:136-260`, 125 lines): pattern, Markowitz in C (`src/sparse.c:7-32`),
  reorder, derivative plan with flat destinations, `ch_sparse_create` with a Mozart op list
  (`src/sparse.c:34-86`; kinds 0, 1, 2 with one comment at `:65`).
- Generated TS1 (`build/generated/ts1.c`, 1,543 lines): `ch_build_ts1()` builder calls and `ch_init_ts1(state)`;
  no kernels (`tools/configure.py:184-217`).
- Rates: `ch_rate.p[32]` positional slots (`include/chem.h:35-41`), the slot table in DESIGN.md; formulas
  read `p[0] * exp(p[2]/t) * pow(t/p[3], p[1])` (`src/rates.c:44-45`); the Python `FIELDS` lists define the
  same order (`tools/configure.py:20-33`) and `type=RATE_NAMES.index(kind)` must match the C enum order
  (`:168`, `include/chem.h:28-31`).
- Integrator `ch_solve` (`src/rosenbrock.c:34-119`, 86 lines): generic loop to 6 stages; the controller is
  a field of the parameters struct, `p->controller(error, p, p->controller_context)` (`:94`,
  `include/chem.h:169-179`); the linear backend per state (`src/state.c:72-83`), with pivot inverses kept
  for the stage solves (`src/sparse.c:137,208-211`); clamp in `ch_post_solve` (`:117`).
- Error norm floors with a conditional (NaN propagates, `src/rosenbrock.c:22`); CH_NAN/CH_INF (`:96-97`).
- Benchmark `tools/ts1_bench.c`: `schedule(dynamic,1)` unless `OMP_SCHEDULE` (`:137`), one shared parameter
  struct (`:153,166`), rates before the warm-up (`:157`).

## Readability
- Dense: 2-space indent, several statements per line (`src/rosenbrock.c:99-103`, `src/state.c:20-23`);
  common count 2.3 % comments, median function 7, 4 over 60 and 2 over 100.
- Positional rate slots: `p[4] == 0.0 ? ... : ...` selects the branch (`src/rates.c:67`), `p[5]`, `p[6]` are
  prepared values (`:32-36`).
- `ch_model_rate_group` depends on a hard-coded `order[]` table by enum position (`src/inspect.c:18-22`).
- Status codes are negative ints (`include/chem.h:17-22`) and builder calls return ids or errors
  (`:75-86`), so `< 0` checks mix the two.
- Documents: README 113, DESIGN 211 (ownership, layout, numerical methods, tuning table, an extension
  section per exercise with the slot table), DECISIONS 132, RESULTS 126, `validation/README.md`.

## Maintainability
- Tests: Python/ctypes (`tests/ffi.py`); `tests/run_cases.py` walks `CASE_LIST` and calls `REGISTRY[id]`
  or prints `FAIL case port not yet implemented`; N/A from a dict (`tests/run_cases.py:18-34`); modules
  register ids with `@case`/`@group`, where a group function dispatches on the name inside
  (`tests/cases_solver.py:15-33`).
- `bin/test_native_cpu extensions` runs before the cases: a C test that installs Doolittle through
  `ch_state_set_linear_solver` and a counting controller, checks the call counts against the statistics,
  and compares the solution with the default backend (`tests/native_cases.c:33-65`). The only entry
  read so far that tests its extension hooks as such.
- The generator runs on the three examples in the build and the results are compared with hand-built
  mechanisms (`Makefile` EXAMPLE_OBJECTS; README).
- Error handling: negative statuses with strings, `ch_allocate` overflow guard, `goto memory_error`;
  `ch_model_finalize` leaves a partly built model on failure (DESIGN: free it and rebuild);
  `ch_update_rates` allocates a per-call buffer (`src/state.c:135`).
- Build: `lib/libchem.so` and the benchmark from the same objects with `-fPIC`, a flags stamp,
  `-Wall -Wextra -Wpedantic`, `GROUP_SIZE` make variable, `CH_GROUP_SIZE` env.

## Exercises
- E1: public `ch_rosenbrock_parameters` passed to `ch_solve`; `ch_rosenbrock_default(method)` cases;
  DESIGN: "Pass the structure directly to ch_solve"; generic loop. Ease 5.
- E2: `parameters.controller(error, params, context)` returns the factor (`include/chem.h:169`); no h; the
  parameters struct is shared by all groups in the benchmark (`tools/ts1_bench.c:153,158,166`); DESIGN warns
  that per-state counters need per-state contexts. Ease 3.
- E3: `ch_state_set_linear_solver(state, factor, solve, context)` per state (`include/chem.h:204-207`,
  `src/state.c:72-83`), factor receives the `ch_sparse` plan, cells, the packed matrix and a scratch array;
  public CSR getters; used by Rosenbrock, backward Euler and constraint init; tested in
  `tests/native_cases.c`. Ease 5.
- E4: C: enum (`include/chem.h:28-31`), slots and defaults (`src/rates.c:10-27`), prepare, calculate
  (`:40-84`), the range check `rate->type > CH_LAMBDA` (`src/model.c:106`), the `type_order` table
  (`src/inspect.c:18-22`); Python: `RATE_NAMES` position, `FIELDS`, the unit branch; tests. DESIGN lists
  these sites. Ease 3.
- E5: `tools/configure.py CONFIG build/generated/x.c --name x`, `ch_build_x()`, `ch_init_x(state)`, a
  Makefile rule; the benchmark hard-codes `ch_build_ts1` (`tools/ts1_bench.c:16-17`); no kernels; ordering
  in C only; tests generate three examples. Ease 4.
