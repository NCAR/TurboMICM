# astra-xhigh: Judge 2's reading notes

Final commit fbd365c. Worktree ~/NSF-NCAR/tier1-work/verify/astra-xhigh/clang/astra-xhigh.

## Shape
- Public structs throughout: `chem_mechanism` (immutable description with optional `forcing`/`jacobian`
  kernel pointers, `include/chem.h:107-122`), `chem_model` (`:144-156`), `chem_state` with per-state
  `factor`/`solve_linear` pointers (`:164-175`), `chem_lu_plan` with its map and op list (`:124-131`).
- `chem_model_create` (`src/model.c:53-209`, 156 lines): species maps, custom labels, pattern, Markowitz
  at run time unless `preserve_order` (`:122-137`), constraints, the LU plan, derivative jobs with flat
  indices (`:170-200`), then `m->forcing = s->forcing && !s->n_constraints ? s->forcing : chem_forcing`
  (`:201-202`): generated kernels are installed without any check against the model's ordering or fill.
- LU: `chem_lu_create` emits an explicit operation list for Mozart or Doolittle (`src/linear.c:25-76`), and
  `chem_lu_factor` replays it with a switch on four op kinds (`:78-103`); kinds 1 to 3 have no comment.
- Generated TS1 (`generated/ts1.c`, 16,341 lines): data, one `forcing` kernel (547 blocks) and one
  `jacobian` kernel with literal flat indices (`tools/generate.py:59-100`), from a Python copy of the
  Markowitz ordering (`:13-34`, with a 64-bit mask for the unsigned cost) and of the fill (`:70-75`).
- Integrator `chem_ros_solve(state, dt, params, controller, context)` (`src/rosenbrock.c:34-120`, 87
  lines): generic loop to 6 stages; the controller receives h and the error and returns the next h
  (`include/chem.h:224-226`); result struct with status and stats; the clamp is a separate call,
  `chem_state_clamp`, that the benchmark makes after each solve (`src/bench.c:205`).
- Error norm: custom `maximum` propagates NaN (`src/rosenbrock.c:7,25`); CHEM_NAN/CHEM_INF (`:99-100`).
- Benchmark options: `--threads`, `--group-size`, `--schedule`, `--chunk`, `--generic` (`src/bench.c:225`);
  CSV parser accepts quoted fields (`:38-63`).

## Readability
- Compact style: many `if (...) { a; b; }` lines (`src/rosenbrock.c:38,56-57,101-110`), 4-space indent;
  common count 2.1 % lines with a comment, median function 7, 5 functions over 60 and 2 over 100.
- The rate union has named fields (`include/chem.h:27-55`); `chem_ros_method` names MICM's DA sets
  `CHEM_RODA3`/`CHEM_RODA4` (`:213`).
- `process.c` holds the owned process model, the rate-type order and queries (`src/process.c`), not the
  kernels; the generic kernels are in `model.c` (`src/model.c:211-252`), the rates in `state.c`
  (`src/state.c:167-188`).
- Documents: README 98, DESIGN 174 (extension points per exercise, a parameter table, how to replace the
  integrator), DECISIONS 118, RESULTS 132, WORK 27.

## Maintainability
- Tests in C with a generated registry: `TEST(prefix, Suite, Name)` macros are scanned by
  `tools/test_registry.py`, which rejects an id not in `tests/cases.txt` or registered twice
  (`tools/test_registry.py:8-17`); `tests/main.c` walks the case list, prints `FAIL not yet ported` for an
  id with no function, else runs that id's function (`tests/main.c:13-24`). No N/A; 680 PASS.
- Generated-kernel checks: the build generates robertson, rates, user_rates, chapman and TS1 with kernels
  (`Makefile`); README names "generated-kernel regressions"; `tools/check_ts1.py` runs nine
  thread/group/schedule/evaluator combinations against the 100-cell reference, `--generic` included.
- Logic in two places: Markowitz and fill in Python (`tools/generate.py:13-34,61-75`) and in C
  (`src/linear.c:128-155`, `src/model.c:122-137`, `src/linear.c:33-37`); no fingerprint when the kernels are
  installed (`src/model.c:201-202`). A mismatch would give wrong kernels silently for a new mechanism;
  the shipped TS1 agrees with MICM to 1e-10.
- Error handling: statuses, `chem_model_create_checked` with messages (`src/model.c:266-290`), `goto fail`
  frees; `chem_lu_create` rejects a missing structural diagonal (`src/linear.c:34`).
- Build: per-source objects, `build/compiler.conf` stamp, `-Wall -Wextra -Wpedantic`; `MICM`/`REFERENCE`
  overridable.

## Exercises
- E1: public `chem_ros_parameters` passed to `chem_ros_solve`; sets as cases of
  `chem_ros_parameters_default` (`src/parameters.c`); generic loop; DESIGN: "return another initialized
  structure; the stage loop does not need a new method-specific branch". Ease 5.
- E2: `chem_step_controller(h, error, params, context)` returns the next h, passed per call with its
  context (`include/chem.h:224-226`, `src/rosenbrock.c:97`); the previous error and h can live in the
  caller's per-group context. Acceptance/rejection stay in the loop (documented). Ease 5.
- E3: `chem_state_set_linear_backend(state, factor, solve, context)` (`include/chem.h:199-205`,
  `src/state.c:75-94`), used at `src/rosenbrock.c:69,87` and by backward Euler and constraint init; the
  plan's map, row_start, column and diagonal are public. Ease 5.
- E4: C: enum and union member (`include/chem.h:21-55`), defaults, prepare if needed, calculate
  (`src/rates.c:8-24,26-43,50-100`); the order function returns `(int)type` for a new type
  (`src/process.c:64-69`); Python: a branch in `rate_initializer` (`tools/configure.py:120-158`) and the unit
  lists (`:95-108`); a test. `CHEM_LAMBDA` + `chem_model_set_lambda` exist. Ease 4.
- E5: `tools/generate.py CONFIG generated/x.c --symbol x` (`--no-kernels` for the generic evaluator),
  a Makefile rule, `chem_model_create(&model, &x, CHEM_LU_MOZART)`; the benchmark hard-codes
  `ts1_mechanism` (`src/bench.c:12`); kernels for any kinetic mechanism, with the duplicated ordering and
  no check; tests generate four other configurations. Ease 3.
