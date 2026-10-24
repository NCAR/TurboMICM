# sol-max: Judge 2's reading notes

Final commit 001fb64. Worktree ~/NSF-NCAR/tier1-work/verify/sol-max/clang/sol-max.

## Shape
- Builder API: `chem_model_species`, `chem_model_reaction`, `chem_model_reorder` (diagonal Markowitz at run
  time, `src/reorder.c:26-69`), `chem_model_finalize` (stable sort of reactions by rate type, custom labels,
  pattern, fill, LU ranges and update triples, derivative jobs, `src/model.c:168-308`, 141 lines).
- `struct ChemModel` and `struct ChemState` are public with every internal field (`include/chem.h:94-122,142-149`);
  the benchmark and the generated loader write fields directly (`src/bench.c:127-130,168`, `tools/configure.py:202,233`).
- Generated code: `tools/configure.py` emits chunked forcing and Jacobian kernels with literal indices
  (`tools/configure.py:160-192`; `generated/ts1.c`: 18 force chunks and 24 Jacobian chunks of 32 blocks,
  21,407 lines) and a loader `chem_load_ts1(int specialized)` that builds the model through the public builder
  and, if `specialized`, installs the kernels in `m->force_kernel`/`m->jacobian_kernel`
  (`tools/configure.py:193-236`; `generated/ts1.c` tail). The benchmark calls `chem_load_ts1(1)` (`src/bench.c:103`).
- The kernel indices come from a Python copy of the C ordering and fill (`markowitz`, `tools/configure.py:118-133`;
  fill `:152-157`); the loader checks only `m->nnz != <count>` (`tools/configure.py:230`). DESIGN.md says the
  indices "are verified against the finalized model": the verification is the nonzero count.
- Reaction order: Python sorts by `KINDS` position (`tools/configure.py:101`), C sorts stably by enum value
  (`src/model.c:186-191`); the generated kernels index rates by the Python order (`k[ri*s+c]`), so the two
  lists must stay parallel (`tools/configure.py:10-12`, `include/chem.h:20-24`).
- Integrator `chem_solve` (`src/rosenbrock.c:114-192`, 79 lines): generic stage loop to 6 stages; a
  `ChemSolver {parameters, controller, controller_context, factor, solve}` (`include/chem.h:189-198`) is built
  per group in the benchmark (`src/bench.c:135,143`).
- Layout `values[variable * stride + cell]` with `stride == cells` (`DESIGN.md`, `include/chem.h:144`).

## Readability
- Dense: several statements per line throughout (`src/rosenbrock.c:9,37,73,165-170`, `src/model.c:34-35,47`);
  common count 1.7 % lines with a comment, median function 9 lines, 3 functions over 100 lines.
- Rate fields are named but overloaded: `B stores k_r` for REVERSIBLE and `A stores the scaling factor` for USER
  (`src/rates.c:68,70`), mirrored by a rename map in the generator (`tools/configure.py:203`).
- One `CHEM_BRANCHED` type with a `branch` flag (`include/chem.h:33`, `src/rates.c:48-58`); JPL19 shares the
  Troe branch with conditionals (`src/rates.c:34-47`).
- Error norm floors with a conditional that propagates NaN, with a comment saying so (`src/rosenbrock.c:24-26`);
  the integrator returns CHEM_NAN/CHEM_INF (`:162-163`).
- Documents: README 41 lines, DESIGN 49 lines of long paragraphs (one per exercise), DECISIONS 16 bullets,
  RESULTS 76. No file map; the prose names functions.

## Maintainability
- Tests: Python/ctypes (`tests/api.py`); `tests/run.py` looks each of the 680 ids up in `CASES`, calls its
  function, prints FAIL for an unregistered id, N/A on a `NotApplicable` exception (`tests/run.py:35-55`);
  ids register through `@case(...)` with full ids, several ids per function where MICM's policy variants
  share one C layout (`tests/solver_cases.py:21`).
- The generator runs in the build on three examples, `configs/validation.yaml`, and TS1 (`Makefile`
  GENERATED_OBJECTS); README: the examples are compared with code-built mechanisms with generic and generated
  kernels, and generated TS1 kernels are exercised.
- Build: per-compiler object directories, static and shared libraries, `.SILENT`, flags stamp,
  `-Wall -Wextra -Wpedantic`, `-fPIC` on all objects.
- Allocation on the solve path: `chem_lu_mozart` mallocs a reciprocal buffer per call and on failure writes NaN
  into the diagonal (`src/linear.c:17-18`); custom constraints malloc per cell (`src/process.c:120`).
- Error handling: statuses, `chem_resize` overflow guard (`src/model.c:12-18`), `chem_model_free` complete;
  `chem_model_finalize` returns CHEM_MEMORY mid-way leaving partial arrays (freed later by `chem_model_free`).

## Exercises
- E1: public `ChemRosenbrock` (`include/chem.h:175-183`), generic loop (`src/rosenbrock.c:115-116,136`),
  sets as cases of `chem_rosenbrock_parameters` (`src/parameters.c`); DESIGN: "fill a ChemRosenbrock table ...
  no integration-loop change is needed". Ease 5.
- E2: `ChemController(error, parameters, context)` (`include/chem.h:189`, `src/rosenbrock.c:160`) returns the
  factor; no h or acceptance passed; the solver struct is per group in the benchmark, so per-group context is
  possible; DESIGN: acceptance/rejection stay in `src/rosenbrock.c`. Ease 3.
- E3: `ChemSolver.factor`/`.solve` `(model, cells, stride, matrix[, rhs])` (`include/chem.h:190-197`), used
  at `src/rosenbrock.c:37,80,98,152`; `chem_lu_doolittle` as the second implementation; DESIGN says how. Ease 5.
- E4: C: enum, fields, `chem_rate` defaults, `chem_rate_value`, the range check in `chem_model_reaction`
  (`src/model.c:92`, documented), per-type preparation in finalize if needed; Python: `KINDS`/`C_KINDS`
  (parallel positional lists whose order is the reaction sort order), `fields`, `rename`, units; tests.
  `CHEM_LAMBDA` + `chem_model_lambda` for C callbacks. Ease 3.
- E5: one command, a Makefile target, `chem_load_NAME(1)`; the benchmark hard-codes `chem_load_ts1(1)`;
  generated kernels for any mechanism; ordering and fill duplicated in Python with an nnz-only check;
  tests generate four other configurations. Ease 3 (duplicated ordering) to 4 (user path); decide with the
  agent's standard.
