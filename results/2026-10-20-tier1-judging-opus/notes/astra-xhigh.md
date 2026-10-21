# astra-xhigh: notes for the judge

Facts only; paths are at `fbd365c` unless marked MICM's. Not for the wiki.

- **Tracing the TS1 path:** The timed path runs through 17 functions in 7 files: `main` (`src/bench.c:144`), `chem_ros_solve`
  (`src/rosenbrock.c:34`), `chem_state_forcing`/`chem_state_jacobian` and `external_evaluate`
  (`src/external.c:148-169`), the generated `forcing`/`jacobian`, `chem_alpha_minus_jacobian`
  (`src/constraints.c:103`), `chem_state_factor`/`chem_state_solve_linear` (`src/state.c:83-94`),
  `chem_lu_factor`/`chem_lu_solve` (`src/linear.c:78,105`), `axpy`, `chem_normalized_error`, and
  `chem_ros_controller` (`src/rosenbrock.c:8,14,28`), then `chem_state_clamp` and `chem_stats_add`
  (`src/state.c:157,190`).
- Two hops sit in files whose names do not suggest them: the forcing/Jacobian dispatch is in
  `src/external.c:166-169`, and the Rosenbrock diagonal shift is in `src/constraints.c:103-112`.
- The test port added wrapper calls, backend hooks, and `if (constrained)` branches to the hot loop
  (`git diff 342a206 fbd365c -- src/rosenbrock.c`; now `src/rosenbrock.c:50-64,68-69,83-87`).
- **Naming and size:** Uniform `chem_` prefix. Single letters change meaning between files: `s` is the mechanism in
  `chem_model_create` (`src/model.c:53`) and the state in `chem_state_create` (`src/state.c:14`).
  LU operation kinds are bare integers 0–3 with no enum (`src/linear.c:57-66,84-100`; `include/chem.h:125`).
- Long functions: `chem_model_create`, 157 lines, does species maps, custom names, pattern, Markowitz,
  constraints, LU plan, jobs, and kernel choice (`src/model.c:53-209`); `chem_ros_solve` 87
  (`src/rosenbrock.c:34-120`); bench `main` 84 (`src/bench.c:144-227`); `emit` 102 (`tools/generate.py:8-109`).
- Dense formatting: several statements per line, and lines over 120 columns, most in `src/external.c`
  (12 lines), `src/process.c` (9), `src/array.c` and `src/model.c` (8 each). `src/factor.c` and
  `src/sparse.c` write `if(` (22 each) where the other files write `if (`.
- **Logic kept in two places (Python and C):** Diagonal Markowitz: `src/linear.c:128-155` and `tools/generate.py:13-33`. Symbolic fill and row-major
  element numbering: `src/linear.c:33-47` and `tools/generate.py:60-75`. Jacobian job order:
  `src/model.c:170-200` and `tools/generate.py:62-69,94`. Forcing/Jacobian arithmetic: `src/model.c:211-252`
  and `tools/generate.py:83-99`. Rate defaults: `src/rates.c:8-24` and `tools/configure.py:126-155`.
  Air density: `src/rates.c:6` and `tools/generate.py:103`.
- Nothing ties the generated literal indices to the run-time plan: `chem_model_create` installs the
  mechanism's kernels whenever there are no constraints (`src/model.c:201-202`), and the generated
  mechanism carries no nnz or permutation (`tools/generate.py:104-108`). For TS1 they agree: my recount
  from `configs/ts1.json` following the C gives 2,381 nonzeros, and the generated Jacobian's largest
  index is 2,380.
- **Dead, development-only, or unreachable code:** The runner's N/A branch (`tests/main.c:20`) never runs: the registry writes a null reason for every case
  (`tools/test_registry.py:21`). `--implemented` is a development filter (`tests/main.c:7-8,17`).
- `tools/import_parameters.py` (a development helper, its line 3) and `tools/capture_reference.py` are not
  called by the Makefile; `--no-kernels` (`tools/generate.py:115`) is documented (`README.md:86`) but unused.
- **Error handling:** C returns `chem_status` codes (`include/chem.h:11-15`); allocations are checked (`src/model.c:71,82,184`;
  `src/state.c:25-28`; `src/linear.c:14-19,31,41,49`), and failed builds free partial state (`src/model.c:207-208`).
- `chem_state_create` ignores `chem_state_set_custom`'s status for mechanism parameters (`src/state.c:35-38`)
  but checks it for defaults (`src/state.c:39-43`). A failed `ts1_final.csv` write exits 1 with no message
  (`src/bench.c:117-118,220`); CSV and solver failures are reported (`src/bench.c:68,83,88,111,196-198,213`).
- `integer()` rejects `-` only as the first character (`src/bench.c:138`); `strtoumax` skips spaces, so a
  quoted `" -1"` parses as `SIZE_MAX`. Low impact; found by reading, not run.
- There is no zero-pivot check in `chem_lu_factor` (`src/linear.c:87`), as in MICM's
  `include/micm/solver/lu_decomposition_mozart_in_place.inl:183`.
- The Python reader raises `ValueError` on unknown keys, types, species, units, and non-finite numbers
  (`tools/configure.py:33-34,44-45,65-66,86-88,92-93,113-117`); the translator hard-codes TS1's 210/547
  (`tools/translate_ts1.py:54`).
- **Documentation against code:** `DESIGN.md:17` ("Scalar operations keep their reference order") and `DECISIONS.md:6-7`: on the TS1 path
  the reactions keep the configuration's order (`src/model.c:83-85`), and MICM type-sorts them (MICM's
  `include/micm/solver/solver_builder.inl:290-351`), so the summation order differs. `DESIGN.md:133-135`
  says that the low-level path keeps the supplied order.
- The README, the DESIGN parameter table, and the benchmark defaults agree (`README.md:41-48`;
  `DESIGN.md:91-99`; `src/bench.c:13-15,146-148,176`).
- **Build hygiene:** `-O3 -DNDEBUG -std=c11 -Wall -Wextra -Wpedantic` (`Makefile:5,27`); a FORCE'd `build/compiler.conf`
  rebuilds objects when `CC` or the flags change (`Makefile:90-93`).
- Machine-specific defaults: `CONDA=$(HOME)/miniforge3/bin/conda` (`Makefile:6`), the Homebrew libomp
  path for clang on macOS (`Makefile:17-18`), and `clang++` hard-coded for `micm-ref` (`Makefile:82`).
- `make test-cpu` needs the conda Python for the registry (`Makefile:7,45-46`), builds the benchmark first
  (`Makefile:42`), and compiles every test source in one command (`Makefile:50-52`); the runner opens
  `tests/cases.txt` and `tests/reference/*.txt` relative to the working directory (`tests/main.c:9`;
  `tests/test_tutorials.c:7`).
- **Possible defects:** Latent: `preserve_order` with generated kernels. `chem_model_create` skips the Markowitz order when the
  flag is set (`src/model.c:122`) but still installs the mechanism's kernels (`src/model.c:201-202`), whose
  indices assume that order (`tools/generate.py:13-34`). This is not reachable from `generate.py` output,
  which never sets the flag (`tools/generate.py:104-108`), or from the tests (`tests/chapman_fixture.h:26-27`);
  external composition clears the kernels (`src/external.c:114-115`). High confidence in the logic, from
  reading; not executed. None found on the TS1 path as built.
- **Records:** The run log says the context was compacted twice (`admin/runs.md:200` on TurboMICM `main`); the
  transcript and the digest show four compactions (21:52, 22:21, 22:50, 23:17 MDT). The page says four.
- Tier-1-Implementations' astra-xhigh citations check out. One wording point: `tests/test_tutorials.c:31-34`
  compares the generated mechanisms with the entry's C ports of the tutorials' hand-built mechanisms
  (`examples/tutorials.c`), not with MICM's code.
