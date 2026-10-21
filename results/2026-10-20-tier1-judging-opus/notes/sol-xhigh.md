# sol-xhigh: notes for the judge (readability, maintainability)

Entry paths at `171f54a`; MICM paths under `micm/include/micm/`. Read-only, nothing built or run.
Function sizes from `judging/funcs.py`, as in `metrics.json`.

**Tracing the TS1 path**
- Timed path: `main` (`src/bench.c:215`) → `ck_solve` (`src/solver.c:189`) → `ck_solve_raw`
  (`src/solver.c:44`) → `src/process.c`: 21 functions, 7 in `src/solver.c`, 9 in `src/process.c`,
  3 in `src/constraints.c`, 1 each in `src/external.c` and `src/model.c`. Four are hooks with
  nothing to do for TS1 (`src/constraints.c:34,68,155-156`; `src/external.c:111-120`). Setup adds ~25.
- The plan the kernels replay sits in uncommented, terse fields of `struct ck_model` (`n`, `nr`,
  `nc`, `nnz`, `njobs`, `nops`, `map`, `divides`, `update_start`; `src/internal.h:25-46`).
  Comment lines: 39 of 3,787 non-blank library and benchmark lines.
- Naming: `ck_` prefix on every public name; kernels carry MICM's names (`ck_alpha_minus_jacobian`,
  `ck_normalized_error`); one-letter locals for model, state, workspace, parameters, job
  (`m`, `s`, `w`, `p`, `d`; `src/solver.c:44-57`; `src/process.c:55-63`).

**Large or duplicated functions**
- 133 functions, median 12 lines; 11 over 60, 6 over 100: `ck_model_create_with_order` 297
  (`src/model.c:248-544`: validation, ordering, rate sort, binding, custom names, constraints,
  external models), `ck_rosenbrock_parameters_init` 192 (`src/parameters.c:6-197`), `main` 140
  (`src/bench.c:215-354`), `ck_solve_raw` 133 (`src/solver.c:44-176`), `read_cells` 117
  (`src/bench.c:72-188`), `symbolic` 101 (`src/model.c:88-188`).
- The same three `free` calls appear five times in the ordering block
  (`src/model.c:320-322,333-335,346-348,360-362,374-376`). Near-duplicates: `ck_linear_solve`
  and `ck_linear_solve_separate` (`src/process.c:181-209,238-271`); `ck_state_set` and
  `ck_state_set_custom` (`src/state.c:110-133`); the two Doolittle loops (`src/process.c:140-177`).

**Logic kept in two places**
- The reaction pattern is built twice in C, in input indices for the ordering
  (`src/model.c:326-356`) and permuted in `symbolic` (`src/model.c:98-108`). The species order is
  computed only in C; the generator writes input order (`tools/configure.py:193-200`).
- Rate defaults in Python (`tools/configure.py:16-19,161-164`) and C (`src/rates.c:24-44`) agree;
  the generator writes every field (`tools/configure.py:209`), so the C ones serve C-built
  descriptors. The 1e-3 absolute tolerance is in both (`tools/configure.py:197`; `src/state.c:66`).
- Custom-parameter names are composed in C (`src/model.c:428-434`) and must match the importer's
  keys from MICM's source (`tools/import_ts1.py:54`); a mismatch fails with `CK_UNKNOWN`
  (`src/state.c:77-84`). `rate_priority` is indexed by the `ck_rate_type` enum
  (`src/model.c:239-243`; `include/chem.h:39-53`); DESIGN.md says to update both (`DESIGN.md:93-95`).

**Dead or redundant code**
- `if (s)` is always true after the check at `src/solver.c:191` (`src/solver.c:198`);
  `v[count == 1 ? 0 : c]` matters only when `cells == 1`, as `count != s->cells` is rejected first
  (`src/state.c:112-119,124-131`); the reactant and product pointer check runs twice
  (`src/model.c:268-274,396-397`); `CK_SINGULAR` is never returned (`include/chem.h:23`).
- The public utilities in `src/matrix.c`, `src/verify.c`, and `src/system.c`, and
  `ck_factor_separate`, are called only by the tests and each other, not by the solver or benchmark.

**Error handling**
- Constructors check every allocation and unwind through one `fail:` label (`src/model.c:541-543`;
  `src/state.c:90-92`); `ck_allocate` guards size overflow (`src/model.c:4-9`). Kernels are `void`
  and check nothing (`src/process.c:5,108,181`); a zero pivot is not detected and surfaces as
  `CK_NAN` or `CK_INF` from the error norm (`src/solver.c:134-143`).
- The benchmark checks every status (`src/bench.c:269-273,284-289,301-306,322-330`) and rejects
  bad CSV and arguments (`src/bench.c:103-129,150-174,224-251`), but ignores the return of
  `ck_rosenbrock_parameters_init` (`src/bench.c:294`).

**DESIGN.md and README.md against the code**: none found to disagree. Checked the build, test, and
benchmark usage (`README.md:10-60`), defaults and tuning table (`DESIGN.md:58-60,153-159`), timing
(`DESIGN.md:141-151`), LU and ordering (`DESIGN.md:66-71`), extension recipes (`DESIGN.md:73-99`).

**Build hygiene**
- `-std=c11 -Wall -Wextra -Wpedantic`, `OPTFLAGS ?= -O3 -DNDEBUG`, and a stamp that rebuilds on a
  compiler, flag, or `GROUP_SIZE` change (`Makefile:4,23,81-86`); machine-specific defaults
  (`Makefile:5,12-13`). `results/` is ignored (`.gitignore:9-12`), so no log behind RESULTS.md is kept.
- One-time scripts sit outside Make: `tools/import_testdata.py` reads a fixed relative MICM path
  (`tools/import_testdata.py:8`); `tools/import_parameters.py` emits one-line C
  (`tools/import_parameters.py:10-18`) where `src/parameters.c` is clang-formatted.
- Five test files include `src/internal.h` and read model fields (`tests/test_configuration.c:3,62-72`).

**Possible defects**
- None found on the TS1 path (high confidence: every function on it read and compared with MICM's).
- The clamp maps NaN and -0.0 to +0.0 (`y > 0 ? y : 0`, `src/solver.c:186`); MICM's
  `std::max(y, 0.0)` keeps both (`util/matrix.hpp:378-384`). Accepted states are finite, so TS1 is
  unaffected (certain, by reading).
- Scaling: ordering and fill use dense n × n arrays and `map` is n × n `size_t`
  (`src/model.c:91,132,602`): 349,448 bytes for TS1's 209 variables, O(n³) ordering time. MICM
  replaced this dense form for large mechanisms (`solver/linear_solver.inl:18-26`).
