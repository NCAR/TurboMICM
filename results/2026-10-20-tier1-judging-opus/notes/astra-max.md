# astra-max: notes for the judge

Facts only, from the final commit `0dc7f67`, read-only. Paths are relative to the entry root.

## Tracing the TS1 path

- Setup: `main` (`tools/ts1_bench.c:112`) → `ch_build_ts1` (generated; 1,147 builder calls: 210 species, 222 species properties, 168 parameters, 547 reactions; `build/generated/ts1.c:3-1157`) → `ch_model_finalize` (`src/model.c:136-260`) → `ch_markowitz`, `ch_sparse_create` (`src/sparse.c:7,34`). Per group: `ch_state_create` (`src/state.c:4`), `ch_init_ts1` → `ch_state_set` → `set_column` → `ch_model_species_index`, a linear name search (`src/model.c:267-271`), 377 times per group, before the timer.
- Timed path: `ch_solve` (`src/rosenbrock.c:34-119`) plus 11 callees in 4 files: `ch_force`, `ch_jacobian`, `ch_alpha_minus_jacobian`, `ch_factor_state`, `ch_linear_state` (`src/state.c`), `ch_lu_factor_work`, `ch_lu_solve_work` (`src/sparse.c`), `axpy`, `ch_normalized_error`, the controller through a function pointer (`src/rosenbrock.c:94`), `ch_post_solve` (`src/constraints.c:212`). Each factor and solve call first tests for a replacement callback (`src/state.c:77,81`).
- No mechanism-specific code: the same functions run for TS1 and for the hand-built test mechanisms (`tools/configure.py:184-217`; `DESIGN.md:119-122`).

## Size and form (judge's metrics, `metrics.json`)

- 182 C functions, median 7 lines; 4 over 60 lines, 2 over 100: `ch_rosenbrock_default` 175 (a coefficient table, `src/parameters.c:7`), `ch_model_finalize` 125 (`src/model.c:136`), benchmark `main` 96 (`tools/ts1_bench.c:112`), `ch_solve` 86 (`src/rosenbrock.c:34`). 2.3 % of lines carry comments.
- Dense statement style: several statements per line are common (`src/rosenbrock.c:99-103`; `src/state.c:12-18`); 23 lines of C over 120 columns, longest 148.
- `ch_model_finalize` does six jobs in one function: numbering, first pattern, ordering, renumbering and second pattern, sparse plan, derivative list (`src/model.c:139-255`).

## Naming

- `nr` is the reaction count in `ch_model` (`src/internal.h:52`) and the reactant count in `ch_reaction` and `ch_derivative` (`src/internal.h:34,47`); `ch_model_finalize` also declares a local `nr` (`src/model.c:192`).
- `ch_model_species_count` returns the 209 state variables, not the 210 declared species (`src/model.c:262`); the declared count is `ch_model_builder_species_count` (`src/species.c:138`).
- The post-solve clamp for every solve, `ch_post_solve`, is in `src/constraints.c:212-219`.
- `tools/tutorials.c` is a test program (the eight tutorial cases) kept under `tools/` (`Makefile:84-86`).

## Logic kept in two or three places

- Rate type numbers: the generator emits `.type=<integer>` from the position in `RATE_NAMES` (`tools/configure.py:17-19,204`), which must match the C enum order (`include/chem.h:28-31`); nothing checks the two agree.
- Rate parameter slots and defaults: Python `FIELDS` (`tools/configure.py:20-33`), C `ch_rate_default` and `ch_rate_calculate` (`src/rates.c:10-27,40-84`), and the DESIGN.md table (`DESIGN.md:194-206`).
- C declarations repeated by hand for ctypes: struct layouts (`tests/ffi.py:13-50`; `tests/arrays.py:10-11`) and 116 `bind(...)` prototypes (`tests/ffi.py`, 86; `tests/arrays.py`, 30), with a second `bind` wrapper (`tests/arrays.py:6-7`).
- Not duplicated: the species order is computed only in C (`src/sparse.c:7-32`); the generator writes builder indices (`tools/configure.py:188,201-202`).

## Dead or unused

- `CH_SINGULAR` (`include/chem.h:19`) has a message (`src/model.c:24`) but no code returns it; a zero pivot gives inf or NaN, reported later as `CH_NAN`/`CH_INF` (`src/rosenbrock.c:96-97`). MICM's in-place Mozart has no numerical check either.
- `na()` in the test harness is never called (`tests/harness.py:24-26`).
- In the state path, `ch_lu_factor_work`'s `inverse` argument is replaced by `inverse_diagonal` (`src/sparse.c:137,153`), so the scratch row passed for it is unused there (`src/state.c:78`).
- `tools/import_parameters.py` is a maintainer tool, not used by the build (`tools/import_parameters.py:1-4`).

## Error handling

- Status codes and NULL returns throughout (`DESIGN.md:7-9`); allocation results are checked (`src/state.c:19-22`; `src/model.c:124,145,183,229,242`; `tools/ts1_bench.c:141,147`).
- `ch_update_rates` allocates and frees a parameter buffer on every call (`src/state.c:135-147`); outside the benchmark's timer, but a host that refreshes rates every step would allocate each time.
- A failed `ch_model_finalize` leaves the model renumbered and its reactions rewritten (`src/model.c:184-208` run before allocations at `:227-242`); DESIGN.md says to free it and build a new one (`DESIGN.md:16`).
- The benchmark rejects a cell with nonpositive temperature or pressure without a message (`tools/ts1_bench.c:103`); other input errors print one (`:63,76,84,92,101`).
- `ch_solve` clamps negatives even after an error return (`src/rosenbrock.c:115-117`), as MICM's `Solve` does (MICM `include/micm/solver/solver.hpp:93-95`).

## Documentation against the code

- Consistent on: group-size precedence (`DESIGN.md:141`; `tools/ts1_bench.c:114-131`), scheduling (`DESIGN.md:143`; `:137`), step-size defaults (`DESIGN.md:76-77`; `src/rosenbrock.c:40-42`), tolerances (`DESIGN.md:77-78`; `src/state.c:10`; `tools/configure.py:193`), no allocation in the solve loop (`DESIGN.md:43`).
- DESIGN.md says `src/model.c` "computes LU fill" (`DESIGN.md:38-41`); the fill and the operation list are in `src/sparse.c:34-86`, called from `src/model.c:227`.
- DESIGN.md says small inner loops use `restrict` and OpenMP SIMD (`DESIGN.md:35`); the Doolittle branch has `omp simd` without `restrict` (`src/sparse.c:157-178`).
- DECISIONS.md says every ID "is registered explicitly in `tests/cases.txt`" (`DECISIONS.md:93-94`); `cases.txt` is the supplied list, and the registrations are decorators and two loops in `tests/cases_*.py` (`tests/harness.py:4-21`; `tests/cases_matrix.py:417-424`; `tests/cases_sparse_arrays.py:255`).

## Tests

- Each id runs the function registered for it (`tests/run_cases.py:26-28`); ids sharing a function share a check (`tests/cases_matrix.py:19-20`). No case runs TS1 (no `ts1` in `tests/`); the gate runs were the agent's own, by hand (`RESULTS.md:43-48`).

## Build hygiene

- A flags file rebuilds objects when `CC` or flags change (`Makefile:46-49`); objects depend on all three headers (`Makefile:51`); generated C depends on the configuration and the generator (`Makefile:62,76-81`).
- `make cpu` needs the local conda Python, since TS1's C is generated (`Makefile:62`). Apple clang uses Homebrew libomp at a fixed `/opt/homebrew` path (`Makefile:16-19`); `micm-ref` hard-codes `clang++` (`Makefile:95`). The library is named `.so` on macOS too (`Makefile:12,55-57`; `README.md:24-25`). No target or fast-math flags (`Makefile:29`).

## Looks like a defect?

- Markowitz cost wraps around: `(rows - 1) * (cols - 1)` is `size_t` (`src/sparse.c:22`), and 6 of the 209 TS1 species have an empty column (products only), so `cols - 1` wraps. MICM computes the same expression in `size_t` (MICM `include/micm/solver/linear_solver.inl:68`), whose comment says it never underflows. A Python replay of the entry's code on `configs/ts1.json` gives 2,381 entries after fill, the count the wiki gives for MICM's ordering, so on TS1 the two agree. Sure of the arithmetic (replayed); not a departure from MICM.
- Scale: `ch_markowitz` and `ch_sparse_create` use dense n × n maps and an O(n³) pivot search (`src/sparse.c:10,14-29,37,46`); fine at n = 209, quadratic memory for large mechanisms.
- No defect found on the TS1 path.
