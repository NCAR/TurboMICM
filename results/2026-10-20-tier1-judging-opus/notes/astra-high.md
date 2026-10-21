# astra-high: notes for the judge

Facts only, at `4055d7c`; read-only, nothing built or run.

## Tracing the TS1 path
- Setup: `main` (`src/ts1_bench.c:106-186`) → generated `ts1_create` (`build/ts1.c:3-5697`) → `chem_mechanism_create` (`src/mechanism.c:27-119`) → `chem_mechanism_reorder` → `chem_mechanism_pattern`, `chem_markowitz` (`src/mechanism.c:156-208`; `src/sparse.c:179-200`) → per group `chem_state_create`, `ts1_initialize`, `chem_solver_create` → `chem_sparse_create` (`src/state.c:11-37`; `src/rosenbrock.c:26-81`; `src/sparse.c:38-84`).
- Timed path: about 20 functions in 8 files: `chem_update_rates` → `chem_rate_evaluate` (`src/state.c:82-106`; `src/rates.c:73-119`), plus two hooks that loop over empty model and constraint lists (`src/external.c:160-168`; `src/constraints.c:127-139`); `chem_solve` (`src/rosenbrock.c:146-225`) → `chem_initialize_constraints` (returns at `src/constraint_initialization.c:48`), `chem_forcing` → `chem_add_kinetic_forcing`, `chem_constraints_forcing` (`src/state.c:107-145`; `src/constraints.c:149-168`), `chem_solver_negative_jacobian`, `chem_add_mass_diagonal`, `axpy`, `chem_normalized_error` (`src/rosenbrock.c:82-140`), and through function pointers set at `src/rosenbrock.c:76-77` the controller and `default_factor`/`default_solve` → `chem_lu_factor`/`chem_lu_solve` (`src/sparse.c:85-178`).
- The kernels have one code path for all mechanisms; there is no generated arithmetic (`DESIGN.md:192-193`).

## Naming, size, density
- Uniform `chem_` prefix and `Chem` types. Short locals carry different types: in `src/rosenbrock.c` `s` is a `ChemSolver` in `chem_solve` (`:146`) and a `ChemState` in `chem_normalized_error` (`:128`); in `src/state.c` it is a `ChemState` (`:82`).
- 253 functions in the C library and benchmark; 10 exceed 50 lines: `chem_rosenbrock_parameters` 175 (coefficient tables, `src/rosenbrock_parameters.c:8`), `chem_mechanism_create` 93 (`src/mechanism.c:27`), `main` 81 (`src/ts1_bench.c:106`), `chem_mechanism_create_with_models` 81 (`src/external.c:54`), `chem_solve` 80 (`src/rosenbrock.c:146`) (brace-depth count, `judging/funcs.py`).
- Comment lines are 138 of 3,555 non-blank lines (3.9 %). 43 library lines exceed 120 characters, the longest 201 (`src/system.c`); many statements share a line (e.g. `src/rosenbrock.c:204-207`; `src/system.c:67-81`). Tests are denser: `tests/test_sparse_functions.c` has 52 lines over 120 characters, the longest 253.

## Logic in two places
- Partial derivatives are computed twice: the dense `kinetic_jacobian_add` used by tests (`src/state.c:146-191`) and the solver's derivative records (`src/rosenbrock.c:48-75,82-108`).
- The Jacobian sparsity relation is encoded in `chem_mechanism_pattern` (`src/mechanism.c:156-180`) and again in `chem_kinetic_has_element` (`src/mechanism.c:223-234`).
- The Chapman regression has a test-only translation of MICM's hard-coded solver (`tests/chapman_oracle.h`; `DECISIONS.md:210-217`).
- The species order is computed only in C (`src/mechanism.c:181-208`); Python does not repeat it.

## On the TS1 path but inactive, or unused by it
- Constraint and external-model branches run but do nothing for TS1 (`src/rosenbrock.c:109-113,184-187,194-197`; `src/state.c:139-144`); DAE support added a per-term test inside the forcing and Jacobian loops (`src/state.c:121,127`; `src/rosenbrock.c:96,102`).
- `chem_jacobian`, `chem_sparse_multiply`, `chem_kinetic_has_element` are used only by tests (grep of `src/` and `tests/`). The stage buffer is sized for six stages; ROS3 uses three (`src/rosenbrock.c:45`).
- Each group builds its own pattern, fill, map and records (`src/rosenbrock.c:34-75`): 313 copies at 10,000 cells. Estimate from the array sizes: about 1.8 MiB per 32-cell group, 341 KiB of it the 209 × 209 `size_t` map.

## Error handling
- `ChemStatus`, 19 codes, with `chem_status_string` (`include/chem.h:12-17`; `src/rates.c:8-16`); allocations checked through `chem_alloc` (`src/internal.h:50-54`) with `goto fail` cleanup (`src/rosenbrock.c:79-80`; `src/sparse.c:82-83`); overflow checks before sizing (`src/state.c:18-22`; `src/rosenbrock.c:40-41`).
- `chem_lu_factor` returns `CHEM_SINGULAR` on a zero or non-finite pivot (`src/sparse.c:91,121`).
- The benchmark exits 1 with no message when setup fails (`src/ts1_bench.c:129,133-134`) or a CSV value does not parse (`:83-85`); unknown columns and missing rows print a message (`:73,77`).
- `tools/configure.py` raises `ValueError` on malformed input (e.g. `tools/configure.py:29-35,120-121`); `DESIGN.md:19-20` counts 49 such checks in `tools/test_configuration.py`.

## Documents against the code
- `DESIGN.md:82` says `src/sparse.c` "adds symbolic fill once"; it runs once per solver, so once per group (`src/rosenbrock.c:34-36`; `src/ts1_bench.c:130-135`).
- `README.md:42-43` describes the timed steps without the rate recomputation at `src/ts1_bench.c:154`.
- `DECISIONS.md:9-10` still says the layout "will allow independent groups" (future tense).
- `DECISIONS.md:301-302` says the TS1 importer was rerun and compared with the committed JSON; no make target or test runs `tools/translate_ts1.py` (grep of `Makefile`, `tools/`, `tests/`).
- `WORK_LOG.md` is kept beside the required documents.

## Build hygiene
- A stamp file rebuilds everything when `CC` or the flags change (`Makefile:65-72`); every object depends on every header (`Makefile:43`); `-Wall -Wextra -Wpedantic` (`Makefile:25`); `make cpu` also builds the benchmark (`Makefile:42,62`).
- The clang branch on Darwin hard-codes Homebrew's libomp prefix (`Makefile:13-16`); `CONDA` defaults to `~/miniforge3` (`Makefile:6`).

## Tests
- 672 cases are `RUN(id, call)` sites in 24 C modules (`tests/test.h:12`). A few ids share one call: five ConstraintSet ids run `set_basic(2)` (`tests/test_constraints.c:280-282,286`), and pairs or triples share `arrhenius_test`, `taylor_test`, `markowitz_test`, `bad_builder`, `metadata(0)` (`tests/test_rates.c:180-181,194-195`; `tests/test_sparse.c:143-144,153`; `tests/test_sparse_matrix.c:148-149`; `tests/test_external.c:365,369`): 16 ids on 6 shared calls among 299 literal id sites (regex scan; ids built in loops not covered).
- Tutorial references were captured by running MICM's own tutorial binaries from a copied build (`DECISIONS.md:270-273`); the prompt allows running MICM's tests for reference.

## Possible defects
- `load_cells` writes a CSV custom column only to the first parameter with that label (`src/ts1_bench.c:72,90`), where `chem_state_set_custom` and `chem_state_fill_custom` write every parameter with the label (`src/state.c:73-79,216-218`). TS1's 142 user-defined and 13 surface labels are all distinct (read from `mechanisms/ts1.json`), so TS1 is unaffected. Sure of the code path; matters only for a mechanism whose reactions share a label.
- `positive_integer` rejects '-' only as the first character (`src/ts1_bench.c:24`); `strtoull` skips leading spaces and negates, so " -5" passes as a huge count and setup then fails without a message. Certain from C11 7.22.1.4; harmless.
- None found in the numerics: the clang result is bit-identical to MICM's in every timed run.

## Tier-1-Implementations citations
- All 109 line references in its astra-high section resolve to the code they describe (each cited line dumped and read). One is loose: `src/state.c:11` labeled `ChemState` is `chem_state_create`; the type is at `include/chem.h:99-106`.
