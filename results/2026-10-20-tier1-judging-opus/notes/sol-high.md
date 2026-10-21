# sol-high: notes for the judge (not for the wiki)

Facts met while writing the entrant page; paths are at `1d944f4`.

## Tracing the TS1 path
- Setup: `main` (`src/bench.c:148`) → `chem_mechanism_create` + `valid_data` (`src/mechanism.c:65,6`) → `chem_state_create` (`src/state.c:6`) → `read_cells` (`src/bench.c:84`) → `chem_calculate_air_density`, `chem_update_rates` → `chem_rate_evaluate` (`src/state.c:150,162`; `src/rates.c:32`).
- Timed path: 13 functions in 5 files: `chem_solve`, `integrate`, `axpy`, `chem_normalized_error` (`src/rosenbrock.c:133,28,9,14`); `chem_initialize_constraints`, returning at once (`src/constraints.c:119-123`); `chem_forcing` → `chem_add_forcing`, `chem_negative_jacobian` → `chem_subtract_jacobian` (`src/process.c:7,12,47,52`); `chem_alpha_minus_jacobian`, `chem_factor`, `chem_linear_solve` (`src/linear.c:202,8,29`); `chem_post_solve_clamp` (`src/state.c:204`). Factor and solve are reached through function pointers in the state (`src/rosenbrock.c:58,82`).
- The 9 C files the TS1 path touches hold 1,642 lines, plus 253 in `include/chem.h` and `src/internal.h`; much of that serves other features (`src/constraints.c` is on the path only for its early return). No macros or generated code beyond the data tables of `build/ts1.c`.
- 5 of the 21 library files and their headers (`array.c`, `sparse.c`, `system.c`, `jacobian.c`, `builder.c`; about 1,055 lines with headers) exist for test ports and tutorial programs; no solver source includes their headers (`#include` lines of `src/*.c`).

## Naming and size
- Prefixes `chem_` and `Chem` throughout `include/chem.h`. Locals are one or two letters (`m`, `s`, `d`, `p`, `r`, `n`, `nc`; e.g. `src/rosenbrock.c:28-34`). Rate parameters are positional, `double p[24]`, with the layout per type only in a comment (`include/chem.h:20-25`); `src/rates.c:27-98` indexes them by number.
- metrics.json: 188 functions, median 6 lines; 6 over 60, 3 over 100: `chem_rosenbrock_parameters` 176 (a data table), `chem_mechanism_create` 105, `integrate` 104; then `chem_builder_build` 80, `main` 79, `chem_rate_evaluate` 67. Comment lines 3.3 %.
- Dense style: several statements per line (`src/state.c:141-142`, `src/mechanism.c:174-176`, `tests/supplemental.py:14,21`); `src/array.c` and `src/sparse.c` often omit spaces after commas (62 of 296 lines in `array.c`).

## Logic kept in two places
- Diagonal Markowitz: Python (`tools/configure.py:167-205`, used for TS1) and C (`src/linear.c:150-180`, called only from `tests/test_linear.py:142`, `tests/test_dae_integration.py:139`, `tests/test_regression.py:17`).
- Kinetic sparsity: `tools/configure.py:176-182`, `src/mechanism.c:79-87`, and again `chem_kinetic_nonzero` (`src/process.c:118-130`).
- Rate defaults: `tools/configure.py:48-60` and `chem_rate_default` (`src/rates.c:100-123`).
- Every public struct is mirrored by hand in ctypes (`tests/chemffi.py:12-76`, also `tests/builderffi.py:10-17`, `tests/systemffi.py:9-13`); a field added to `include/chem.h` needs the same edit there.
- Forward/back substitution loops appear twice (`src/linear.c:29-54` and `:122-148`); the gas constant is written twice (`src/rates.c:7`, `src/state.c:155`).
- Tests pass rate types as integers, `rate(9)`, `rate(10,...)` (63 calls, e.g. `tests/test_configs.py:28-35`); inserting an enum value would shift them.

## Adding a rate type (DESIGN.md vs code)
- `DESIGN.md:144` says "extend the enum and default factory in `rates.c`"; the enum is in `include/chem.h:6-10`. Not mentioned: `valid_data` rejects any type above `CHEM_CALLBACK` (`src/mechanism.c:30`), and `rate_family` indexes a 12-entry table by type without a bound check (`src/mechanism.c:213-221`). A value appended after `CHEM_CALLBACK` fails mechanism creation until line 30 changes; a value inserted before it shifts `CHEM_CALLBACK` and the integers the tests use. The generator emits `CHEM_<kind>` names, so configure.py kind strings must equal the enum names (`tools/configure.py:236`).
- `DESIGN.md:113` says `CHEM_GROUP_SIZE` can be set "at build"; the Makefile has no variable for it, so it must ride in `OPTFLAGS` (`Makefile:5,23`). The generator's `--no-reorder` is likewise not reachable from the Makefile (`Makefile:47`).
- Otherwise README.md and DESIGN.md match the code where checked: groups, scheduling, warm-up, CSV handling, clamp, extension hooks (`DESIGN.md:111-142`; `README.md:51-64`).

## Error handling
- `ChemStatus` codes (`include/chem.h:84-86`) from every create and solve function; allocation failures are checked and cleaned up (`src/mechanism.c:73,76,90,116-117,140,165-168`; `src/state.c:34-37`). Descriptor validation covers indices, duplicate names (O(n²), `src/mechanism.c:11-20`), finite yields, parameter ranges.
- The solver checks every pivot and every new value for finiteness (`src/linear.c:12`; `src/rosenbrock.c:99-101`); rate constants must be finite (`src/state.c:197`).
- Accessors do not check NULL (`src/state.c:88,144-149`). Functions taking `const ChemState *` write the state's scratch row (`src/process.c:16,56`), so two threads must not share one state (each group has its own).
- Benchmark: CSV errors name the row or column (`src/bench.c:96-112`); a failing group stops only its own steps and is reported after the timer (`src/bench.c:200,208-211`); usage errors exit 2 (`src/bench.c:223-225`). Python raises `ConfigurationError` on unknown species, types, or units (`tools/configure.py:111-137`).

## Build hygiene
- Fixed `CFLAGS` with `-Wall -Wextra -Wpedantic` (`Makefile:23`); `build/compiler-flags` with a FORCE rule rebuilds when CC or flags change (`Makefile:52-57`).
- The 14 library sources are compiled separately into each of five outputs (bench, `libchem.so`, `test_rates`, `rate_tutorial`, `cpu_tutorials`; `Makefile:50,64,69,86,91`); there is no static library.
- `test-cpu` prerequisites are spread over five rules (`Makefile:65,71,83,87,92`). Homebrew libomp and `~/miniforge3` paths are hard-coded (`Makefile:7,13-14`).
- `build/ts1.c` depends on `configs/ts1.json` but not on `tools/translate_ts1.py` (`Makefile:45-47`); the JSON is regenerated by hand (`README.md:83-84`) and a test checks retranslation (`tests/test_configs.py:97`). The translator hard-codes TS1's temperature and pressure (`tools/translate_ts1.py:58`; values equal MICM's `test/benchmark/ts1_mechanism.hpp:6501-6502`).

## Test fidelity
- Each of the 680 ids runs its own check: 26 `Suite.Test` names shared by two executables are all registered by full id (`tests/test_analytical.py:113-114`, `tests/test_aerosol.py:172`, `tests/test_chapman.py:141-142`, `tests/test_sparse_arrays.py:130`); checked by reading `tests/cases.txt` and the registrations.
- Standard/vector pairs share one function (`tests/test_solvers.py:24`); analytical variants with the same C representation run once at the tightest bound (`tests/test_analytical.py:98-110`; `DECISIONS.md:74-79`).
- `PYTHON_CASES` merges 21 dicts with `{**a, **b}` (`tests/run_cases.py:28`), so a duplicate key would silently override; a static scan of the literal keys found no duplicates (generated keys not checked).

## Possible defects (none confirmed by running)
- `factor_doolittle` scans the dense index map in a triple loop, about n³/3 ≈ 3 million lookups per factorization for TS1 (`src/linear.c:59-71`): correct, but slow if chosen. Read, not run.
- `chem_rate_by_name` returns a writable `ChemRate *` into the mechanism (`src/mechanism.c:205-211`); writing an out-of-range `type` there would make `rate_family` read past its table (`:220`). API misuse, low likelihood; read only.
- The overflow guard in `chem_state_create` (`src/state.c:17`) does not bound the `6·n·cells` stage allocation (`:32`) when `nnz + reactions + parameters + 8 < 5n`; unreachable in practice. Read only.
- The 2.03e-13 gate difference (not bit-identical) was not explained: rate formulas, accumulation orders, LU order, and the Markowitz rule all match MICM's text; species order should match (`ts1_mechanism.hpp` declaration order is the phase order). Not investigated further.
