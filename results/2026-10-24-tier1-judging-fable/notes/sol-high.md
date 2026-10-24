# sol-high: Judge 2's reading notes (facts with citations; provisional ratings at the end)

Final commit 1d944f4. Worktree ~/NSF-NCAR/tier1-work/verify/sol-high/clang/sol-high.

## Shape
- One generic path for every mechanism: `ChemMechanismData` descriptors (generated C data,
  `tools/configure.py:215-245`) -> `chem_mechanism_create` builds the dense pattern, the fill,
  the row/diagonal maps, the LU operation list and the derivative jobs (`src/mechanism.c:65-169`).
  No mechanism-specific kernels: forcing and Jacobian are run-time loops over reactions and
  derivative jobs with a per-reaction scratch copy of the rates (`src/process.c:12-46,52-100`).
- TS1 path: `src/bench.c:187` `chem_rosenbrock_parameters(3)` -> `chem_solve` (`src/rosenbrock.c:133-138`)
  -> `integrate` (`src/rosenbrock.c:28-131`, 104 lines) -> `chem_forcing`, `chem_negative_jacobian`,
  `s->linear.factor` / `s->linear.solve` (function pointers, `src/internal.h:29`), `chem_normalized_error`.
- Layout: species-major within a group, `y[i*nc + c]` (`src/rosenbrock.c:21`, `src/process.c:21`);
  default group 32 (`src/bench.c:11-13`), shrunk for small inputs (`src/bench.c:166-171`).
- Benchmark: warm-up solve per group, then the timed loop over STEPS with `omp parallel for schedule(static)`
  (`src/bench.c:188-206`); rates computed before the warm-up (`src/bench.c:191`) and held through the run.

## Readability
- Dense style: several statements per line (`src/rosenbrock.c:110`, `src/mechanism.c:174-176`, `src/rates.c:109`);
  2,292... (common count: 3,034 library+benchmark lines, median function 6 lines, 3.3 % lines with a comment,
  the highest share of the nine).
- Names follow MICM's (`chem_alpha_minus_jacobian`, `chem_normalized_error`, `new_function`, `reject_last`/`reject_more`).
- Rate parameters are positional: `ChemRate.p[24]` with a comment key (`include/chem.h:20-25`); formulas read
  `p[0] * exp(p[2] / t) * pow(t / p[3], p[1])` (`src/rates.c:27-30`), so a reader needs the key for every type.
- The integrator matches MICM's structure line for line (stage loop with `new_function`, `c/h` coefficients,
  `m` and `e` sums, the accept/reject rules with `reject_last`/`reject_more`: `src/rosenbrock.c:55-127`).
- One indentation slip at `src/rosenbrock.c:123`.
- Documents: README (91 lines), DESIGN (155) with an "Extension points" section, DECISIONS (114) detailed.

## Maintainability
- Tests: Python through ctypes (`tests/chemffi.py`), 21 modules each exporting a `CASES` dict keyed by
  `Suite.Test`; the runner walks `tests/cases.txt` (680 ids) and, per id, calls the Python check, or the C
  check by name in `bin/test_rates --list`, or prints `FAIL case port pending` (`tests/run_cases.py:41-64`).
  Several ids may map to one function through the `case` decorator (`tests/test_solvers.py:10-14,24`),
  e.g. Standard and Vector variants share one check. 5 N/A with reasons (`tests/not_applicable.json`).
- Supplemental checks after the 680: reader verification, config-vs-hand-built equivalence for Robertson
  and a parameterized mechanism, a TS1 stiff-cell regression and the 100-cell TS1 fixture at group sizes
  1, 7, 32, 128 against MICM's reference to 1e-10 (`tests/supplemental.py:10-36`).
- The generator runs on four other configurations in the build (`Makefile`: build/no_user.c, user.c,
  robertson.c, parameterized.c).
- Error handling: status codes; `memory_error` cleanup in `chem_mechanism_create` (`src/mechanism.c:165-168`);
  input validation of descriptors (`src/mechanism.c:6-63`); CSV errors reported (`src/bench.c:96-112`);
  statuses checked per group (`src/bench.c:208-211`). NaN: the error norm floors with `fmax` but the
  integrator inspects every candidate and error value (`src/rosenbrock.c:98-101`), so a NaN state returns
  CHEM_NONFINITE.
- Logic in two places: diagonal Markowitz ordering in Python (`tools/configure.py:167-205`, build time)
  and in C (`src/linear.c:150-180`, run time, offered through the API); fill only in C.
- Build: one compile line for the library and benchmark (`Makefile`, bin/ts1_bench_cpu), `build/compiler-flags`
  stamp; `-Wall -Wextra -Wpedantic`. Verified without warnings? Not recorded for Sol (results/2026-10-15-sol-verification).
- Hidden coupling on the rate enum: `valid_data` rejects `rate->type > CHEM_CALLBACK` (`src/mechanism.c:30`)
  and `rate_family` indexes a 12-entry table by type (`src/mechanism.c:213-221`), so a type appended to
  `ChemRateType` (`include/chem.h:6-10`) is rejected or reads past the table unless both change.

## Exercises (Judge 2's own reading)
- E1: public struct `ChemRosenbrock` (`include/chem.h:104-114`), generic stage loop to 6 stages
  (`src/rosenbrock.c:30,61`), sets as `case`s of `chem_rosenbrock_parameters` (`src/parameters.c:6-181`);
  the benchmark picks 3 at `src/bench.c:187`. DESIGN says how. Ease 5 (fill a struct or add a case).
- E2: hook `ChemStepController(error, order, safety, factor_min, factor_max, context)` (`include/chem.h:93-94`,
  called at `src/rosenbrock.c:102-104`) returns the factor only; it receives neither h nor the acceptance,
  and the context lives in the parameter struct shared by all threads (`src/bench.c:187,202`). A PI controller
  needs a new signature and per-state storage. Ease 3.
- E3: `ChemLinearSolver {factor, solve}` per state, `chem_set_linear_solver` (`include/chem.h:95-98,209`,
  `src/linear.c:181-186`), called at `src/rosenbrock.c:58,82` and `src/backward_euler.c:44,48`; receives
  (mechanism, cells, packed matrix) with `chem_jacobian_index` for the storage. Doolittle shows the pattern.
  Ease 5.
- E4: C: enum, `chem_rate_evaluate`, `chem_rate_default`, `valid_data` range, `rate_family` table
  (5 sites, 3 files); Python: `RATE_PARAMETERS`, the kind map, the unit scaling (3 sites); tests 2.
  `CHEM_CALLBACK` with function/context exists for C-side rates. Ease 3.
- E5: `tools/configure.py CONFIG OUT.c --symbol NAME` + a Makefile rule; the benchmark hard-codes
  `ts1_data` (`src/bench.c:14`), so a new mechanism runs through it only by generating with `--symbol ts1_data`
  or editing one line; generic kernels, no generated code; reader remaps indices consistently when it reorders
  (`tools/configure.py:200-204`); the tests run the generator on four configurations. Ease 4.
