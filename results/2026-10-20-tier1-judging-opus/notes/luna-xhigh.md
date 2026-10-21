# luna-xhigh: notes for the judge

Facts from the final commit `ada55d5`, read-only. Paths are in the entry unless marked MICM.

## Tracing the TS1 path
- The run-time path is two files, `src/ts1_bench.c` and `src/chem.c`, plus the generated tables in `build/generated/ts1.c`. `main` (`src/ts1_bench.c:256-394`) calls `chem_integrate_group` for the warm-up (`src/chem.c:1545-1557`) and `chem_integrate_group_steps` for the timed run (`:1516-1529`); both reach `chem_integrate_group_steps_with_parameters` (`:1061-1267`), which calls `chem_rosenbrock_parameters_valid` (`:323-353`), `chem_symbolic_build` (`:444-620`) and through it `chem_markowitz_order` (`:376-442`), `chem_allocate` (`:1053-1059`), `chem_calculate_rates` (`:622-638`) with `chem_rate_value` and `chem_falloff` (`:52-151`), then per step `chem_add_forcing`, `chem_build_jacobian`, `chem_factor`, `chem_solve`, `chem_normalized_error`, and `chem_step_size_factor` (`:640-1051`), and `chem_symbolic_destroy` at the end (`:355-374`): 19 library functions, 14 of them `static`.
- One `chem_` prefix throughout; the integrator has three thin wrappers (`src/chem.c:1516-1557`). `include/chem.h` (502 lines) holds the mechanism descriptor, solver parameters, state, objects, dense and sparse LU, backward Euler, and the integrators' prototypes together.

## Large functions and logic kept in two places
- Large: `chem_integrate_group_with_constraints`, 233 lines (`src/chem.c:1282-1514`); `chem_integrate_group_steps_with_parameters`, 207 (`:1061-1267`); `chem_constraint_set_initialize`, 199 (`src/constraints.c:319-517`); `chem_symbolic_build`, 177 (`src/chem.c:444-620`); `chem_integrate_group_backward_euler_with_parameters`, 169 (`src/backward_euler.c:83-251`); the benchmark's `main`, 139 (`src/ts1_bench.c:256-394`); and `tools/test_cases.py`'s `main`, 274 lines, one `if`/`elif` chain over suites (`tools/test_cases.py:55-328`). `src/chem.c` is 1,744 lines: rates, five parameter sets, symbolic analysis, kernels, both Rosenbrock integrators, and finite-difference Jacobian checks.
- Forcing twice: species-major sparse `chem_add_forcing` (`src/chem.c:640-677`) and cell-major `chem_accumulate_forcing` (`:706-742`). Jacobian twice: `chem_build_jacobian` (`:841-880`) and dense cell-major `chem_accumulate_negative_jacobian` (`:755-792`). Error norm twice: `chem_normalized_error` (`:1023-1046`) and `chem_calculate_normalized_error` (`:809-839`). Rate constants twice: `chem_calculate_rates` and backward Euler's `calculate_rates` through `chem_evaluate_rate` (`src/backward_euler.c:67-81`).
- The Rosenbrock step loop is written twice: sparse in `chem_integrate_group_steps_with_parameters` and dense per cell in `chem_integrate_group_with_constraints` (`src/chem.c:1282-1514`).
- Four LU implementations: sparse Mozart and Doolittle (`src/chem.c:882-979`), dense Mozart and Doolittle (`src/linalg.c:12-44`), coordinate sparse (`src/linalg.c:119-203`), and a partial-pivoting dense LU for constraint initialization (`src/constraints.c:227-275`).
- The species order is computed only in C (`src/chem.c:376-442`); the generator emits the input order (`tools/generate.py:172-173`). No Python and C duplication there.
- `config/tutorials/rate_constants_*` is identical to `config/rate_constants_*` (checked with `diff -r`); only the latter is referenced (`Makefile:15-18`; `grep` finds no use of `config/tutorials`).

## Dead or unused code
- `chem_species.flags` (`include/chem.h:70`) is never read; `CHEM_SPECIES_THIRD_BODY` (`include/chem.h:26`) appears only in a test initializer (`tests/test_processes.c:72`).
- `chem_symbolic.reaction_count` is written (`src/chem.c:452`) and never read.
- `if kind == "BranchedRateConstant": pass` (`tools/generate.py:267-268`).
- `chem_alpha_minus_jacobian` (`src/chem.c:794-807`) is called only by a test (`tests/test_rosenbrock_utilities.c:32`).

## Error handling
- Public APIs return 0, −1 for invalid input, and −2 for allocation failure (`include/chem.h:274-276,378-382`). The integrator returns −1 for every failure, allocation or convergence, and prints the reason to stderr from inside the library (`src/chem.c:1091,1127,1132,1157-1158,1207`); the step-limit message hard-codes "1000 internal steps" although the limit is a parameter (`:1126-1127`).
- Sizes are checked for overflow before allocation (`src/chem.c:449-464,1054-1055`); allocation failures unwind through `goto cleanup` (`:1086-1098,1253-1266`) and `chem_symbolic_destroy` (`:480-487`).
- A zero or non-finite pivot stops the run, reporting the cell and pivot (`src/chem.c:893-898,1155-1160`); a non-finite error norm stops it too (`:1206-1210`).
- Benchmark: bad arguments exit 2 (`src/ts1_bench.c:215-224,259-278`); thread and cell counts are parsed as doubles and truncated (`:264,271-272`); unknown CSV columns are an error (`:313-317`); a failed group exits 1 (`:347-353,370-376`). The CSV reader does not handle quoted fields (`:71-87`).

## Documentation versus code
- DESIGN.md says the cell loops use `#pragma omp simd` "where applicable" (`DESIGN.md:112-114`); the flux and derivative loops carry none (`src/chem.c:652-660,855-862`). The transcript shows the agent removing two pragmas, on the rate-constant and derivative loops, after clang's "loop not vectorized" warnings (09:40).
- README says `tools/translate_ts1.py` "keeps the TS1 source-to-configuration translation reproducible" (`README.md:68-69`); Make never runs it, and `config/ts1.json` is checked in (`Makefile:47-49`; `DESIGN.md:9-10`).
- DESIGN.md's recipe for a new rate type (`DESIGN.md:143-145`) does not say that rate parameters are positional in one 24-slot array shared by all kinds (`include/chem.h:26,193`; `tools/generate.py:26-34`): Troe uses slots 5–12, JPL19 adds 18–23, branched 13–16, surface reuses 0.
- Checked and matching: one group per thread and static scheduling, per-call symbolic build and allocation, no Python at run time, the clip rule (`DESIGN.md:5-8,105-130`; `DECISIONS.md:25-33,60-63`).

## Build hygiene
- Clang on macOS gets Homebrew's libomp include and library paths hard-coded (`Makefile:22-25`); GCC gets `-fopenmp` (`:31-32`).
- Each of the 23 test programs has its own hand-written compile and link rules (`Makefile:71-244`), with no pattern rule; every object depends on all six headers (`Makefile:12,53`).
- `make test-cpu` runs the benchmark twice and writes `ts1_final.csv` and `build/ts1_ours_*.csv` into the tree (`tools/test_cases.py:23-32`), both ignored (`.gitignore:3,5`). No warnings with either compiler (verification).
- From 12:04 to 12:16 the agent's GCC runs used `CC=gcc`, Apple clang on this machine; the third turn reran with `gcc-15` (transcript; commit `ada55d5`).

## Apparent defects
1. **Three-stage error weight.** `e[1] = -2.907955871680469821793236208017` (`src/chem.c:197`) is not MICM's `-0.29079558716805469821718236208017e+01` (MICM `include/micm/solver/rosenbrock_solver_parameters.hpp:227`): as doubles they differ by 7.7e-14, 174 ulp. Certain: compared the literals and their double values in Python. The same value was hard-coded in `2811c0c`. It enters only the error estimate (`src/chem.c:1198`); the gate passes at 1.77e-13.
2. **Unlabeled user-defined reactions.** `emit` sets the parameter index and the scaling factor only for a labeled reaction (`tools/generate.py:271-275`), so an unlabeled `USER_DEFINED` or `PHOTOLYSIS` reaction gets index 0 and scale 0, a zero rate, without an error. Fairly sure, from reading; not exercised: all 142 TS1 user-defined reactions and the tutorial's have labels (checked in `config/ts1.json` and `config/rate_constants_user_defined/reactions.json`).
3. **Test mappings.** MICM's four `DiagonalMarkowitzReorder` cases reorder a random 50 × 50 pattern and assert less fill (MICM `test/unit/solver/test_linear_solver_policy.hpp:411-470`); here they run the diagonal-matrix solve (`tests/test_lu.c:263-264,270,274,280`). The 52 LU and linear-solver cases test the dense routines of `src/linalg.c`, not the sparse `chem_factor` used by TS1 (`tests/test_lu.c:42-53,143-152`). Certain, from reading.
4. **Departures, not coding errors:** the clip keeps negatives below −atol (`src/chem.c:1239-1246`; documented, `DECISIONS.md:60-63`), and the species order differs from MICM's because every species gets a diagonal before ordering (`src/chem.c:489-490`): a replay of both procedures on `config/ts1.json` gives 2,359 nonzeros after fill here and 2,381 for MICM's order.
