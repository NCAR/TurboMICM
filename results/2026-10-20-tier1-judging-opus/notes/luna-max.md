# luna-max: notes for the judge

Facts only, at `da2d722`; paths relative to the entry. Not for the wiki.

## Tracing the TS1 path
- Three C files and the generated table: `src/benchmark.c`, `src/solver.c`, `src/kinetics.c`, `build/ts1_model.c`. About 30 functions: 13 in `src/benchmark.c` (CSV parsing, timing, output); 7 for setup (`chem_rosenbrock_default`, `chem_solver_create`, `valid_mechanism`, `chem_make_pattern`, `chem_effective_rate_constant`, `chem_rate_constant`, `falloff`); 10 in the timed path (`chem_solver_step`, `solve_group`, `work_vector`, `evaluate_rhs`, `default_factor_callback`, `factor`, `default_solve_callback`, `solve_lu`, `chem_normalized_error`, `default_step_factor`).
- The whole step lives in one 154-line function, `solve_group` (`src/solver.c:644-797`): attempt loop, cell loop, stage loop, error norm and step control.
- One header, `include/kinetics.h` (308 lines), declares rates, mechanism, state, species, solver and backward Euler; constraints, matrices and system objects have their own headers.

## Size and naming
- Largest functions: `solve_group` 154 lines; `chem_constraint_set_project` 155 (`src/constraint.c:470-624`); `main` 114 (`src/benchmark.c:259-372`); in the tests, `test_analytical_integration_case` 239 (`tests/test_cpu.c:1419-1657`) and `run_case` 176 (`:4153-4328`). `tests/test_cpu.c` is one 4,357-line file with 76 functions.
- `chem_` prefix throughout. Terse fields in `struct chem_solver`: `n`, `nr`, `ncustom` (`src/solver.c:21-23`). The benchmark's step loop reuses the variable `cell` as the step counter (`src/benchmark.c:336-341`). `--threads` must come before `--group-size` (`src/benchmark.c:270-283`).

## Logic kept in two places
- Forcing: `evaluate_rhs` (`src/solver.c:255-268`) and `chem_mechanism_forcing` (`src/kinetics.c:181-204`). Jacobian: sparse and negative in `evaluate_rhs` (`src/solver.c:271-286`), dense and positive in `chem_mechanism_jacobian_dense` (`src/kinetics.c:206-239`), which backward Euler uses (`src/backward_euler.c:75-76`).
- Mechanism validation: `valid_mechanism` (`src/solver.c:606-642`) and `kinetic_tables_valid` (`src/kinetics.c:156-179`).
- Elimination three times: sparse `factor` (`src/solver.c:326-356`), dense `chem_lu_factor_dense` (`src/matrix.c:558-577`), and the test's `test_factor_callback`, a copy of `factor` without `fma` (`tests/test_cpu.c:2171-2200`).
- Fill computed in Python to choose the order (`tools/generate_c.py:186-189`) and again in C for the pattern (`src/solver.c:118-127`).
- Air density computed twice per rate constant (`src/kinetics.c:28,121`).

## Dead or unreachable code
- `allocation_error` is never read or written (`src/solver.c:35`); `ncustom` is written (`:817`) and never read.
- In `run_case`, the branches at `tests/test_cpu.c:4201-4202` (`solver_lambda_rate_constant`) and `:4261-4262` (`integrated_reaction_rates`) cannot be reached: earlier branches take those suites (`:4186-4188,4250-4253`). The N/A block at `:4181-4185` is repeated at `:4318-4321`. For `grouped_views`, the result of `test_matrix_views(1)` is overwritten (`:4273-4275`).
- Used only by tests: `chem_rate_store_layout` (`src/kinetics.c:127-154`), the `chem_sparse_matrix` API (`src/matrix.c:241-522`), `include/version.h`.

## Error handling
- `chem_solver_step` returns -1 for every integration failure: bad parameters, too many steps, step too small, non-finite error, zero or non-finite pivot (`src/solver.c:666-678,698-707,767-768,997-1002`). The `chem_solver_state` codes (`include/kinetics.h:134-139`) come back only for a failed constraint initialization (`src/solver.c:969-970`).
- A failed group leaves other groups advanced and its own cells partly advanced; there is no rollback (`src/solver.c:775-776,997-1002`). Backward Euler, by contrast, restores the cell on failure (`src/backward_euler.c:99,108`).
- The main allocations have no overflow guards (`src/solver.c:837-839`); `chem_make_dense_pattern` and `chem_solver_set_constraints` do check (`:209-214,899-908`). Allocation failures return NULL or -1 throughout.
- With constraints, `evaluate_rhs` mallocs and frees a gradient on every call (`src/solver.c:288-295,321`).
- The thread count is read once, at solver creation (`src/solver.c:823`); a later `omp_set_num_threads` changes neither the group count nor the team size.
- The CSV reader rejects unknown columns, short rows and bad numbers (`src/benchmark.c:160-203`).

## Documentation against code
- DESIGN.md:72 lists OpenMP scheduling as overridable by "compiler/runtime OpenMP setting"; the code fixes `schedule(static)` (`src/solver.c:980`).
- DESIGN.md:47-48: "factors `alpha * I - J` in place". `factor` copies the negative Jacobian into a separate buffer and factors the copy (`src/solver.c:330-332`).
- The comment "one sparse LU workspace per solver group" (`src/solver.c:34`) describes a per-thread allocation indexed by thread (`:839,655`).
- RESULTS.md is from `2bef5f9`; `da2d722` changed the LU to `fma` (`src/solver.c:345-378`, `src/matrix.c:569-593`) without updating it. The verified clang result of `da2d722`, 2.78e-13 (BCARY, cell 4179), equals the value RESULTS.md reports from before the change.
- `bin/test_cpu` alone prints PASS with "performed by the Python test driver" for the eight programs and FlowTube (`tests/test_cpu.c:4162-4171,4263-4266`); only `tests/run_cases.py` fills those rows in (`tests/run_cases.py:132-144`).

## Build hygiene
- `CC = clang` is set unconditionally, so `CC` from the environment is ignored; `make CC=gcc-15` works (`Makefile:1`).
- Object rules list `include/kinetics.h` but not `include/constraint.h`, which it includes (`include/kinetics.h:6`; `Makefile:61-89`); no generated dependency files.
- `SDKROOT_FLAG` is defined after the rules that use it, which works because recipes expand late (`Makefile:97,99`).
- The benchmark links object files, not `lib/libchem.a` (`Makefile:33,96-97`).
- `make micm-ref` compiles MICM's reference driver against `../micm/include` (`Makefile:149-150`); `cpu` and `test-cpu` do not use it.
- gcc-15: three `-Wint-in-bool-context` warnings (`src/solver.c:837-838`; `src/benchmark.c:307`); the tests add 38 `-Wmissing-field-initializers` (verification README).
- `tests/run_cases.py` loads `tools/generate_c.py` by path (`tests/run_cases.py:15-20`); `bin/test_cpu` opens `../reference/micm_cases.txt` relative to the working directory (`tests/test_cpu.c:4332`).

## Possible defects
- Stage reuse (high confidence, by reading `src/solver.c:703-728` against MICM `rosenbrock.inl:165-168`; not run): a stage's forcing is copied into the next stage only inside `if (new_function_evaluation[stage])`, or from stage 0 to stage 1. With two consecutive stages flagged 0, the second starts from whatever its work vector holds: the cell's previous attempt, or uninitialized `malloc` memory on first use (`src/solver.c:838`). No built-in set has such a pair (`src/solver.c:414-416,464,487-490,528-531,561`), but DESIGN.md:88-90 invites user-defined sets.
- Test coverage of the benchmark's LU: the 52 LU and linear-solver ids test only the dense routines (`tests/test_cpu.c:2558-2609,4313-4316`); the sparse `factor` and `solve_lu` are exercised only through integration cases and the 100-cell TS1 comparison.
- The six tutorial-program ports check finiteness and signs, not MICM's printed values (`tests/test_programs.c:39-184`).
- Scale: `chem_make_pattern` uses a dense n × n byte matrix and an n × n `size_t` position map (`src/solver.c:102,136`), 8n² bytes for the map (349 KB at n = 209).
