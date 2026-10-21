# luna-high: notes for the judge

Facts bearing on readability and maintainability, cited at `d3c2480`. Not for the wiki.

## Tracing the TS1 path

- Three files from driver to kernels: `main` (`src/main.c:199-391`) → `chem_solver_workspace_init_group` (`src/solver.c:68-132`), which calls `chem_jacobian_plan_init` (`src/mechanism.c:433-486`) and `make_pattern` (`src/solver.c:254-303`) → `chem_rosenbrock3_solve_group_with_workspace` (`src/solver.c:767-945`) → `chem_rate_constants` (`src/mechanism.c:243-295`, with `chem_reaction_rate` and `chem_surface_rate`), `chem_forcing_with_constants` → `chem_forcing_impl` (`src/mechanism.c:297-354`), `chem_negative_jacobian_with_plan_pattern` (`src/mechanism.c:516-580`), `factor` and `solve_factored` (`src/solver.c:305-344`). About a dozen library functions; the mechanism arrives through one `extern` (`src/main.c:14`).
- The benchmark's integrator is a separate, three-stage-only copy of the Rosenbrock loop (`src/solver.c:805-938`) beside the general one (`src/solver.c:558-680`); reading the general path does not show what TS1 runs.
- Comments: 10 comment lines in the 4,667 lines of C library and benchmark (`src/mechanism.c:284,539`; `src/solver.c:287,892`; `include/chemistry.h:140,193,197,257,475,606`). The header's 146 public function declarations are otherwise uncommented.
- Positional parameters: each reaction's `p[24]` has per-kind meanings (`include/chemistry.h:51`; `tools/generate_c.py:147-198`); `p[7]` is the custom-parameter index for user-defined and surface rates (`src/mechanism.c:153,234`) and Troe's N (`src/mechanism.c:162`). The group loop's coefficients are named `aa`, `cc`, `mm`, `ee` (`src/solver.c:796-803`).

## Size

- 199 functions, median 10 lines; 10 over 60 lines and 5 over 100 (judge's metrics): `rosenbrock_solve_impl` 202 (`src/solver.c:489-690`), `main` 193 (`src/main.c:199-391`), the group solver 179 (`src/solver.c:767-945`), `project_constraints` 135 (`src/solver.c:353-487`), backward Euler 124 (`src/solver.c:972-1095`). The test program's `main` is 645 lines (`tests/test_core.c:4801-5445`).
- `src/main.c` repeats the same seven-pointer `free` list at eight exits (`src/main.c:255,284-285,297,304,325,355,370,388`).

## Logic kept in two or more places

- Rosenbrock loop and step controller: twice (`src/solver.c:558-680,805-938`); the three-stage coefficients twice (`src/solver.c:237-249,796-804`).
- Normalized error: three times (`src/solver.c:642-655,890-907`; `src/matrix.c:68-102`), the last called only by the tests (`tests/test_core.c:3227`).
- Markowitz ordering in Python (`tools/generate_c.py:40-86`) and in C (`src/linear.c:77-138`, called only by the tests, `tests/test_core.c:4282`). They differ: the C routine rejects a pattern with a missing diagonal (`src/linear.c:93-97`); the Python one keeps such species, with a wrapped cost (`tools/generate_c.py:62-67`).
- Symbolic fill in `make_pattern` (`src/solver.c:288-292`) and again in a test (`tests/test_core.c:4298-4309`).
- Two Jacobian formulas: `chem_jacobian_impl` divides the rate by the concentration (`src/mechanism.c:386-387`), behind the public `chem_jacobian` (`src/mechanism.c:417-422`); the solvers use the plan's product of the other reactants (`src/mechanism.c:552-561`).
- Arrhenius and Taylor-series formulas twice: with pressure rebuilt from density in `chem_reaction_rate` (`src/mechanism.c:133-138,193-203`) and with the given pressure in `chem_rate_constants` (`src/mechanism.c:261-274`).
- The gas constant as the literal 8.31446261815324 in seven places (`src/main.c:278`; `src/constraints.c:71`; `src/solver.c:957`; `src/mechanism.c:135,199,209,236`). Rate-kind tables with both upper- and lower-case keys (`tools/generate_c.py:8-37`). The Makefile's source list four times (`Makefile:39,41,44,46`).

## Dead or unreachable code

- `tools/generate_c.py:202-205` converts dict-form reactants and products, but dict terms would already fail at `term.get` (`tools/generate_c.py:125-126`); only lists reach it, and it skips them. By reading.
- `src/main.c:312`: after `src/main.c:281-288` has returned on any bad cell, this loop has nothing to find.
- 5 of the 667 ids in `implemented` match no case id, for example `phase_species:PhaseSpecies.SetDiffusionCoefficient` (`tools/test_cases.py:7-599`).

## Error handling

- Status codes throughout, no `abort`; allocation failures in `chem_solver_workspace_init_group` return `CHEM_SOLVER_STATUS_BAD_ARGUMENT` (`src/solver.c:85-86,115-127`), and the status enum has no out-of-memory value (`include/chemistry.h:650-658`).
- Validation inside the hot loops: every forcing evaluation re-checks each term's species index (`src/mechanism.c:307-308,317,326,332`), every Jacobian evaluation each plan entry (`src/mechanism.c:547-551`).
- `status` holds per-cell input errors and per-group solve results in one array (`src/main.c:252,276,312-318,336`). In the warm-up's parallel loop each thread reads `status[c]` for its cells while threads 1 to 9 write `status[1..9]`, inside group 0's range (`src/main.c:312,315`): a data race in C11 terms, harmless while every value is 0; after a failed group it can skip group 0's warm-up. By reading.
- `src/main.c:242` divides by `nr`, the custom-parameter count, while `src/main.c:317` guards `nr == 0`: a mechanism without custom parameters divides by zero. By reading.
- Test exit codes collide: 52 at `tests/test_core.c:4802` and `:5403`; 70 at `:5442`, `70+status` at `:4808`, `70 + backward_euler_status` at `:5410`.
- The CSV reader's line buffer is 262,144 bytes (`src/main.c:69-76`); a longer line becomes a malformed row, which is reported.

## Docs against code

- DESIGN.md: the factorization "can be replaced behind a linear-solver interface" (`DESIGN.md:112-113`); `factor` and `solve_factored` are static and called directly (`src/solver.c:305-344,591,624,854,877`).
- DESIGN.md: controllers "separate from the method coefficient tables, so either can be changed independently" (`DESIGN.md:110-112`); the controller is inline in both loops (`src/solver.c:654-677,907-935`) and the benchmark's loop has its own coefficients (`src/solver.c:796-804`), so a new parameter set or controller reaches TS1 only by editing that loop.
- The prompt asks for TS1 "read by your configuration reader from your translation"; TS1's JSON is the entry's own format (`tools/import_ts1.py:113`), read by the generator, not by `tools/read_config.py` (`Makefile:36-37`).
- README.md says "Implemented cases print `PASS`" (`README.md:20-23`); the lines come from a list once one program has passed (`Makefile:48-51`; `tools/test_cases.py:640-778`), and 301 of the 662 carry the default reason "C rate formula checked" (`tools/test_cases.py:652-654`), 184 of them sparse-matrix cases.

## Build hygiene

- No objects or archive: each binary compiles all its sources in one command (`Makefile:39-41,44-46`). `make test-cpu` depends on `cpu` (`Makefile:48`), so the tests need Python and build the benchmark, which they do not run.
- `mechanisms/ts1.json` (18,861 lines) is tracked with no rule to rebuild it from the header (`Makefile:36-37`); `bin/` and `generated/` are ignored (`.gitignore:4-5`). The verification found no warnings with either compiler.

## Possible defects

1. CAMP index mismatch (known; `tools/read_config.py:83-86`, `tools/generate_c.py:100-102,126`). Reproduced again for these notes with copies of the two scripts in scratch: species A, B, C, D with A + B → C and D → A gave the species table D, A, B, C and the reactions D + A → B and C → D. Sure.
2. A mechanism without custom parameters gets empty initializers for `custom_labels` and `default_custom_parameters` (`tools/generate_c.py:215-219`): with `-std=c11 -Wpedantic`, clang warns (C23 empty initializer, zero-size array) and gcc-15 stops with "zero or negative size array". Checked with `-fsyntax-only` on the generated file in scratch. Sure. MICM's three tutorial configurations all have labeled reactions, so none hits it; such a mechanism would also meet the division at `src/main.c:242`.
3. Constrained single-cell path: the accepted-step projection uses `jac` as scratch (`src/solver.c:636-638`), and `project_constraints` zeroes it and refills it with the constraint Jacobian only (`src/solver.c:388-393`). If the attempt is then rejected, the retry copies that `jac` into `a` (`src/solver.c:584-588`); the kinetic Jacobian is computed once per step (`src/solver.c:571-573`). Constrained solves with a rejected step only; not TS1. By reading, not run; fairly sure.
4. The rejected-step counter counts every rejection (`src/solver.c:675,933`), changed at 17:26 (`fa93b85`) from MICM's rule (MICM `include/micm/solver/rosenbrock.inl:276-279`). Counter only.
5. The `status` race in the warm-up loop (above). Sure of the race; it matters only after a failure.
