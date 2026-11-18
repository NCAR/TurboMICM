# The test sample in the three Terra entries

Judge 1's twelve cases (wiki Tier-1-Judging, "The test sample"), traced by Judge 1 itself in
terra-max, terra-xhigh and terra-high at their final commits, read-only, in
`~/NSF-NCAR/tier1-work/verify/<run>/clang/<run>/`. MICM's assertions for each case are those
summarized in `results/2026-10-20-tier1-judging-opus/tests/T1.md` and `T2.md`, and the grades
follow the rules stated there (F, P, S, M, N/A; "How the grades were applied").

## terra-max (9071941)

**Runner.** `make test-cpu` builds 33 programs and runs `tests/run_cases.py`
(`Makefile:194-195`). It runs 27 C programs and one Python ordering test; each prints
`<id> PASS|FAIL <reason>` from a check run for that id, and the runner raises on a line whose
id is unknown, repeated, or malformed (`tests/run_cases.py:19-28`). Eight tutorial programs
are compared with MICM's printed output, timings and thread lines masked (`:30-52`); four core
programs fail the run on a nonzero exit but report no id (`:54-58`); the four SFINAE ids are
N/A (`:60-66`); an id no program reports is FAIL, "case has not yet been ported" (`:68-71`).
Some checks serve several ids: the four Markowitz ids share one Python check
(`tests/reference_ordering.py:50-55`), and `RegressionRosenbrock.Jacobian` and `VectorJacobian`
one function (`tests/reference_chapman.c:273-274`).

| # | Case | Grade | Check | Id printed at | Explanation |
|---:|---|---|---|---|---|
| 1 | Robertson | F | `tests/reference_nonlinear_integration.c:41-80` (advance :25-39) | `:208-212` | MICM's mechanism, conditions, air density 1e6, rtol 1e-10, atol 1e-12, twelve outputs, criterion (:75) and reference table (:52-64). Five parameter sets with Doolittle (:208-209) and the three-stage set with the default in-place Mozart (:211), as MICM pairs them; one cell has no width variants, and the integrator factors only in place. |
| 2 | RejectedStepAlpha | P | `tests/reference_rosenbrock.c:69-88` | `:101` | MICM's mechanism, h_start 1e3, 1e12 everywhere, T 298.15, solve(1e3); steps − accepted ≥ 2 and equal counts (:86). Both runs factor in place: the default Mozart plan against Doolittle (:82-84). |
| 3 | LU, random, vector | P | `tests/reference_lu.c:79-117`, `check` :119-135 | `:147-150` | MICM's 10 × 10 pattern with p = 0.5 and the diagonal, lognormal(−2, 2) values copied to every cell (:39,56-61), groups of 5 and 1-4 cells (:125-131), L·U against A at 1e-12 (:91-99), with an LCG and Box-Muller. The routine checked is `chem_mozart_linear_solver`, the library's searching Mozart (`src/linear.c:306-345`), not the precomputed plan `chem_factor` that the solver runs (`src/linear.c:167-199,303-305`); no test calls `chem_factor` directly. |
| 4 | Markowitz reorder | P | `tests/reference_ordering.py:32-47` | `:51-55` (four ids) | A random 50 × 50 pattern with p = 0.5 (:36-39), reordered by the generator's minimum-fill `reorder_mechanism`, the ordering TS1 uses (`tools/generate.py:12-57`); the reordered pattern must equal the original under the permutation (:44-47), and the fill must fall strictly (:47). The fill is counted by the test's own elimination (:23-29), not by the library's symbolic plan. |
| 5 | Chapman, factored αI − J | F | `tests/reference_chapman.c:149-160` (pattern :125-140) | `:263` | −1 at the 23 positions, then the library's in-place Mozart with alpha 2; the 23 values at 0.01, diagonals as reciprocals (:156-157). |
| 6 | Regression, vector Jacobian | P | `tests/reference_chapman.c:89-112` (oracle :59-69) | `:274` | One 3-cell state, all seven rates user-set, lognormal(−2, 4) from a seed-739 minstd (:91-96); all 81 entries of each cell against a hand-written Jacobian at MICM's tolerance (:103-105). It reads `chem_evaluate`, whose reaction terms come from `add_reaction_terms`, a per-cell dense helper (`src/rates.c:280-323`); the integrator uses `chem_negative_jacobian` (`src/rates.c:234-278`; `src/rosenbrock.c:173`). Same code as the Jacobian id, no width variants. |
| 7 | Troe, all arguments | F | `tests/reference_rates.c:98-114` | `:164-166` | MICM's T, density and eight parameters (:103-104), MICM's expression (:110-111), tolerance 0 (:112), through `chem_rate_calculate`. |
| 8 | Surface, missing diffusion | F | `tests/reference_rate_store.c:130-137` | `:145` | A surface rate with diffusion 0 and MW 0.025 must give CHEM_MISSING_DIFFUSION (`src/rates.c:145`). The rate carries the species' properties; state construction does not validate rates, so evaluation is where the error appears. |
| 9 | State, invalid concentrations | P | `tests/reference_state.c:84-104` | `:160` | A 3-cell state of four species and ten reactions; the bulk setter must return CHEM_INVALID and leave every value 0 (:95-101). Each entry carries three values and only the first name is unknown (:92-94), so only the unknown-name path runs; MICM's one-value entries are absent. |
| 10 | Backward Euler, Oregonator | F | `tests/reference_nonlinear_integration.c:156-203` | `:214-215` | MICM's mechanism, initial values, absolute tolerances (:160-161), relative tolerance 1e-9, MICM's backward-Euler defaults (`src/backward_euler.c:8`), 12 outputs of 30τ in 18,000 substeps that must each succeed and advance (:25-39,192), relative error below 5e-3 against MICM's table (:171-183,195); default Mozart and Doolittle (:214). |
| 11 | DAE, conservation | F | `tests/reference_external.c:92-108` (model :61-83) | `:364` | A_GAS and AEROSOL.A_AQ, a transfer at k = 0.1 and the algebraic constraint A_GAS + A_AQ − 1 in the A_AQ row (:68-76), no reactions, the four-stage DA set, 0.9 and 0.1, 298 K, 101325 Pa; 20 solves of 10 s, each CHEM_OK with the sum within 1e-4 (:101-105). |
| 12 | User-defined tutorial | F | `tests/tutorial_program.c:24-50` | `tests/run_cases.py:30-34` | MICM's tutorial configuration, byte-identical to MICM's files, through the reader and generator (`Makefile:152-153`); MICM's run (:30-46); stdout compared byte for byte with `tests/expected/user.txt`, whose rows equal MICM's documented table. |

Faithful: 7 of 12 (1, 5, 7, 8, 10, 11, 12).

## terra-xhigh (288f274)

**Runner.** `make test-cpu` builds the test programs and runs `tools/test_cpu.py`
(`Makefile:169-171`). It runs 25 C programs, each printing `<id> 1|0` from a check run for that
id; a repeated id raises (`tools/test_cpu.py:10-16`), and so does any reported id that is not in
the case list (`:64-65`). `test_solver`'s exit status is ANDed into
`BackwardEuler.CanCallSolve` (`:17-18`); the two by-hand tutorials and six tutorial programs are
compared with stored output, nanoseconds masked and the OpenMP tutorial's concurrent lines cut
(`:19-40`); `openmp_unit`'s exit status decides `OpenMP.OneSolverManyStates` (`:41-42`); four
SFINAE ids are N/A (`:43-48`); an id no program reports is FAIL (`:60-62`). Some checks serve
several ids: the Markowitz check prints both prefixes' ids (`tests/test_lu_cases.c:154-159`),
and every state id runs at widths 1 and 4 (`tests/test_state_cases.c:103`).

| # | Case | Grade | Check | Id printed at | Explanation |
|---:|---|---|---|---|---|
| 1 | Robertson | F | `tests/test_analytical.c:124-160` | `:183-185` | MICM's mechanism, conditions, air density 1e6, rtol 1e-10, atol 1e-12, twelve outputs, criterion (:156) and reference table (:130-143). Five parameter sets × Mozart and Doolittle × the single-cell integrator and the grouped one through a grouped matrix of width 1-4 (:145-154,183-184): 50 runs. |
| 2 | RejectedStepAlpha | P | `tests/test_rosenbrock_cases.c:60-77` | `:94` | MICM's mechanism, h_start 1e3, 1e12 everywhere, T 298.15, solve(1e3); steps − accepted ≥ 2 and equal counts (:75). The two runs are the single-cell and the grouped integrator with the same Mozart on the same storage (:73-74); the library factors a copy of the Jacobian in every path (`src/rosenbrock.c:292-300`). |
| 3 | LU, random, vector | F | `tests/test_lu_cases.c:21-67` | `:131-141` | MICM's 10 × 10 pattern with p = 0.5 and the diagonal (:31-34), lognormal(−2, 2) values drawn once and copied to the five blocks of a library sparse matrix of width 1-4 (:36-45), each block factored by `chem_linear_factor_in_place`, which runs the integrators' `factor` (`src/rosenbrock.c:287-331,365-395`); L·U against A at 1e-12, with L zero above and U zero below the diagonal (:60-63). |
| 4 | Markowitz reorder | P | `tests/test_lu_cases.c:104-124` | `:154-158` | MICM's 50 × 50 pattern (:109-112), `chem_diagonal_markowitz_order`, a full permutation and sum check (:114-120), and a strict > on the library's symbolic fill (:121-122), at widths 1-4. The C function is not the ordering TS1's tables come from: the generator applies its own Python copy of the rule (`tools/generate.py:17-45,141`), which no case id reaches. |
| 5 | Chapman, factored αI − J | F | `tests/test_chapman.c:24-49` | `:225` | −1 at the 23 positions and 2 on the diagonal, `chem_linear_factor` with Mozart, the 23 values at 0.01 with diagonals as reciprocals (:29-48). |
| 6 | Regression, vector Jacobian | P | `tests/test_chapman.c:121-138` (oracle :94-106) | `:233-238` | Three cells, `chem_jacobian`, which runs the integrators' `chem_jacobian_add` (`src/chem.c:190-221`; `src/rosenbrock.c:583`), against a hand-written Jacobian at the 23 positions and MICM's tolerance (:132-135). Inputs are uniform on [0.01, 3.01) from a seed-739 LCG (:127-128), not lognormal(−2, 4); the width loop fills grouped matrices that the check never reads (:123,129). |
| 7 | Troe, all arguments | F | `tests/test_cases.c` (`falloff`) | `:659` | MICM's T, density and eight parameters, MICM's expression, exact equality. |
| 8 | Surface, missing diffusion | F | `tests/test_cases.c` (`missing_surface`) | `:672` | A surface reaction with diffusion 0 and MW 0.025; `chem_validate_mechanism` must return CHEM_MISSING_DIFFUSION with a message naming the species and the coefficient (`src/chem.c:28-34`). |
| 9 | State, invalid concentrations | F | `tests/test_state_cases.c:65-66` (fixture :9-13) | `:102-103` | The 3-cell foo/bar/baz/quz state with ten rate constants and quux and corge; the bulk setter gets FUU and bar with one value each, MICM's two kinds of invalid input, and must fail, at widths 1 and 4. |
| 10 | Backward Euler, Oregonator | F | `tests/test_stiff.c:78-106` | `:111-118` | MICM's mechanism, initial values, absolute tolerances, rtol 1e-9, MICM's backward-Euler defaults, 12 outputs of 30τ in 18,000 substeps that must each succeed and advance (:102,`advance`), relative error below 5e-3 against MICM's table; both LU algorithms, single-cell and grouped (:110-111). |
| 11 | DAE, conservation | F | `tests/test_external_constraints.c` (`basic`, kind 2) | `:246-247` | A_GAS and AEROSOL.A_AQ, a transfer at 0.1 with the algebraic row for A_AQ, the four-stage DA set, 0.9 and 0.1, 298 K, 101325 Pa; 20 solves of 10 s, each succeeding with the sum within 1e-4. |
| 12 | User-defined tutorial | F | `tests/test_tutorial.c` | `tools/test_cpu.py:20-25` | MICM's tutorial configuration, byte-identical to MICM's files, and a hand-built copy, checked against each other; MICM's run; stdout compared with `tests/tutorial_expected.txt`, whose rows equal MICM's documented table. |

Faithful: 9 of 12 (all but 2, 4, 6).

## terra-high (a8a2ea0)

**Runner.** `make test-cpu` builds 57 targets and runs `tools/test_cpu.py`
(`Makefile:214-216`). Its `CASES` table maps each id to one of 47 programs (or the tutorial
script) with a fixed reason (`tools/test_cpu.py:9-293`); each program runs once, and every id
mapped to it is PASS if it exited 0 (`:317-323`). Five ids are N/A (`:295-301`), and an id not in
the table is FAIL, "not yet ported" (`:313-314`). No id selects a check: all the ids mapped to
one program share its exit status, so, as for luna-high, the grades below describe the check in
the mapped program that corresponds to the case, and M would mark a case with none. The
programs stop at the first failed `assert`. Three tutorial and Robertson fixtures are generated
from the neighbouring MICM tree (`Makefile:73-74,147-150`).

| # | Case | Grade | Check | Id mapped at | Explanation |
|---:|---|---|---|---|---|
| 1 | Robertson | P | `tests/test_nonlinear.c:23-57` (advance :12-22) | `tools/test_cpu.py:172` (`build/test_nonlinear`) | MICM's Robertson configuration, generated from MICM's own file (`Makefile:73-74`), with k = 0.04, 3e7, 1e4 set directly; rtol 1e-10, atol 1e-12, twelve outputs, criterion (:51) and reference table (:25-37). Five sets × dense and sparse × Mozart and Doolittle, 20 runs, all factored in place; the integrator's separate-factor option (`factor_out_of_place`, `include/rosenbrock.h:28`), MICM's default storage, is not run. |
| 2 | RejectedStepAlpha | F | `tests/test_rejected_factorization.c:13-49` | `:50` (`build/test_rejected_factorization`) | MICM's mechanism and rates at 298.15 K, h_start 1000, 1e12 everywhere, solve(1000); each of dense and sparse × Mozart and Doolittle run in place and with separate factors (:28-34); steps − accepted ≥ 2, equal counts, and final states within 1e-14 (:37-43). |
| 3 | LU, random, vector | S | `tests/test_doolittle.c:13-79` | `:232-237` (`build/test_doolittle`, every `lu_decomposition_*` id) | One program for all 32 LU ids: dense 3 × 3, 10 × 10 random, 6 × 6 diagonal and 5 × 5 matrices through the solver's sparse and dense Mozart and Doolittle, in place and copied, from five initial values, L·U against A at 1e-12 (:20-46). Its 10 × 10 matrices are diagonally dominant with uniform entries, a new matrix per block (:58-65), not lognormal(−2, 2) values copied to every block, and there is no vector ordering. |
| 4 | Markowitz reorder | P | `tests/test_reorder.c:22-45` | `:248-250` (`build/test_reorder`) | A 50 × 50 pattern with p = 0.5 and the diagonal (:28-31), stored at widths 0-4, `chem_permutation_init` with a full permutation check (:32-40) and a strict < on fill (:41). The fill is counted by the test's own elimination (:9-17), and the C ordering is not the one TS1's tables come from: the benchmark's generator applies `reference_markowitz_order` in Python (`tools/generate.py:41-63,196-198`; `Makefile:4`). |
| 5 | Chapman, factored αI − J | F | `tests/test_chapman_regression.c:31-57` | `:72-74` (`build/test_chapman_regression`) | −1 at the 23 positions and 2 on the diagonal, `chem_dense_factor`, the integrator's dense Mozart, the 23 values at 0.01 with diagonals as reciprocals (:40-49). |
| 6 | Regression, vector Jacobian | P | `tests/test_nonlinear.c:271-323` | `:162` (`build/test_nonlinear`) | Three cells one at a time, `chem_jacobian`, the benchmark's `chem_accumulate_jacobian` after a clear (`src/chem.c:22-49`), against a hand-written Jacobian over all 81 entries at MICM's tolerance (:303-322). Inputs are deterministic, 0.05-3.2 (:304-305), not lognormal(−2, 4); no widths; the same check backs the Jacobian id. |
| 7 | Troe, all arguments | F | `tests/test_falloff.c:30-38` | `:227` (`build/test_falloff`) | MICM's T, density and eight parameters, MICM's expression, exact equality (:36-38). |
| 8 | Surface, missing diffusion | F | `tests/test_rate_store.c:63-69` | `:78-81` (`build/test_rate_store`) | A store descriptor with no diffusion coefficient; `chem_rate_store_init_diagnostic` must return CHEM_RATE_STORE_MISSING_DIFFUSION for entry 0, with MICM's message, and leave the store empty. |
| 9 | State, invalid concentrations | F | `tests/test_state.c:29-51` | `:181-184` (`build/test_state`) | A 3-cell foo/bar/baz/quz state with ten rates; the bulk setter gets MICM's four one-value entries with FUU and must not return OK (:50-51). |
| 10 | Backward Euler, Oregonator | F | `tests/test_nonlinear.c:135-181` | `:171` (`build/test_nonlinear`) | MICM's mechanism, initial values, absolute tolerances, rtol 1e-9, MICM's backward-Euler defaults (`src/backward_euler.c:10`), 12 outputs of 30τ in 18,000 substeps, each asserted OK (:169-170), relative error below 0.005 against MICM's table (:172-176); dense and sparse × Mozart and Doolittle. |
| 11 | DAE, conservation | F | `tests/test_external_constraints.c:76-98` (callbacks :8-29) | `:159` (`build/test_external_constraints`) | A transfer at k = 0.1 and the custom constraint y0 + y1 − 1 in the second row, through the solver's callback and custom-constraint interfaces; the four-stage DA set; 0.9 and 0.1; 20 solves of 10 s, each OK with the sum within 1e-4 (:90-97), dense and sparse, for three initial splits. |
| 12 | User-defined tutorial | F | `examples/rate_constants.c` | `:66` (`tools/test_tutorial_programs.py`) | MICM's own tutorial configuration, generated from MICM's tree (`Makefile:149-150`); MICM's run (`examples/rate_constants.c:41-64`); stdout compared byte for byte with `tests/reference/tutorial_custom.expected`, whose rows equal MICM's documented table (`tools/test_tutorial_programs.py:18-44`). The id shares the script's exit status with the seven other tutorial ids. |

Faithful: 8 of 12 (2, 5, 7, 8, 9, 10, 11, 12).

## Summary

| # | Case | terra-max | terra-xhigh | terra-high |
|---:|---|---|---|---|
| 1 | Robertson | F | F | P |
| 2 | RejectedStepAlpha | P | P | F |
| 3 | LU, random, vector | P | F | S |
| 4 | Markowitz reorder | P | P | P |
| 5 | Chapman, factored αI − J | F | F | F |
| 6 | Regression, vector Jacobian | P | P | P |
| 7 | Troe, all arguments | F | F | F |
| 8 | Surface, missing diffusion | F | F | F |
| 9 | State, invalid concentrations | P | F | F |
| 10 | Backward Euler, Oregonator | F | F | F |
| 11 | DAE, conservation | F | F | F |
| 12 | User-defined tutorial | F | F | F |
| | **Faithful of 12** | **7** | **9** | **8** |

How each runner turns checks into the 680 lines: terra-max and terra-xhigh print one line per id
from a check run for that id, by the program that runs it (some checks serve several ids);
terra-high looks each id up in a table and reports the mapped program's exit status.
