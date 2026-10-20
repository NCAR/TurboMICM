# Results

## Build and configuration checks

- `make cpu CC=clang`: passed on the provided macOS machine.
- `make cpu CC=gcc-15`: passed on the provided macOS machine.
- `tests/test_core.c`: passed with Clang and GCC 15, including rate formula
  checks across Arrhenius, branched, lambda, reversible, surface, Taylor, Troe,
  ternary activation, tunneling, and user-defined rates; also forcing and small
  analytic integration checks, C reaction construction and value-copy checks,
  and density-parameterized reaction multipliers.
- `tests/test_codegen.py`: passed with both compilers through `make test-cpu`,
  checking rate-type ordering, parameter defaults and indexes, derived branched
  parameters, surface parameters, density multipliers, and reaction counts in
  generated C.
- The C core test ports the README tutorial's Foo/Bar/Baz mechanism and checks
  all ten printed states at the tutorial's three-significant-digit precision.
- The C core test ports `test_solver_results` with three grouped cells and
  checks the printed solver state and counters. Its 18 accepted steps, 36
  function calls, 18 Jacobian updates/decompositions, and 54 solves match MICM.
- The C core test also runs the grouped chemistry for ten 200-second intervals
  and checks all 90 printed concentrations from `test_multiple_grid_cells` at
  the tutorial's three-significant-digit precision.
- Repeated-reactant Jacobians now match for both separate unit terms (`B + B`)
  and one stoichiometric term (`2 B`), covering a coefficient case that the
  previous negative-Jacobian implementation mishandled.
- JSON and YAML tutorial mechanisms normalize to 7 species; the input reactions
  expand to 10 and 7 generated reactions, respectively (JSON and YAML).
- TS1 import contains 210 phase species, 547 reactions, and the expected initial
  defaults. The generated mechanism compiles with both compilers.
- `make test-cpu CC=clang` and `make test-cpu CC=gcc-15`: both emit all 680 case
  statuses. Forty-six cases pass; 634 cases still report incomplete ports,
  so both targets exit nonzero.

## TS1 final-state comparison

The seed-1 first cell, 100 steps at `DT = 30 s`, passed the supplied MICM
comparison utility with both compilers. The 100-cell seed-1 run also passed:

| Run | Compiler | C elapsed time | MICM elapsed time | Solver steps | Largest relative difference |
| --- | --- | ---: | ---: | ---: | ---: |
| First cell | Clang | 173.842 ms | 10.767 ms | 900 | `0` |
| First cell | GCC 15 | 165.693 ms | 10.767 ms | 900 | `5.69e-15` |
| 100 cells, 4 threads | Clang | 4,371.516 ms | 354.051 ms | 3,600 | `9.70e-13` |
| 100 cells, 4 threads | GCC 15 | 4,526.347 ms | 354.051 ms | 3,600 | `9.70e-13` |

The 10,000-cell seed-2026 run, 100 steps, and 10 threads completed in
300,799.575 ms; the reference completed in 27,457.982 ms. The comparison failed
2,185 of 2,110,000 values across 26 cells; the largest relative difference was
`5.33e-4` for HO2NO2 in cell 3974. A single-cell run of that input also fails,
so the discrepancy is not caused by cell grouping. This result is a regression
case for further solver arithmetic work. The implementation is about eleven
times slower than the reference at this scale, and 634 of the 680 inventory cases
remain unported.
