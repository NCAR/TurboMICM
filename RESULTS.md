# Results

## Build and configuration checks

- `make cpu CC=clang`: passed on the provided macOS machine.
- `make cpu CC=gcc-15`: passed on the provided macOS machine.
- `tests/test_core.c`: passed with Clang and GCC 15, including rate formula
  checks across Arrhenius, branched, lambda, reversible, surface, Taylor, Troe,
  ternary activation, tunneling, and user-defined rates; also forcing and small
  analytic integration checks.
- JSON and YAML tutorial mechanisms normalize to 7 species; the input reactions
  expand to 10 and 7 generated reactions, respectively (JSON and YAML).
- TS1 import contains 210 phase species, 547 reactions, and the expected initial
  defaults. The generated mechanism compiles with both compilers.
- `make test-cpu CC=clang` and `make test-cpu CC=gcc-15`: both emit all 680 case
  statuses. Thirty-one rate-constant cases pass; 649 cases still report
  incomplete ports, so both targets exit nonzero.

## TS1 final-state comparison

The seed-1 first cell, 100 steps at `DT = 30 s`, passed the supplied MICM
comparison utility with both compilers. The 100-cell seed-1 run also passed:

| Compiler | C elapsed time | MICM elapsed time | Solver steps | Largest relative difference |
| --- | ---: | ---: | ---: | ---: |
| Clang | 142.861 ms | 10.430 ms | 900 | `1.13e-14` |
| GCC 15 | 165.924 ms | 11.009 ms | 900 | `2.58e-15` |

For 100 seed-1 cells, 100 steps, and four threads, Clang completed in
4,416.783 ms; the comparison run took 347.269 ms. All 21,100 state values
passed, with a largest relative difference of `9.70e-13`.

The 10,000-cell seed-2026 run, 100 steps, and 10 threads completed in
300,799.575 ms; the reference completed in 27,457.982 ms. The comparison failed
2,185 of 2,110,000 values across 26 cells; the largest relative difference was
`5.33e-4` for HO2NO2 in cell 3974. A single-cell run of that input also fails,
so the discrepancy is not caused by cell grouping. This result is a regression
case for further solver arithmetic work. The implementation is about ten times
slower than the reference at this scale, and 649 of the 680 inventory cases
remain unported.
