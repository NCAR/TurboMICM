# Results

## Build and configuration checks

- `make cpu CC=clang`: passed on the provided macOS machine.
- `make cpu CC=gcc-15`: passed on the provided macOS machine.
- `tests/test_core.c`: passed with Clang and GCC 15, including Arrhenius,
  forcing, and small analytic integration checks.
- JSON and YAML tutorial mechanisms normalize to 7 species; the input reactions
  expand to 10 and 7 generated reactions, respectively (JSON and YAML).
- TS1 import contains 210 phase species, 547 reactions, and the expected initial
  defaults. The generated mechanism compiles with both compilers.
- `make test-cpu`: emits all 680 case statuses. Two Arrhenius calculation cases
  pass; 678 cases still report incomplete ports, so the target exits nonzero.

## TS1 final-state comparison

The seed-1 first cell, 100 steps at `DT = 30 s`, passed the supplied MICM
comparison utility with both compilers:

| Compiler | C elapsed time | MICM elapsed time | Solver steps | Largest relative difference |
| --- | ---: | ---: | ---: | ---: |
| Clang | 247.144 ms | 10.383 ms | 900 | `2.24e-15` |
| GCC 15 | 347.589 ms | 10.636 ms | 900 | `2.42e-15` |

These timings are a one-cell correctness check, not the requested 10,000-cell
performance run. The current C implementation is substantially slower than the
reference because it uses per-cell integration and repeats dense workspace
allocation and sparse-pattern factor work for each solver call. The 10,000-cell,
100-step benchmark and all 680 case ports remain incomplete.
