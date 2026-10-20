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
comparison utility with both compilers. A broader run over all 100 supplied
cells also passed:

| Compiler | C elapsed time | MICM elapsed time | Solver steps | Largest relative difference |
| --- | ---: | ---: | ---: | ---: |
| Clang | 235.379 ms | 10.430 ms | 900 | `2.24e-15` |
| GCC 15 | 322.778 ms | 11.009 ms | 900 | `2.42e-15` |

For 100 seed-1 cells, 100 steps, and four threads, Clang completed in
6,113.087 ms; the comparison run took 338.929 ms. All 21,100 state values
passed, with a largest relative difference of `9.7e-13`.

These timings are not the requested 10,000-cell performance run. The current C
implementation is substantially slower than the reference because it uses
per-cell integration and scalar sparse-matrix work. Workspaces and symbolic
fill are reused between calls. The 10,000-cell, 100-step benchmark and all 680
case ports remain incomplete.
