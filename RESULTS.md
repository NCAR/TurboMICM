# Validation results

All 680 reference cases and supplemental checks pass on both compilers. The
full TS1 gate and all additional seed/group/schedule comparisons pass.

## Case results

Both compiler builds use `OPTFLAGS='-O3 -DNDEBUG'` and the Makefile's corresponding
OpenMP flags. The ported cases keep original or tighter numerical bounds.

| Compiler | PASS | FAIL | N/A |
| --- | ---: | ---: | ---: |
| Apple Clang 21 | 680 | 0 | 0 |
| GCC 15.2 | 680 | 0 | 0 |

The 680 cases comprise 31 rate-constant tests, 11 rate-store cases, 22 state
cases, 32 LU decomposition
cases, 20 linear-solver cases, 2 backward Euler cases, 42 constraint
definition/set cases, 12 constraint initialization cases, 12 DAE integration
cases, 8 Rosenbrock utility/regression cases, 10 process cases and 7 Jacobian
verification cases, all 117 dense/vector matrix cases, 14 species/phase/system
cases, 5 version cases, callback-rate override, OpenMP state-isolation, and all
190 sparse matrix cases, 27 grouped-view/view-category cases, 20 external-model
constraint cases, 18 aerosol integration cases and two Terminator integration
cases, 13 Chapman unit cases, one Chapman integration case, 12 Chapman
regression cases, and two flow-tube regression cases, plus all 32 analytical
integration cases, eight solver-builder cases and
eight complete tutorial programs. Tutorial comparisons preserve every printed
result and counter; measured durations and thread announcement interleaving are
validated and normalized as documented in DECISIONS.md. Additional
checks pass on both compilers: all three tutorial configurations versus
independent hand-built mechanisms (one and seven cells, ten integration
intervals), 49 configuration validation checks, independent formulas for three
additional generated rate types, all 379 TS1 initial fields by name, and a TS1
interval-clamp regression. Builds have no compiler warnings.

## TS1 final-state gate

MICM was rebuilt with `make micm-ref` using Clang++ and the same optimization
flags. All comparisons below use the required absolute tolerance `1e-30` and
relative tolerance `1e-10`, through the supplied `compare_final.py`.

| Input and run | C configuration | Checked values | Largest relative difference | Result |
| --- | --- | ---: | ---: | --- |
| 100 cells, seed 1, 100 steps, 30 s | Clang, 1 thread, groups of 32 | 21,100 | 9.7e-13 | PASS |
| Same | GCC, 1 thread, groups of 32 | 21,100 | 9.7e-13 | PASS |
| Same | GCC, 4 threads, groups of 7 | 21,100 | 9.7e-13 | PASS |
| 10,000 cells, seed 2026, 100 steps, 30 s | Clang, 10 threads, groups of 32 | 2,110,000 | 0 | PASS |
| Same | GCC, 10 threads, groups of 32 | 2,110,000 | 1.43e-13 | PASS |

Final extended runs also checked the following settings on both compilers,
all for 100 steps of 30 seconds after warm-up:

| Cells / seed | Group / threads / OpenMP schedule | Clang max relative difference | GCC max relative difference |
| --- | --- | ---: | ---: |
| 127 / 17 | 1 / 1 / static | 0 | 4e-14 |
| 127 / 17 | 7 / 4 / dynamic,1 | 0 | 4e-14 |
| 127 / 17 | 32 / 10 / guided,2 | 0 | 4e-14 |
| 127 / 17 | 1000 / 3 / static,2 | 0 | 4e-14 |
| 129 / 73 | 3 / 10 / static,1 | 0 | 3.92e-14 |
| 129 / 73 | 64 / 4 / dynamic,3 | 0 | 3.92e-14 |
| 129 / 73 | 129 / 1 / static | 0 | 3.92e-14 |

A group request larger than the cell count is capped to that count. The 127-cell
runs compare 26,797 fields each; the 129-cell runs compare 27,219 fields each.

The 10,000-cell MICM reference used ten threads, one group per thread. The C
benchmark used 313 groups, each with independent adaptive solver state. All
runs include an untimed warm-up interval before the timed intervals.

The retained one-cell TS1 regression checks the public solver's post-interval
nonnegative concentration clamp. It catches a difference that the 100-cell
sample alone did not expose.

## Diagnostic timings

Single runs on the supplied laptop; these are development measurements, not
idle-machine medians for ranking. No compiler optimization flags were tuned.

| Run | Clang C | GCC C | MICM Clang++ |
| --- | ---: | ---: | ---: |
| 100 cells, seed 1, 100 steps, 30 s, one thread | 516.17 ms | — | 1,316.91 ms |
| 10,000 cells, seed 2026, 100 steps, 30 s, ten threads | 9,625.25 ms | 10,602.78 ms | 27,421.51 ms |

The final GCC validation separately measured MICM at 26,035.22 ms. Clang C was
2.85 times faster than its paired MICM run. These measurements are single runs;
the evaluation harness determines the median on an idle machine.

For the full-size C run: 281,700 internal/accepted steps, zero rejected steps,
563,400 function calls, 281,700 Jacobian updates and decompositions, and 845,100
linear solves. MICM's larger groups report 9,000 internal/accepted steps,
18,000 function calls, 9,000 Jacobian updates and decompositions, and 27,000
solves. These counts sum over groups, not cells.

The final source audit found no forbidden extensions, target flags, extended
floating types or operating-system APIs in the C implementation. The case list
matches the supplied inventory exactly, and the TS1 translation reproduces the
committed mechanism. All C and generated C use the Makefile compiler and flags.
Linux execution was not available in this environment; the portability audit
and both compiler runs were performed on the supplied macOS machine.

Reproduce the extended numerical comparisons with `tools/validate_benchmarks.py`
as described in README.md. Logs and full comparison files are retained under
`build/validation/final-clang/` and `build/validation/final-gcc/`.
