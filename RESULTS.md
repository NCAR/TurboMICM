# Validation progress

The complete acceptance gate has **not** been met: 336 of the 680 case ports
remain pending. The test runner reports these as FAIL, not N/A.

## Current case results

Both compiler builds use `OPTFLAGS='-O3 -DNDEBUG'` and the Makefile's corresponding
OpenMP flags. The ported cases keep original or tighter numerical bounds.

| Compiler | PASS | Pending FAIL | N/A |
| --- | ---: | ---: | ---: |
| Apple Clang 21 | 344 | 336 | 0 |
| GCC 15.2 | 344 | 336 | 0 |

The 344 cases comprise 31 rate-constant tests, 11 rate-store cases, 22 state
cases, 32 LU decomposition
cases, 20 linear-solver cases, 2 backward Euler cases, 42 constraint
definition/set cases, 12 constraint initialization cases, 10 DAE integration
cases, 8 Rosenbrock utility/regression cases, 9 process cases and 7 Jacobian
verification cases, all 117 dense/vector matrix cases, 14 species/phase/system
cases, 5 version cases, callback-rate override and OpenMP state-isolation cases. Additional
checks pass on both compilers: all three tutorial configurations versus
independent hand-built mechanisms (one and seven cells, ten integration
intervals), and a TS1 interval-clamp regression.

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
| 10,000 cells, seed 2026, 100 steps, 30 s, ten threads | 8,624.21 ms | 10,293.17 ms | 27,188.29 ms |

For the full-size C run: 281,700 internal/accepted steps, zero rejected steps,
563,400 function calls, 281,700 Jacobian updates and decompositions, and 845,100
linear solves. MICM's larger groups report 9,000 internal/accepted steps,
18,000 function calls, 9,000 Jacobian updates and decompositions, and 27,000
solves. These counts sum over groups, not cells.

Remaining validation includes the full case inventory, external/aerosol models,
all solver regression and tutorial output checks, additional seeds and tuning
values, and a final portability audit. See WORK_LOG.md for the completion ledger.
