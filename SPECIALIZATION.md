# Specialization

## What is fixed

For a generated mechanism, the C source fixes species and reaction indexes,
reactant multiplicity, product yields, rate-constant parameters, third-body
order, species ordering, Jacobian sparsity, symbolic fill, and each LU and
triangular-solve operation. The reaction and matrix loops become individual C
statements. TS1's diagonal Markowitz order reduces its generated sparse LU
pattern to 2,360 positions for 209 state species.

The code keeps temperature, pressure, air density, concentrations, custom rate
parameters, step size, and tolerances as run-time values. The rate constants
are computed when the solver is created because temperature and pressure stay
fixed during a chemistry run. A scenario can update user-defined rates between
steps; photolysis parameters remain cell data scaled by the factor at the
current step.

The Rosenbrock stages, coefficients, error estimate, tolerances, controller,
acceptance rules, and state clamping continue to use the shared solver. The
generated linear algebra is exact no-pivot LU with the same `fma` update order
as the generic implementation. Constraints select the generic path when they
change the sparse pattern.

## Benchmark

Clang 21, `-O3 -DNDEBUG`, 10,000 cells, 100 steps of 30 seconds, seed 2026,
default solver tolerances, 18 OpenMP threads. Each implementation uses one
untimed warm-up step. The 18-thread figures are medians of three runs; the
one-thread figures are single runs.

| Threads | Generated solver | Tier 1 C | MICM | Generated / C | Generated / MICM |
|---:|---:|---:|---:|---:|---:|
| 18 | 5,307 ms | 15,063 ms | 12,479 ms | 2.84x faster | 2.35x faster |
| 1 | 52,212 ms | 183,852 ms | 134,931 ms | 3.52x faster | 2.58x faster |

All three report 16,200 internal steps on 18 threads and 900 on one thread.
The generated 18-thread final state passed `compare_final.py` against MICM for
all 2,110,000 values; the largest relative difference was `2.78e-13`.

## Kernel ablations

The generator can omit one generated callback at a time with its
`--without-*-kernel` options. The Makefile target `make bench-ablation` builds
those experimental variants. Each row below is one 18-thread run of the same
10,000-cell workload; the complete generated solver's 5,307 ms is the median
above.

| Disabled generated kernel | Elapsed | Change from complete solver |
|---|---:|---:|
| None | 5,281 ms | baseline run |
| Rate setup and custom-rate update | 5,375 ms | within run-to-run variation |
| Reaction and Jacobian updates | 8,022 ms | 51% slower |
| Sparse LU factorization | 8,560 ms | 61% slower |
| Triangular solves | 7,742 ms | 46% slower |

The rate kernel is used while the solver is created, before the benchmark's
timed region. Its disablement therefore does not measurably affect the primary
benchmark, which matches the timing contract. The same kernel also updates
custom rates inside a timed scenario step. The reaction, factor, and solve
kernels run for every internal solver step and account for the measured
improvement. No tested specialization was dropped; `make bench` emits all four
kernels.
