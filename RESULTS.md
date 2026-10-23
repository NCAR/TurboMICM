# Results

## Tier 1 regression and specialized checks

On this M5 Pro, `make test-cpu` passes with Clang 21 and GCC 15. Each run
reports the full 680-case inventory: 662 pass and 18 C++-only cases are marked
N/A. The Tier 1 generic benchmark continues to build with both compilers.

`make test-specialized` passes with both compilers. It checks the TS1 solver
against the generic solver over 100 cells and 100 steps, and checks generated
forcing, negative Jacobians, and complete short solver runs against the generic
path for the rate-constant YAML and user-defined-rate JSON tutorial mechanisms.
The TS1 specialized final state also passes the MICM final-state comparator for
the 100-cell seed-1 case; the largest relative difference was `1.32e-14`.

## TS1 benchmark

Clang 21, `-O3 -DNDEBUG`, 10,000 cells, 100 steps, `DT=30 s`, seed 2026. These
are local single-run measurements, not the experimenter's median-of-five
ranking measurements. `elapsed_ms` is the program's timed region; wall time
includes process startup and output.

| Threads | Generic Tier 1 | Specialized | MICM |
| ---: | ---: | ---: | ---: |
| 18 | 21,337 ms | 7,762 ms | 12,736 ms |
| 1 | 305,344 ms | 90,423 ms | 134,178 ms |

All three programs reported the same work: 900 internal steps, 900 accepted,
zero rejected, 1,800 function calls, 900 Jacobian updates and decompositions,
and 2,700 solves on one thread. At 18 threads the corresponding totals were
16,200, 16,200, zero, 32,400, 16,200, 16,200, and 48,600. The generated C file
is 1,124,511 bytes. The final clean `make bench` build took 11.84 s with Clang.

The clean `make bench` build took 30.83 s with GCC 15. A clean GCC verification
sequence passed both `make test-cpu` and `make test-specialized`.

## Day and night scenarios

The table reports the final Clang 21 specialized runner with 32-cell groups,
`--rtol 1e-6 --atol 1e-15`. Each row passed `score.py --level 10` against the
provided MICM snapshots, with no NaN, infinity, or value below `-1e-20`. GCC 15
also passed all five disclosed scenarios with the same settings.

| Scenario | MICM agreement, overall digits | Worst species digits | Reference agreement, overall digits |
| --- | ---: | ---: | ---: |
| diurnal-30min | 14.37 | 13.82 | 4.748 |
| diurnal-10min | 14.29 | 13.69 | 5.052 |
| diurnal-5min | 14.24 | 13.69 | 5.214 |
| polluted-10min | 14.08 | 13.46 | 5.398 |
| clean-10min | 14.34 | 13.80 | 4.907 |

Reference scores are for scale, because the supplied high-accuracy solution is
not the gate reference. The final generated LU schedule was used for all
reported Clang scores. No undisclosed scenarios were available for local
verification.

## Build and portability

The generated source is recreated from `mechanisms/ts1.json` by `make`; no
generated C is hand-edited. C targets use `$(CC)` and `$(OPTFLAGS)`, ISO C11,
libm, and CPU OpenMP. The compiler-specific OpenMP include and link flags are
selected by the Makefile. Clang and GCC builds both pass the generic and
specialized test targets.
