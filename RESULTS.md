# Results

## Correctness gates

With both Clang 21 and GCC 15, `make test-cpu` passes 671 of the 680 Tier 1
MICM CPU and tutorial cases. Nine cases are `N/A` because they cover C++
compile-time traits or CUDA-only behavior; neither compiler reports a failure.
With both compilers, `make test-specialized` passes the generated-versus-generic
solver comparison for TS1, Robertson, and both tutorial configurations.

The 100-cell TS1 comparison checked 21,100 values and passed for both
compilers. The largest relative difference was `3.14e-14` with Clang and
`4.61e-14` with GCC. The 10,000-cell, 18-thread Clang comparison checked
2,110,000 values and passed at `1e-10` relative tolerance; the largest relative
difference was `2.78e-13` for BCARY in cell 4179.

## TS1 benchmark

Clang 21, `-O3 -DNDEBUG`, seed 2026 from
`../reference/make_cells_csv.py`, 10,000 cells, 100 steps of 30 seconds. The
benchmark performs one untimed warm-up step. The 18-thread figures are medians
of three runs; one-thread figures are single runs.

| Threads | Generated solver | Tier 1 C | MICM | Generated / C | Generated / MICM |
|---:|---:|---:|---:|---:|---:|
| 18 | 5,307 ms | 15,063 ms | 12,479 ms | 2.84x faster | 2.35x faster |
| 1 | 52,212 ms | 183,852 ms | 134,931 ms | 3.52x faster | 2.58x faster |

All three implementations report 16,200 internal steps on 18 threads and 900
on one thread. The measured specialized-kernel ablations are in
[`SPECIALIZATION.md`](SPECIALIZATION.md).

The clean `make bench` build took `2.22 s` with Clang and `9.66 s` with GCC 15.
The generated TS1 C source was 1,453,684 bytes (28,434 lines). Both builds use
the Makefile's default `-O3 -DNDEBUG` flags.

## Day-night scenarios

The five disclosed scenarios ran with 18 threads, 32-cell groups, relative
tolerance `1e-6`, and absolute tolerance `1e-15`. All five pass the
`score.py --level 10` gate with Clang and GCC. Every run has zero NaN or
infinite values and zero values below `-1e-20`. The table shows the final
Clang run's time and score against the supplied MICM result, then its agreement
with the supplied high-accuracy reference solution.

| Scenario | Elapsed | MICM agreement, overall (worst species) | Reference agreement, overall (worst species) |
|---|---:|---:|---:|
| diurnal-30min | 999 ms | 14.37 (ISOP 13.80) | 4.748 (ISOP 3.610) |
| diurnal-10min | 1,935 ms | 14.29 (ISOP 13.69) | 5.052 (ISOP 3.824) |
| diurnal-5min | 3,131 ms | 14.23 (ISOP 13.68) | 5.214 (ISOP 4.328) |
| polluted-10min | 2,227 ms | 14.04 (ISOP 13.42) | 5.398 (ISOP 4.257) |
| clean-10min | 1,809 ms | 14.35 (ISOP 13.81) | 4.907 (ISOP 3.766) |

The reference scores show the expected difference between the MICM gate
tolerances and the supplied higher-accuracy reference. The specialized solver
tracks MICM within the required 10-digit gate on every scenario.

An additional 48-hour winter case used 1,000 cells, 900-second steps, 65°N,
day 355, a 01:00 start, and doubled NOx. Against a fresh MICM run with the same
32-cell groups and tolerances, it passed at 14.43 digits overall and 13.85 for
the worst species (ISOP), with no invalid or below-floor values.
