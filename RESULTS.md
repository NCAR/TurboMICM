# Results

## Correctness gates

`make CC=clang test-cpu` and `make CC=gcc-15 test-cpu` each pass 671 of the
680 MICM CPU and tutorial cases. Nine cases are marked N/A because they exercise
C++ compile-time traits or CUDA-only behavior; neither compiler reports a
failure.

The 10,000-cell TS1 comparison uses seed 2026, 100 intervals of 30 seconds,
and `-O3 -DNDEBUG` for both implementations. With both one and eight threads,
the strict comparison checked 2,110,000 values (10,000 cells by 209 species)
and passed; the largest relative difference was `2.78e-13` for BCARY in cell
4179. The comparison tolerance is `1e-30 + 1e-10 * |MICM|`.

## Benchmark

Single-run elapsed times from this machine are shown below. C and MICM reported
the same timed solver statistics at each thread count. With one thread they
reported 900 internal and accepted steps, zero rejected steps, 1,800 function
calls, 900 Jacobian updates and decompositions, and 2,700 linear solves. With
eight threads, the aggregate counts were 7,200 accepted steps, zero rejections,
14,400 function calls, 7,200 Jacobian updates and decompositions, and 21,600
linear solves.

| Threads | C elapsed (ms) | MICM elapsed (ms) | C / MICM |
|---:|---:|---:|---:|
| 1 | 185,107.482 | 145,414.379 | 1.273x |
| 8 | 34,479.378 | 31,198.866 | 1.105x |

These are single measurements, not medians. MICM completed this benchmark
faster on this machine; the experiment's ranking uses separate repeated
measurements on an idle machine.
