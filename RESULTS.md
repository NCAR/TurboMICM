# Validation status

All reference case IDs are covered. Final benchmark and portability checks are in progress.

| Check (`-O3 -DNDEBUG`) | Apple clang 21 | GCC 15.2 |
|---|---:|---:|
| Reference cases passing | 675 | 675 |
| Failing reference cases | 0 | 0 |
| C++-only compile-time cases (N/A) | 5 | 5 |
| TS1, 100 cells, 100 x 30 s, group size 31, 2 threads | PASS, bit-identical | PASS, max relative 3.27e-14 |

All eight tutorial programs match the reference's printed numerical results and
statistics. The configuration tutorial's elapsed time and OpenMP's interleaved
startup logging are nondeterministic and excluded from textual equality. The
rate tutorials and Robertson configuration additionally compare generated and
C-built mechanisms at full precision: rate constants, forcing and Jacobians
(1e-13 relative), and integrated states (1e-12 relative).

A full-size development run used 10,000 cells, 100 timed steps of 30 s after one
warm-up, input seed 2026, all available threads, and identical `-O3 -DNDEBUG`
flags for C and the clang++ MICM reference. C used groups of 128 cells.

| Implementation | elapsed_ms |
|---|---:|
| Portable C | 8903.025 |
| MICM reference | 25381.683 |

All 2,110,000 final species values were bit-identical. These are single local
measurements, not the evaluation's median-of-five ranking. The final-state gate
uses `1e-30 + 1e-10 * abs(reference)` for every cell and species.
