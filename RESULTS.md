# Validation status

Implementation has started from an empty repository. No correctness or performance
claim has yet been established. Required validation includes all 680 reference
cases, the complete TS1 reference comparison, both requested compilers, tutorial
configuration equivalence, and benchmark timing with matching compiler flags.

Initial clang validation (`-O3 -DNDEBUG`): the C TS1 benchmark with one thread,
100 cells, 100 timed steps of 30 s passes the supplied full-state comparison
(21,100 values; largest relative difference 9.7e-13). Timed integration: 608 ms.
This is a development measurement, not a ranking claim. The 680-case suite and
GCC validation remain outstanding.

Current development checkpoint:

| Check | clang | gcc-15 |
|---|---:|---:|
| Implemented reference cases passing | 140 | 140 |
| Remaining reference cases (reported FAIL) | 540 | 540 |
| N/A | 0 | 0 |
| TS1, 100 cells, 100 x 30 s, groups of 31, 2 threads | PASS (bit-identical) | PASS (max relative 3.27e-14) |

All eight tutorial programs match the reference's printed numerical results and
statistics. The configuration tutorial's elapsed time and OpenMP's interleaved
startup logging are nondeterministic and excluded from textual equality. The
rate tutorials and Robertson configuration additionally compare generated and
C-built mechanisms at full precision: rate constants, forcing and Jacobians
(1e-13 relative), and integrated states (1e-12 relative).
