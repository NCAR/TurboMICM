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
