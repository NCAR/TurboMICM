# Verification results

Flags: `-O3 -DNDEBUG`, C11, OpenMP; no architecture or floating-point tuning.
Verified on the supplied Apple M5 system with Apple Clang 21 and GCC 15.2.

| Check | Clang | GCC 15 |
|---|---:|---:|
| Applicable reference cases | 675 PASS | 675 PASS |
| Failures | 0 | 0 |
| C++ compile-time-only cases | 5 N/A | 5 N/A |
| TS1, 100 cells, seed 1, 100 steps, 30 s | PASS | PASS |
| TS1, 10,000 cells, seed 2026, 100 steps, 30 s | PASS | PASS |
| Largest relative difference in the 10,000-cell gate | 2.03e-13 | 2.14e-13 |

Both `make test-cpu CC=clang` and `make test-cpu CC=gcc-15` exit successfully.
Each prints all 680 supplied case IDs exactly once. The five N/A cases contain
only C++ template/concept assertions; individual reasons are retained in
`tests/not_applicable.json`. All numerical and runtime behavior cases are ported.
The full builds produce no compiler warnings.

Coverage includes the rate formulas, constraints and initialization, all five
Rosenbrock methods, backward Euler, analytical and conservation DAE mechanisms,
all dense/sparse/grouped views, both LU algorithms and separate factors, state
ownership and named updates, system metadata, parameterized species, diagnostics,
solver construction, aerosol/external models, Chapman and flow-tube regressions,
and OpenMP shared-mechanism execution. The tests retain the source quantities,
inputs, output times, and tolerance bounds. Independent test-only Python provides
the legacy Chapman comparison; no test execution loads MICM or MUSICA.

All eight plain-program tutorials run through C library code and match retained
reference concentration tables and solver statistics. Comparisons exclude only
measured timing values and the unsynchronized OpenMP thread startup messages.
The eight-thread unit case requires bit-identical results across independent
states sharing one mechanism.

Generated YAML/JSON mechanisms match independently assembled mechanisms through
ten integration steps. The checks cover every supplied configuration example,
condition-dependent species, Robertson, and deterministic TS1 translation.
The TS1 gate compares all 2,110,000 temperature, pressure, and species values in
the 10,000-cell result with `atol=1e-30`, `rtol=1e-10`. Both final compiler runs
pass `../reference/compare_final.py` against the same retained MICM result.

## Local timing

Input: 10,000 cells, seed 2026, 100 timed steps of 30 s after one warm-up.
These are individual local observations, not the experimenter's median-of-five
ranking measurements. Reference and C programs use the same OPTFLAGS.

| Program | Threads/groups | elapsed_ms |
|---|---|---:|
| C, Clang | OpenMP default (10 threads), groups of 32 | 8,083.676 |
| C, GCC 15 | OpenMP default (10 threads), groups of 32 | 9,634.764 |
| MICM reference, Clang++ | 10 threads/groups | 27,409.392 |
| MICM supplied baseline, Clang++ | 1 thread | 137,383 |

The final Clang observation is approximately 3.39 times faster than the local
10-thread MICM reference. The one-thread MICM baseline is the supplied median,
not a new local measurement.

Reference statistics for 10 groups: 9,000 internal steps, all accepted, 18,000
function calls, 9,000 Jacobian updates/decompositions, 27,000 solves. C statistics
for 313 groups: 281,700 internal steps, all accepted, 563,400 function calls,
281,700 Jacobian updates/decompositions, 845,100 solves. Both use nine accepted
internal steps per group per public solve. Statistics exclude the warm-up.

The source review found only ISO C11, standard library/libm calls, and CPU OpenMP
in production C. No intrinsics, target attributes, relaxed floating-point flags,
architecture tuning, third-party numerical libraries, or Python runtime links
are used. Linux/x86-64 execution has not been measured in this environment.
