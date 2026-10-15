# Verification in progress

Flags: `-O3 -DNDEBUG`, C11, OpenMP; no architecture or floating-point tuning.

| Check | Clang | GCC 15 |
|---|---:|---:|
| Ported reference cases | 593 PASS | 593 PASS |
| Cases still awaiting ports | 82 | 82 |
| C++ compile-time-only cases | 5 N/A | 5 N/A |
| TS1, 100 cells, seed 1, 100 steps, 30 s | PASS | PASS |
| TS1, 10,000 cells, seed 2026, 100 steps, 30 s | PASS | PASS |
| Largest relative difference in the 10,000-cell gate | 2.03e-13 | 2.14e-13 |

The full 680-case gate is **not complete**. All pending cases report FAIL and
make test-cpu exits nonzero. Current passes cover rate constants, linear and
equilibrium constraints, constraint sets, Rosenbrock unit quantities, constraint
initialization, backward Euler basics, callback replacement, two tutorials, all
117 dense/vector-matrix cases, all 190 sparse-matrix cases, the 21 runtime grouped
view cases, and runtime vector-view use. Five compile-time-only cases have
individual N/A reasons in tests/not_applicable.json.

Coverage also includes all 52 LU/linear-solver cases, 22 state cases, 14
species/phase/system cases, five version cases, 11 reaction-rate-store cases,
four process cases, five process-set cases, and seven Jacobian-verification
cases. LU reconstruction and solve residuals retain the original absolute or
four-ULP bounds. Both compilers pass the condition-parameterization configuration
equivalence check. Every outstanding FAIL is explicitly a pending case port.

All 32 analytical integration cases now pass: the six rate-dependent chains and
their fast-equilibrium variants, surface reactions, Robertson, HIRES, E5, and
Oregonator. Method, LU, grouping, time intervals, reference concentrations, and
tolerances are retained in the test fixtures. Both Terminator cases and all 12
equilibrium/conservation DAE integration cases pass, including algebraic
tolerance sensitivity, nonnegativity, coupled constraints, non-unit
stoichiometry, state reordering, and the combined chlorine/Robertson mechanism.

Both rate tutorial programs print exactly the retained reference output.
Generated YAML/JSON mechanisms match independently assembled mechanisms through
ten integration steps. Robertson configuration equivalence and deterministic TS1
translation are also checked. The TS1 gate uses every temperature, pressure, and
species value with `atol=1e-30`, `rtol=1e-10`.

## Preliminary timing

Input: 10,000 cells, seed 2026, 100 timed steps of 30 s after one warm-up.
These are individual local observations, not the experimenter's median-of-five
ranking measurements. Reference and C programs were built with the same OPTFLAGS.

| Program | Threads/groups | elapsed_ms |
|---|---|---:|
| C, Clang | OpenMP default, groups of 32 | 8,167.059 |
| C, GCC 15 | OpenMP default, groups of 32 | 9,699.149 |
| MICM reference, Clang++ | 10 threads/groups | 27,409.392 |
| MICM supplied baseline, Clang++ | 1 thread | 137,383 |

Reference statistics for 10 groups: 9,000 internal steps, all accepted, 18,000
function calls, 9,000 Jacobian updates/decompositions, 27,000 solves. C statistics
for 313 groups: 281,700 internal steps, all accepted, 563,400 function calls,
281,700 Jacobian updates/decompositions, 845,100 solves. Both use nine accepted
internal steps per group per public solve. Statistics exclude the warm-up.

## Remaining work

Port the remaining Chapman, regression, aerosol/external-model,
integrated reaction-rate, solver-builder, OpenMP, and tutorial cases; address
any missing API behavior they identify. Then run the complete suite under both
compilers, finish portability review, repeat representative benchmark gates, and
replace this partial report with the completed evidence.
