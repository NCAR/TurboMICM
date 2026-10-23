# Mechanism specialization

The generator consumes CAMP JSON or YAML, including `camp-files`, using the Tier
1 reader. It is independent of species names and reaction counts. TS1 is an input
configuration, not a branch in the generator. Its 209 dynamic species and 547
reactions produce a 2,381-entry filled matrix and 6,986 elimination updates.

Generation fixes species ordering, sparse coordinates and fill, derivative jobs,
reaction arithmetic, stoichiometric yields, rate families and fixed coefficients,
alpha shifts, elimination operations, and triangular-solve rows. The selected
implementation emits direct reaction/Jacobian code and immutable compact LU
operation tables. The tables share short SIMD kernels, reducing instruction size
compared with full unrolling. Four lower updates share each upper-factor load.

Temperature, pressure, air density, concentrations, custom rate parameters, step
length, and tolerances remain data. Generated rate formulas read current values
at every update. No Python, configuration parser, or third-party runtime is linked
into the C library or production benchmarks.

The generated solver instantiates the shared Rosenbrock algorithm with direct
kernel hooks and a constant species count. The parameter tables, error estimate,
acceptance/rejection decisions, and step-size controller remain MICM's. Three-stage
vector passes are fused while retaining each element's operation order. Reciprocal
pivots replace repeated diagonal divisions with double-precision multiplies.
Normalized error components are computed with SIMD in completed stage storage;
the sum of squares keeps the original strict species/cell order and expression.
The stage storage is rebuilt before every subsequent trial, including retries.

All alternatives below were generated and built with Clang and `-O3 -DNDEBUG`.
Measurements ran sequentially. Each entry is the median of three runs, with one
untimed warm-up, 30-second steps, and the same seed-2026 cells. Every measured
variant passed the final-state comparison. These are local measurements, not the
experimenter's ranking results. Rows are complete combinations; differences
between rows that change several controls do not isolate a single cause.

| Variant | 1 thread, 1,000 cells × 20 steps (ms) | 18 threads, 10,000 cells × 100 steps (ms) |
|---|---:|---:|
| Generic Tier 1 | 895.678 | 3668.458 |
| Direct chemistry and fully unrolled LU, dynamic width | 778.693 | — |
| Same, 16 operations per cell loop | 1388.108 | — |
| Unrolled kernels, constant width 32 | 1003.911 | — |
| Same, reciprocal pivots | 1056.532 | — |
| Grouped workspace, 2-cell tiles, fused stages | 2609.205 | — |
| Grouped workspace, 4-cell tiles, fused stages | 1568.978 | — |
| Grouped workspace, 8-cell tiles, fused stages | 1698.650 | — |
| Packed 2-cell matrix tiles, compact LU | 1501.294 | — |
| Direct chemistry, compact LU, one update per loop | 683.765 | 2926.636 |
| Same, shared reaction-shape kernels | 709.832 | 3240.083 |
| Direct chemistry, compact LU, four updates per loop | 673.311 | 3117.265 |
| Same, SIMD error components and strict scalar sum | 685.188 | 2690.342 |
| Same, narrow plan indices | 710.061 | 2768.705 |
| Narrow indices, SIMD errors, one LU update per loop | 727.745 | 3012.490 |

The selected defaults are direct chemistry, compact LU with four updates,
32-cell full groups, no cell tiling, `size_t` plan entries, reciprocal pivots,
fused stages, and SIMD error components. Smaller code alone did not guarantee
faster execution: reaction-shape sharing reduced kernel duplication but added table
lookups, and narrow indices did not improve the selected combination. Full LU
unrolling and small cell tiles were slower on this workload. They remain optional generator controls, with short
and partial-group coverage in `make test-specialized`.

A symbolic ordering experiment compared stable Markowitz ties with randomized
ties and greedy minimum fill. The best symbolic alternative reduced elimination
updates from 6,986 to 6,314 and entries from 2,381 to 2,293. It was not adopted or
timed: the implementation retains the reader's established stable ordering, and
no runtime speedup is attributed to this symbolic count reduction.

Generated C is 990,603 bytes with the selected controls, identical under both
compilers. Clean `make -j4 bench` took 13.77 seconds with Clang and 14.89 seconds
with GCC, excluding environment creation. README.md describes reproduction;
DESIGN.md lists every retained control. RESULTS.md records the final full-workload
measurements and all scenario gates.

A final interleaved ablation changed one control at a time. On 18 threads, the
selected solver took 2746.376 ms, disabling stage fusion took 2919.862 ms, and
ordinary LU diagonals took 2881.896 ms (medians of three). On the 1,000-cell,
one-thread workload, the respective times were 676.137, 690.752, and 674.274 ms.
Reciprocal pivots therefore helped the primary workload without a measurable
one-thread gain. Every ablation still passed the MICM final-state comparison.
