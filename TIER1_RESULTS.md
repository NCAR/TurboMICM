# Historical Tier 1 results

Retained from the initial repository. Fresh Tier 3 gates and measurements are
in RESULTS.md.

# Validation results

The CPU correctness gate passes with both requested compilers. Builds use
`-O3 -DNDEBUG`, ISO C11, pedantic diagnostics, and the specified OpenMP flags.

| Check | Apple clang 21 | GCC 15.2 |
|---|---:|---:|
| Reference case IDs | 680 | 680 |
| PASS | 675 | 675 |
| FAIL | 0 | 0 |
| N/A: C++ compile-time concepts/categories | 5 | 5 |
| TS1, 100 cells, 100 x 30 s, seed 1 | PASS, bit-identical | PASS, max relative 3.27e-14 |
| TS1, 10,000 cells, 100 x 30 s, seed 2026 | PASS, bit-identical | PASS, max relative 1.43e-13 |
| TS1, 1,000 cells, 100 x 30 s, seed 7415 | PASS, bit-identical | PASS, max relative 8.10e-14 |

The full-size comparison checks every one of 2,110,000 numeric CSV entries by
column name. The gate is `abs(y-reference) <= 1e-30 + 1e-10*abs(reference)`.
Clang full-size runs with group sizes 16, 32, 64, 128 and 256 are bit-identical
to MICM. Small-input checks also exercise group sizes 1, 3 and 31; GCC checks
1, 3, 31, 128 and 511 plus dynamic scheduling with groups of 7. Partial groups
are included. All comparisons pass.

All eight plain-program tutorials match deterministic printed numerical results
and statistics. Elapsed time and interleaved OpenMP startup logging are excluded
as nondeterministic. The three tutorial configurations additionally match
independently C-built mechanisms at full precision: rates, forcing and
Jacobians (1e-13 relative), and trajectories (1e-12 relative).

Ten reader validation tests cover includes/cycles, duplicate and unknown species,
non-finite numbers, reactant quantities, units/repeated terms, extra rate types,
C identifiers and missing surface properties. Malformed benchmark CSV and CLI
inputs are rejected. An additional clang build with
`-O0 -g -fsanitize=address,undefined` passes 675 cases without sanitizer findings.
The translated TS1 JSON is byte-identical when regenerated from the specified
source revision. Linked-library inspection shows only the system C runtime
and OpenMP for the clang benchmark; Python is absent.

## Local benchmark measurements

10,000 cells, seed 2026, one untimed warm-up then 100 steps of 30 s, all ten CPU
threads. The C benchmark and clang++ MICM reference use identical
`-O3 -DNDEBUG` flags. CSV input/output and setup are outside the timed region.

| Implementation / group size | elapsed_ms |
|---|---:|
| C clang, 16 | 8322.540 |
| C clang, 32 (selected default) | 7714.497; earlier run 8020.927 |
| C clang, 64 | 8024.765; repeat 8393.091 |
| C clang, 128 | 8990.938 |
| C clang, 256 | 13263.903 |
| C GCC, 32 | 8968.018 |
| MICM clang++, ten contiguous groups | 26480.626 |

The selected default is about 3.4 times faster than this MICM measurement.
These are local development timings, not the evaluator's median-of-five ranking.
Group size remains configurable for other machines. The default groups shrink
for small inputs to occupy available threads.

With groups of 32 the timed C statistics are: 281700 internal/accepted steps,
zero rejected steps, 563400 function calls, 281700 Jacobian updates and
decompositions, and 845100 solves. MICM's ten groups report respectively
9000, 0, 18000, 9000 and 27000. Counts differ because statistics are summed over
independent cell groups; every final concentration agrees.

## N/A cases

The following source cases contain only C++ compile-time assertions. The C API
uses one view representation and explicit runtime metadata, so it has no
meaningful equivalent type-concept or overload-identity assertion:

- `view_category_sfinae:ViewCategorySFINAE.VectorTypes`
- `view_category_sfinae:ViewCategorySFINAE.MatrixViews`
- `view_category_sfinae:ViewCategorySFINAE.SparseMatrixViews`
- `view_category_sfinae:ViewCategorySFINAE.MatrixTypesExcluded`
- `grouped_views:GroupedView.CategoryTagsAreDistinct`

Runtime views, access/mutability checks, grouped iteration, fill/copy and
reductions are exercised by the passing cases.
