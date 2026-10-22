# Validation results

## Tier 2: current machine and final defaults

Apple M5 Pro, 18 CPU cores; clang 21.0.0 and GCC 15.3.0. All builds use
`-O3 -DNDEBUG`, C11, pedantic diagnostics and the appropriate OpenMP flags.
The pinned conda environment provides build-time Python only.

| Check | clang | gcc-15 |
|---|---:|---:|
| `make test-cpu`: exact 680 reference IDs | 675 PASS, 0 FAIL, 5 C++-only N/A | Same |
| Reader validation | 10 PASS | 10 PASS |
| Generator structural/edit checks | 4 PASS | 4 PASS |
| `make test-specialized`: 70 numerical cases | PASS | PASS |
| TS1, 100 cells, seed 1, 100 x 30 s | PASS, bit-identical | PASS, max relative 3.27e-14 |
| TS1 specialized, 10,000 cells, seed 2026, 100 x 30 s | PASS, max relative 2.31e-14 | PASS, max relative 1.43e-13 |
| Generic TS1, both sizes and both timed thread counts | PASS | PASS |
| Five disclosed scenarios: final states and snapshots | All PASS | All PASS |
| Two additional seasonal/partial-group scenarios | All PASS | All PASS |

The 70 numerical cases cover seven generated solver variants, widths
1/3/32/37/64, and default/tight absolute tolerances. They compare calculated
rates, forcing, every packed Jacobian entry, LU factors, triangular solutions,
six external integration steps and solver statistics against the generic
library. Environmental and custom inputs vary during the runs, including darkness.
TS1 and all three Tier 1 configuration-file mechanisms are included. Extra
fixtures cover all configurable rate types, repeated reactants/products,
a parameterized reactant, an isolated species and zero-order chemistry.
An alternate variant combines optional fusion, row-oriented LU updates,
64-unit batches, seven-unit helpers, a preferred width of three and tiles of two.
The numerical bound is `1e-30 + 1e-10 * abs(reference)`.

`tests/cases.txt` is byte-identical to the supplied reference list. Assertions
remain active under `-DNDEBUG`. Switching directly from GCC back to clang,
without cleaning, also rebuilds and passes both suites. Linked-library inspection
of both specialized binaries shows only the system C runtime and OpenMP.
Source inspection finds no long double, intrinsics, compiler vector extensions,
assembly, target attributes, compiler builtins or operating-system-specific APIs.

## Benchmark and build measurements

10,000 cells, seed 2026, 100 timed steps of 30 seconds after one untimed step.
These are single local runs, not the experimenter's median-of-five measurements.
CSV parsing, construction and output are outside timing. All C runs use 32-cell
groups (313 groups); MICM uses the reference driver's default of one contiguous
group per thread. All final values pass the reference comparator; the full-size
comparison checks 2,110,000 numeric entries by column name.

| Implementation | 18 threads, elapsed_ms | 1 thread, elapsed_ms |
|---|---:|---:|
| Generated C, clang | 2,611.020 | 36,171.686 |
| Generic C, clang | 3,442.730 | 43,871.219 |
| MICM, clang++ | 11,812.435 | 138,888.809 |
| Generated C, gcc | 2,720.704 | 38,006.544 |
| Generic C, gcc | 3,840.910 | 52,159.734 |

The selected clang specialization is 1.32 times faster than the generic C
benchmark on 18 threads and 1.21 times faster on one thread; it is 4.52 and 3.84
times faster than these MICM measurements. C statistics are 281,700 accepted
internal steps, zero rejected steps, 563,400 forcing calls, 281,700 Jacobian
updates/factorizations and 845,100 triangular solves. MICM reports 16,200 internal
steps on 18 threads and 900 on one thread because its group count differs.

A clean `make -j6 CC=clang bench` takes **18.43 s**; the corresponding
`CC=gcc-15` build takes **42.61 s**. Both include Python generation and all three
benchmark/scenario binaries, with the conda environment already present; tests
and reference C++ compilation are excluded. Generated TS1 C is **3,431,254 bytes**.

Exact benchmark form:

```sh
.conda/bin/python ../reference/make_cells_csv.py 10000 2026 --output results/cells_10000.csv
bin/ts1_bench_specialized --threads 18 10000 100 30 results/cells_10000.csv
bin/ts1_bench_cpu --threads 18 10000 100 30 results/cells_10000.csv
# Repeat each with --threads 1; run in separate output directories.
```

SPECIALIZATION.md records individual kernel measurements and the attempted
specializations, including the slower alternatives. DESIGN.md lists every
parameter and the selected defaults.

## Day-night agreement

All runs use consecutive 32-cell groups, `--rtol 1e-6 --atol 1e-15`, no chemistry
warm-up, and the scenario's original cells and photolysis table. Timing includes
interpolation, photolysis scaling, rate updates and integration. Snapshot copies
and CSV output are outside timing.

Entries below give significant digits overall (worst species), scored from all
snapshots by the supplied tool. The required gate is 10 overall and 9 for every
species. Both snapshots and final states pass for both compilers, with zero
nonfinite values and zero concentrations below -1e-20. Against the high-accuracy
references the two compiler results round to the same three-decimal scores.

| Scenario | elapsed_ms, clang / GCC | Against MICM, clang | Against MICM, GCC | Against high-accuracy reference |
|---|---:|---:|---:|---:|
| diurnal-30min | 495.847 / 522.206 | 14.369 (13.810) | 14.366 (13.591) | 4.748 |
| diurnal-10min | 924.748 / 969.100 | 14.287 (13.649) | 14.296 (13.688) | 5.052 |
| diurnal-5min | 1505.977 / 1554.152 | 14.235 (13.631) | 14.243 (13.701) | 5.214 |
| polluted-10min | 1094.477 / 1136.453 | 14.076 (13.430) | 14.083 (13.424) | 5.398 |
| clean-10min | 875.468 / 915.089 | 14.350 (13.842) | 14.345 (13.845) | 4.907 |

The reference-solution agreement intentionally remains MICM's roughly five
significant digits: the method, error estimate and step control were preserved.
The internal step counts match the supplied MICM counts: 48,619; 90,980;
145,891; 106,393; and 84,273 in the table's order.

Two additional scenarios use independently generated MICM results at identical
tolerances and group sizes. `polar` uses 67 cells, seed 6147, latitude 66 N,
day 355, start 22 h, NOx x3, and 144 steps of 1,200 s (48 hours).
`southern` uses 35 cells, seed 961, latitude 35 S, day 172, start 17 h,
NOx x0.2, and 240 steps of 600 s (40 hours). Its snapshot interval of 36 steps
also exercises the final snapshot when the step count is not a multiple of the
interval. Both exercise partial final groups.

| Additional scenario | Against MICM, clang | Against MICM, GCC |
|---|---:|---:|
| polar | 14.495 (13.857) | 14.519 (13.973) |
| southern | 14.467 (14.050) | 14.385 (13.900) |

Reproduce the disclosed gates after building each compiler with
`tools/check_scenarios.sh LABEL 18`. The script checks final states and snapshots
with `score.py --level 10` and separately scores against the high-accuracy
reference. Development logs and CSV outputs are under ignored `results/`.

## Tier 1 results retained from the earlier machine

The following historical measurements describe the original 10-core machine
and GCC 15.2. Current-machine Tier 1 regression and benchmark results are above.

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

### Local benchmark measurements

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

### N/A cases

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
