# Tier 2 validation and performance

Measured on 2026-10-22 on the supplied Apple M5 Pro, 18 cores and 24 GiB,
macOS 26.7.1, with Apple Clang 21.0.0 and Homebrew GCC 15.3.0. MICM is the
read-only e15898ce input. Every implementation uses `-O3 -DNDEBUG`; C adds only
C11, warnings, includes, OpenMP and link flags. These are local measurements;
the experimenter's independent timing determines the ranking.

## Final correctness checks

Both compilers pass all 680 retained CPU cases and all 70 specialized cases.
The latter cover seven mechanisms, five group sizes and two tolerance settings,
including changing runtime inputs, copied states, padded tiles, every rate
family supported by the reader, duplicate terms/labels, isolated species, and
zero reactions. Component and full-solution comparisons enforce a pure relative infinity-norm
error of at most 1e-10, plus per-value comparisons at rtol=1e-10 and atol=1e-30
(the supplied TS1 gate's zero-scale rule). The branched fixture also verifies the reader's
existing conversion of its atom-count parameter to the library's integer type.

Both generic and specialized benchmarks pass all nine supplied 100-cell TS1
settings with each compiler: scalar, irregular/oversized groups, multiple thread
counts, static/dynamic/guided scheduling, and generic evaluation. Both binaries
also pass the full 10,000-cell seed-2026 gate with each compiler. The specialized
maximum relative differences are 2.96e-13 with Clang and 2.75e-13 with GCC,
well within the required 1e-10. No tolerances were relaxed.

All five disclosed day-night scenarios pass on both compilers at group size 32,
rtol=1e-6 and atol=1e-15. Both final states and every snapshot are checked with
the supplied scorer. All values are finite and none are below the negative
floor. The table reports snapshot digits overall, with the worst species in
parentheses. The last column is comparison to the supplied high-accuracy
reference; it reflects the retained method's discretization error, not the
much smaller difference from MICM's three-stage result.

| Scenario | Clang vs MICM | GCC vs MICM | Clang vs high-accuracy reference |
| --- | ---: | ---: | ---: |
| diurnal-30min | 14.377 (13.828) | 14.363 (13.791) | 4.748 (3.610) |
| diurnal-10min | 14.285 (13.656) | 14.292 (13.687) | 5.052 (3.824) |
| diurnal-5min | 14.230 (13.688) | 14.232 (13.537) | 5.214 (4.328) |
| polluted-10min | 14.047 (13.305) | 14.050 (13.432) | 5.398 (4.257) |
| clean-10min | 14.342 (13.793) | 14.335 (13.811) | 4.907 (3.766) |

The additional 48-hour scenario described in DECISIONS.md passes against a fresh
MICM run: snapshot agreement is 14.178 digits overall (13.712 worst species) for
Clang and 14.181 (13.755) for GCC. It exercises 450-second interpolation, a
different season/start time/NOx scale, and a one-cell tail group. An alternative
GCC build also passes all 70 specialized cases, the TS1 gate, and polluted-10min
with tile width 3, default group size 7, chemistry/algebra batches 3/5, function
size 97, and Mozart/division. The scenario explicitly selects 32-cell groups.

The case registry matches the supplied 680-case list byte for byte. Retranscribing
TS1 from its read-only header reproduces configs/ts1.json byte for byte. Source
inspection found no intrinsics, extended vector types, target attributes,
long double, offload, operating-system APIs, or forbidden compiler flags in the
C implementation. The Clang executable links only libomp and the system C
runtime. Portability was tested with the two supplied compilers on this host;
no Linux execution is claimed.

## Final benchmark timings

Input: 10,000 cells, seed 2026, 100 timed steps of 30 seconds after one untimed
warm-up. Runs were sequential, without overlapping compilation or benchmarks.
All timed-run final states were checked against the fresh MICM result. The
specialized 18-thread result is the median of five runs; the other 18-thread
results are medians of three. One-thread results each have one sample.

| Implementation | Threads | Cells per group | Elapsed ms |
| --- | ---: | ---: | ---: |
| Specialized C / Clang | 18 | 32 | 2,012.084 |
| Tier 1 C / Clang | 18 | 32 | 3,122.687 |
| MICM / clang++ | 18 | 555 or 556 | 11,878.097 |
| Specialized C / Clang | 1 | 32 | 27,458.710 |
| Tier 1 C / Clang | 1 | 32 | 38,235.266 |
| MICM / clang++ | 1 | 10,000 | 138,063.983 |

The specialized solver is 1.55 times faster than Tier 1 and 5.90 times faster
than the supplied MICM benchmark at 18 threads. The respective one-thread ratios
are 1.39 and 5.03. MICM uses its supplied default of one group per thread; the C
benchmarks use 32-cell groups. Each cell takes nine internal steps per requested
interval at the default tolerances. Both C implementations report 281,700
internal/accepted steps, zero counted rejections, 563,400 forcing calls,
281,700 Jacobians/factorizations and 845,100 solves over the timed region. MICM
reports the corresponding statistics over its smaller number of groups.

## Build and reproduction

After `make clean`, with the local conda environment already present,
`make -j6 CC=COMPILER bench` took 46.53 s with Clang and 38.36 s with GCC.
The default TS1 generated translation unit is 2,607,627 bytes and 75,172 lines.
It includes descriptor/inspection kernels as well as the packed solver. All
C build products are regenerated from the configuration and compiled through
Make with the selected compiler and flags.

```sh
make env
make bench test-cpu test-specialized CC=clang
.conda/bin/python tools/check_ts1.py --output .runs/gate-cpu
.conda/bin/python tools/check_ts1.py --binary bin/ts1_bench_specialized --output .runs/gate-specialized
.conda/bin/python tools/check_scenarios.py --output .runs/day-night
make bench test-cpu test-specialized CC=gcc-15
```

SPECIALIZATION.md records the performance experiments, including the slower
alternatives. Compact measured results, per-scenario scores, build metadata and
the generated-source checksum are in results/tier2/. Full local CSVs and logs
are retained under ignored .runs/. The following section preserves the Tier 1
measurements on its original machine; it is historical, not a second set of
Tier 2 timings.

# Tier 1 retained results

Measurements were made on the supplied Apple M5 host (10 CPU cores, 16 GB),
macOS 27.0.1, using Apple Clang 21.0.0 and GCC 15.2.0. Every C object used
`-O3 -DNDEBUG -std=c11` and the documented warning/include/OpenMP flags.
MICM was built with `make micm-ref`, using clang++ and the same `-O3 -DNDEBUG`.
No additional optimization, architecture, floating-point, LTO or profiling flags
were used.

## CPU cases

| Compiler | Cases passing | Failing | N/A |
| --- | ---: | ---: | ---: |
| Apple Clang 21 | 680 | 0 | 0 |
| GCC 15.2 | 680 | 0 | 0 |

The registry contains exactly the 680 supplied IDs, in the original output
order. Tests include all original solver parameter sets and matrix-layout
variants, analytical trajectories, dense/sparse operations, rate families,
constraint initialization, external models, finite differences, and concurrent
states. Both flow-tube cases compare all 3,601 supplied rows at atol=1e-10
(the original tests checked ten rows). Chapman regressions exercise both
hand-built and automatically generated kernels against an independent legacy
scalar oracle.

All eight tutorial programs run through `bin/tutorial` and the test harness.
Every deterministic output byte matches the captured MICM output, including
printed concentrations and solver statistics. Only measured elapsed time and
thread announcement order are handled as nondeterministic fields; DECISIONS.md
describes their validation. Each of the three YAML/JSON configuration examples
also matches an independently hand-built mechanism: rate constants, forcing,
Jacobians and ten multi-cell integration intervals, using both LU algorithms
and both generic and generated evaluators.

Both full builds are warning-free. Runtime dependency inspection of the Clang
benchmark shows only libomp and the system C runtime. Source/build inspection
found no nonstandard C APIs, long double, intrinsics, target-specific constructs,
C++/Python runtime dependencies or forbidden compiler flags. Linux execution
was not available; portability is supported by the ISO C11/OpenMP implementation
and compiler checks on this host, rather than a claimed Linux test.

## TS1 numerical gate

All comparisons use the supplied `compare_final.py`, with the unchanged gate
`abs(y - reference) <= 1e-30 + 1e-10 * abs(reference)`. Each run includes one
warm-up and 100 timed steps of 30 seconds. The comparison includes every cell,
209 state species, temperature and pressure.

- Supplied 100-cell seed 1: Clang and GCC pass default execution, threads/groups
  (1,1), (4,7), (1,32), (3,64), (4,101), and generic evaluation (2,13).
  Largest relative differences: 4.25e-14 (Clang), 5.30e-14 (GCC). The final
  Clang rebuild also passes the dynamic and guided variants on this input.
- Fresh 10,000-cell seed 2026 reference: all five Clang all-core runs pass all
  2,110,000 values. Largest relative difference: 2.03e-13.
- Fresh MICM ten-thread and one-thread results agree exactly for seed 2026.

- The Clang one-thread 10,000-cell run passes with the same 2.03e-13 maximum.
- Fresh 137-cell seed 314159: Clang passes all nine settings, including groups
  of 1/7/13/31/32/64/138, irregular tails, a group larger than the domain, generic
  evaluation, and static/dynamic/guided scheduling. Maximum relative difference
  is 7.58e-14. GCC passes the same nine settings with a maximum of 4.46e-14.
- GCC also passes all 2,110,000 values in the 10,000-cell run, with a maximum
  relative difference of 2.14e-13.

TS1 translation is reproducible: running `tools/translate_ts1.py` on the supplied
header reproduces `configs/ts1.json` byte for byte, with 210 declared species
(209 variables and one air-density parameter) and 547 reactions.

## Timings

The input is generated by the supplied `make_cells_csv.py 10000 2026`.
The measured interval excludes configuration/model construction, CSV reading,
rate initialization, the warm-up and final CSV writing. Runs were sequential;
no compilation or other benchmark was launched concurrently.

| Implementation | Threads | Cell group | Runs | Elapsed ms |
| --- | ---: | ---: | ---: | ---: |
| MICM / clang++ | 1 | 10,000 | 1 | 143,929.029 |
| MICM / clang++ | 10 | 1,000 | 1 | 28,084.099 |
| C / Clang | 1 | 32 | 1 | 38,283.845 |
| C / Clang | 10 (default) | 32 | 5 | **7,306.806 median** |
| C / GCC 15 | 10 (default) | 32 | 1 | 7,785.417 |

The five C times were 6,893.911, 7,146.751, 7,306.806, 7,520.881 and 7,473.683 ms.
The median is 19.7 times faster than this one-thread MICM measurement, and 3.84
times faster than the ten-thread MICM measurement. The one-thread C run is
3.76 times faster than the one-thread MICM run. These are local measurements,
not a prediction of the evaluator's idle-machine median or another architecture.
The reference comparisons have one timed sample each; the supplied reference
README separately reports a one-thread median of 137,383 ms for this seed.

Default C timed statistics: 281,700 internal/accepted steps, zero rejected steps,
563,400 forcing calls, 281,700 Jacobian updates and factorizations, and 845,100
linear solves. MICM's one-thread statistics are 900 accepted/internal steps,
1,800 forcing calls, 900 Jacobians/factorizations and 2,700 solves. The difference
is the number of independently controlled groups (313 versus one); each group
uses nine accepted internal steps per requested interval. No tolerance or solver
formula was changed to obtain the timing.

A 20-step pilot on the same 10,000 cells compared group sizes 16/32/64/128/256,
with static and dynamic scheduling. Static times were 1,419/1,431/1,511/1,716/2,605
ms; dynamic times were 1,404/1,457/1,494/1,703/2,581 ms. The small difference
between 16 and 32 and between schedulers did not justify replacing the existing
32-cell/static default. All choices remain runtime parameters.

## Reproduction

```sh
make env
make test-cpu CC=clang
make test-cpu CC=gcc-15
make cpu CC=clang
make micm-ref
.conda/bin/python tools/check_ts1.py --output build/ts1-checks
```

`tools/check_ts1.py` accepts another input/reference pair, cell count, and binary.
It checks scalar, irregular-tail, oversized, static/dynamic/guided and generic
variants. `tools/benchmark.py` runs sequential repetitions, saves stdout and final
CSV files, optionally applies the supplied gate after each run, and records all
measurements and their median in `summary.json`. For example:

```sh
.conda/bin/python tools/benchmark.py --input CELLS_CSV --output build/timings \
  --expected MICM_FINAL_CSV --repeats 5
```

Compact machine-readable evidence is retained in `results/benchmarks.json` and
`results/validation.json`. Full local run logs and CSVs are under ignored `build/`.

Keep outputs inside this checkout. `make_cells_csv.py` writes a named CSV in its
working directory; run it from a directory under `build/` when regenerating inputs.
