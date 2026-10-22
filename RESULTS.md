# Tier 2 validation and measurements

Measured on 2026-10-22 on the supplied 18-core Apple M5 Pro, macOS 26.7.1, with Apple clang 21.0.0 and GCC 15.3.0. All C and the independent clang++ MICM references use `OPTFLAGS='-O3 -DNDEBUG'`. Added flags provide C11, warnings, includes, shared-library support and OpenMP. There are no added optimization, target, floating-point, LTO or profile-guided flags.

## Correctness gates

| Compiler | Tier 1 PASS | N/A | FAIL | Strict generator comparisons | Disclosed scenarios |
| --- | ---: | ---: | ---: | ---: | ---: |
| clang | 675 | 5 | 0 | 42 PASS | 5 PASS |
| gcc-15 | 675 | 5 | 0 | 42 PASS | 5 PASS |

Each CPU run reports all 680 identifiers in the supplied case list, unchanged and in order. The five N/A cases remain the four compile-only `ViewCategorySFINAE` cases and `GroupedView.CategoryTagsAreDistinct`; C has no corresponding template/category machinery. All eight C tutorial programs retain their reference output checks. Generic rates, arrays, all five Rosenbrock tables, backward Euler, sparse solvers, constraints, external models and the other Tier 1 functionality remain covered.

`make test-specialized` generates and compiles TS1, all three configuration examples, the all-rate-family validation mechanism, the source/sink/duplicate-term fixture, and a TS1 variant with a changed kinetic coefficient and an added reaction. Each runs with 1, 2, 4, 7, 32 and 33 cells. Tests compare rates, forcing, Jacobians, LU, triangular solves and three integrated intervals with the independent generic library at strictly `1e-10 * abs(reference)`, requiring exact expected zeros. Temperature, pressure, independently supplied air density and all custom parameters change between intervals. The four-cell backend, every partial tile size, allocated scratch and scratch reuse are exercised. Both compiler suites exit zero.

The supplied 100-cell seed-1 MICM gate passes for both native benchmarks with each compiler: largest relative differences are `3.48e-14` (clang) and `3.43e-14` (GCC). A run with 17-cell adaptive groups and three workers also passes, including a partial final group. The clang specialized executable also passes with an empty environment apart from PATH and an OpenMP scheduling setting, without Python.

Every full-size benchmark below passes the supplied comparison against MICM over 2,110,000 values (10,000 cells × 211 columns). Largest relative differences are `1.26e-13` for both clang C programs and `1.16e-13` for both GCC C programs, within `1e-30 + 1e-10 * abs(MICM)`. The one-worker and 18-worker MICM outputs agree exactly. Static and dynamic scheduling produce byte-identical specialized output in the scheduling check.

Input: 10,000 cells, seed 2026, SHA-256 `43febf31a0a1f8e0acfe02d93ac96d95ecc024df9572547dd1000926a85ad8b5`. The independently computed MICM output SHA-256 is `133290093c64317d520dff78e89622db36dc467110eb0f3d0e4bc6b24f16ba7d`.

## Full-size benchmark

One untimed warm-up, then 100 timed steps of 30 s, MICM's default tolerances, 32-cell adaptive groups for every implementation, and static group scheduling. The native programs report the times, excluding initialization, warm-up and CSV output. Multiple samples use their median; single samples are individual observations. These local measurements are separate from the experimenter's ranking runs.

| Program | Workers | Samples | Elapsed ms | Speedup over MICM |
| --- | ---: | ---: | ---: | ---: |
| MICM, clang++ | 1 | 1 | 126,021.018 | 1.00× |
| Tier 1, clang | 1 | 1 | 36,508.885 | 3.45× |
| Specialized, clang | 1 | 1 | 31,376.901 | 4.02× |
| Tier 1, gcc-15 | 1 | 1 | 39,647.217 | 3.18× |
| Specialized, gcc-15 | 1 | 1 | 28,456.502 | 4.43× |
| MICM, clang++ | 18 | 3 | 9,989.994 | 1.00× |
| Tier 1, clang | 18 | 5 | 2,927.945 | 3.41× |
| Specialized, clang | 18 | 5 | 2,330.483 | 4.29× |
| Tier 1, gcc-15 | 18 | 1 | 3,080.516 | 3.24× |
| Specialized, gcc-15 | 18 | 1 | 2,182.353 | 4.58× |

The primary clang specialized median is 1.26 times faster than Tier 1 (20.4% less time), and 4.29 times faster than MICM. Its five times were 2353.189, 2296.965, 2314.236, 2330.483 and 2348.061 ms. The generic clang times were 2853.358, 2842.138, 2927.945, 2974.044 and 3047.414 ms; MICM's 18-worker times were 9506.383, 9989.994 and 10676.704 ms. Each timed C run passed the full state gate.

With matching groups, all implementations report 281,700 internal/accepted steps, zero rejections, 563,400 forcing calls, 281,700 Jacobian updates/decompositions and 845,100 solves. Statistics count logical adaptive-group operations; execution tiles leave these totals unchanged. The different timings therefore come from execution/layout specialization rather than altered step counts or chemistry.

## Clean builds and generated size

| Compiler | Clean `make -j8 bench` wall time |
| --- | ---: |
| clang | 22.60 s |
| gcc-15 | 41.75 s |

These measurements start after `make clean`, with the private conda environment already installed. They include configuration generation, the generic static/shared library, tutorial runner, generic benchmark, specialized benchmark and scenario runner; they exclude creating the environment and building extra test fixtures. Both clean builds have no warnings.

The complete generated TS1 C file is 2,222,161 bytes and 66,075 lines, identical for both compilers. Generated files are ignored build products. The final clang benchmark and scenario executables are byte-identical to the executables used for the recorded timings and scenario scores. The generator emits recursive Make dependencies; touching an included reaction file was verified to regenerate both generic and specialized C. See `SPECIALIZATION.md` for measurements of retained and discarded variants.

## Day-night results

Every disclosed run uses 32-cell adaptive groups, `--rtol 1e-6 --atol 1e-15`, 18 workers, no warm-up, and the prescribed photolysis interpolation/update. Both final files and every snapshot pass the supplied scorer at level 10 with MICM as reference. Across both compilers there are zero nonfinite values and zero values below `-1e-20` in either output format. The following table gives snapshot significant digits: overall, then the worst species in parentheses.

| Scenario | clang versus MICM | GCC versus MICM | Versus high-accuracy reference (both compilers) |
| --- | --- | --- | --- |
| diurnal-30min | 14.379 (ISOPNO3 13.775) | 14.368 (ISOP 13.796) | 4.748 (ISOP 3.610) |
| diurnal-10min | 14.285 (ISOP 13.710) | 14.286 (ISOP 13.649) | 5.052 (ISOP 3.824) |
| diurnal-5min | 14.242 (ISOPNO3 13.574) | 14.236 (ISOP 13.689) | 5.214 (ISOP 4.328) |
| polluted-10min | 14.065 (ISOP 13.435) | 14.059 (ISOP 13.460) | 5.398 (ISOP 4.257) |
| clean-10min | 14.343 (ISOP 13.809) | 14.359 (ISOP 13.851) | 4.907 (ISOP 3.766) |

For the final-state files alone, the lowest overall agreement with MICM is 14.041 digits and the lowest species agreement is 13.301 digits. Agreement with the high-accuracy references is the expected three-stage-method accuracy, about 4.75–5.40 digits overall, matching MICM's own scale rather than changing its numerical method.

Additional generated scenarios also pass final and snapshot level-10 gates with both compilers. The table gives overall digits, with worst-species digits in parentheses.

| Extra case | clang versus MICM | GCC versus MICM |
| --- | --- | --- |
| 65 cells; 48 h; latitude 55, day 355; NOx ×0.3 | 14.403 (13.871) | 14.403 (13.851) |
| 33 cells; 124 × 700 s; latitude -35, day 355; NOx ×5 | 14.340 (13.887) | 14.310 (13.710) |
| 9 cells; omitted photo column; irregular interpolation/endpoints; snapshots 3, 6, 7 | 14.945 (14.517) | 14.876 (14.342) |

The winter case uses seed 817, 600 s steps from 21 h local time, snapshots every 36 steps. The summer case uses seed 991 from 13 h, snapshots every 30 steps including final step 124. The small input case has photolysis rows `(100,0.25)`, `(260,1)`, `(500,0)`, 120 s steps, and one photolysis column omitted so its initialized default is used. These exercise changed latitude/season/start time/NOx, a 48-hour integration, interpolation between table rows and outside the endpoints, relative paths/comments, partial groups/tiles and a final snapshot outside the regular interval.

## Reproduction

```sh
make env
mkdir -p results
make clean
/usr/bin/time -p make -j8 CC=clang bench
make -j8 CC=clang test-cpu test-specialized > results/clang-tests.log
make micm-ref micm-scenario
.conda/bin/python ../reference/make_cells_csv.py 10000 2026 --output results/cells.csv
.conda/bin/python tools/measure.py bin/micm_ts1_ref results/cells.csv \
  --label micm-18 --threads 18 --repeats 3 --micm
.conda/bin/python tools/measure.py bin/ts1_bench_cpu results/cells.csv \
  --label tier1-clang-18 --threads 18 --repeats 5 \
  --reference results/micm-18/run-1/ts1_final.csv
.conda/bin/python tools/measure.py bin/ts1_bench_specialized results/cells.csv \
  --label specialized-clang-18 --threads 18 --repeats 5 \
  --reference results/micm-18/run-1/ts1_final.csv
.conda/bin/python tools/validate_scenarios.py --label clang-scenarios
make clean
/usr/bin/time -p make -j8 CC=gcc-15 bench
make -j8 CC=gcc-15 test-cpu test-specialized > results/gcc-tests.log
.conda/bin/python tools/validate_scenarios.py --label gcc-scenarios
```

Use `--threads 1 --repeats 1` with the measurement tool for single-worker runs. After the GCC build, run the same two C measurement commands with new GCC labels and the retained MICM CSV reference. The tool saves native stdout, final state, complete comparison output, commands, input/reference hashes and a JSON timing summary under each label. The scenario tool saves final/snapshot scores against both MICM and high-accuracy references. Large outputs, environments, caches and logs remain outside git.

The portability audit of runtime and generated C found only standard C headers and `omp.h`, CPU OpenMP constructs, double precision, and ordinary array arithmetic. No intrinsics, inline assembly, compiler vector types, target attributes/builtins/pragmas, wider floating-point type or nonstandard runtime API is used. The clang native executable links only the system C runtime and OpenMP, with no Python, MICM or MUSICA runtime dependency.
