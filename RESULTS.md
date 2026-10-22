# Tier 2 validation results

Validated on the supplied 18-core Apple M5 Pro with GNU Make 3.81, Apple Clang
21.0.0 and GCC 15.3.0. Both C compilers use `OPTFLAGS=-O3 -DNDEBUG`; the Makefile
adds only C11, warning, include, OpenMP and link flags. All results below use the
retained defaults: 32-cell groups/stride, 32 operations per helper, fusion of one
operation, exact rate caching and the fused three-stage estimate.

## Tests and build products

| Check | Clang | GCC 15.3 |
| --- | --- | --- |
| `make cpu bench` from clean C build | PASS | PASS |
| `make test-cpu` reference cases | 680 PASS, 0 FAIL, 0 N/A | 680 PASS, 0 FAIL, 0 N/A |
| Tier 1 supplemental configuration/tutorial/TS1 checks | PASS | PASS |
| `make test-specialized`, six mechanisms | PASS | PASS |
| Full 10,000-cell final-state gate, 18 and 1 threads | PASS | PASS |
| Five disclosed scenarios, snapshots and final states | PASS | PASS |
| Additional 48-hour seasonal/interpolation case | PASS | PASS |

The six mechanisms are TS1, the two rate-constant tutorials, Robertson, the
additional-rate fixture and the regenerated mutation fixture. There are 60
integration trajectories (six mechanisms, five cell counts, two tolerances),
each checked over 12 intervals. Rates and every final concentration use a strict
`1e-10` relative bound, with no absolute allowance; expected zero requires zero.
Independent factor/solve tests, singular/NaN pivot detection, runtime rate-input
changes, cache hits/invalidation and solver statistics are also checked.
Clang additionally passed variable-stride generation with helper size 17 and
fusion 3 at all tested cell counts. No warnings occurred in either default build.

Clean-build wall time excludes conda provisioning and the optional MICM C++
reference executables. It measures `make -j6 cpu bench CC=...` after `make clean`.
The pinned Python 3.12.12/PyYAML 6.0.3 environment was already provisioned.

| Compiler | Clean build | Generated TS1 C |
| --- | ---: | ---: |
| clang | 20.70 s | 2,048,850 bytes |
| gcc-15 | 26.06 s | 2,048,850 bytes |

Generated C includes the builder, rate expressions, chemistry, sparse algebra
and template instantiation. It is regenerated, not committed or manually edited.
Both source and generated-C audits found no forbidden compiler extensions,
intrinsics, target pragmas, extended precision or OS-specific runtime APIs.
Linux execution was not available; portability was checked from C11 source and
both compiler builds on the supplied macOS machine.

## Full TS1 benchmark

10,000 cells, seed 2026, one untimed warm-up and 100 timed steps of 30 seconds,
MICM default tolerances. Times are single sequential diagnostic runs, not the
experimenter's idle-machine medians. MICM is the supplied reference driver built
with Clang++ and the same optimization flags in both sets of measurements.
It uses its default one group per thread; both C drivers use groups of 32.

| C compiler | Threads | Specialized ms | Tier 1 ms | MICM ms |
| --- | ---: | ---: | ---: | ---: |
| clang | 18 | 3,191.89 | 4,000.28 | 11,270.71 |
| clang | 1 | 41,728.77 | 48,800.11 | 136,514.56 |
| gcc-15 | 18 | 3,086.21 | 4,586.56 | 11,378.00 |
| gcc-15 | 1 | 39,370.98 | 59,143.42 | 137,491.18 |

At 18 threads, the Clang specialized solver is 1.25 times as fast as Tier 1
and 3.53 times as fast as MICM in these runs. Its one-thread factors are 1.17
and 3.27 respectively. SPECIALIZATION.md records the development experiments,
including slower alternatives.

Every comparison checked 2,110,000 values through the supplied `compare_final.py`
with relative tolerance `1e-10` and absolute tolerance `1e-30`. At both thread
counts the maximum relative differences from MICM were:

| Compiler | Generic C | Specialized C |
| --- | ---: | ---: |
| Clang | 0 | 2.31e-14 |
| GCC 15.3 | 1.43e-13 | 1.43e-13 |

Both C drivers report 281,700 internal/accepted steps, zero rejected steps,
563,400 function calls, 281,700 Jacobian updates/decompositions and 845,100
linear solves. MICM reports 16,200 internal steps at 18 threads and 900 at one
thread because its groups differ; statistics count group steps, not cell steps.

The supplied 100-cell seed-1 result also passes for both drivers and compilers,
with four threads, groups of 7 (`dynamic,1`) and 32 (`static`). Maximum relative
difference is zero for Clang and 3.27e-14 for GCC, using the same gate thresholds.

## Day-night results

All runs use 32-cell groups, `--rtol 1e-6`, `--atol 1e-15`, no warm-up and the
prescribed photolysis updates. Every disclosed case has 1,000 cells and four
snapshots (4,000 rows). The supplied `score.py --level 10` passes for both
snapshots and final states on both compilers. All outputs have zero nonfinite
values and zero values below `-1e-20` mol m-3.

Entries below are significant digits, overall / worst species, against MICM's
provided three-stage results. Reference agreement uses the provided high-accuracy
six-stage solutions; both compilers agree on these reference scores to the
displayed precision. Times include photolysis changes, rates and integration,
and exclude input/setup and snapshot output.

| Scenario | Clang ms | GCC ms | Clang vs MICM | GCC vs MICM | vs high-accuracy reference |
| --- | ---: | ---: | --- | --- | --- |
| diurnal-30min | 627.32 | 574.16 | 14.369 / 13.810 | 14.366 / 13.591 | 4.748 / 3.610 |
| diurnal-10min | 1,206.37 | 1,105.17 | 14.287 / 13.649 | 14.296 / 13.688 | 5.052 / 3.824 |
| diurnal-5min | 1,934.79 | 1,781.87 | 14.235 / 13.631 | 14.243 / 13.701 | 5.214 / 4.328 |
| polluted-10min | 1,382.54 | 1,284.93 | 14.076 / 13.430 | 14.083 / 13.424 | 5.398 / 4.257 |
| clean-10min | 1,109.75 | 1,041.75 | 14.350 / 13.842 | 14.344 / 13.844 | 4.907 / 3.766 |

The final-state-only comparisons also pass: the minimum agreement over both
compilers and all disclosed cases is 14.054 digits overall and 13.246 for any
species. Final-state agreement with the high-accuracy solutions ranges from
4.840 to 6.077 digits overall. These lower reference scores reflect the retained
MICM numerical method and requested tolerances, not a specialization change.

The additional case runs 65 cells for 48 hours with 450-second steps, latitude
65 N, day 172, a 21:00 start, NOx scaled by 7 and seed 337. It exercises
photolysis interpolation halfway between tabulated times, a one-cell tail group
and eight snapshots. Against a fresh MICM run at identical settings:

| Compiler | Snapshot overall / worst species | Final overall / worst species |
| --- | --- | --- |
| clang | 14.193 / 13.626 | 14.160 / 13.652 |
| gcc-15 | 14.175 / 13.628 | 14.123 / 13.572 |

All three implementations take exactly 26,546 internal/accepted steps in this
additional case, with no rejected steps. The new MICM executable is built only
from the supplied read-only reference source using `make micm-scenario-ref`.

## Evidence and reproduction

`tools/validate_specialized.py --compiler clang` (or `gcc-15` after rebuilding)
reproduces the full benchmark and all disclosed scenario scores. Build logs,
full CSVs and JSON summaries are in `validation/<compiler>/`; the extra case is
in `validation/seasonal/`. Development comparisons are in `validation/tuning/`.
README.md documents generation, testing and runtime controls.

The sections below preserve Tier 1 measurements from its earlier machine and
compiler version. They are historical records, not the Tier 2 timings above.

# Historical Tier 1 results

All 680 reference cases and supplemental checks pass on both compilers. The
full TS1 gate and all additional seed/group/schedule comparisons pass.

## Case results

Both compiler builds use `OPTFLAGS='-O3 -DNDEBUG'` and the Makefile's corresponding
OpenMP flags. The ported cases keep original or tighter numerical bounds.

| Compiler | PASS | FAIL | N/A |
| --- | ---: | ---: | ---: |
| Apple Clang 21 | 680 | 0 | 0 |
| GCC 15.2 | 680 | 0 | 0 |

The 680 cases comprise 31 rate-constant tests, 11 rate-store cases, 22 state
cases, 32 LU decomposition
cases, 20 linear-solver cases, 2 backward Euler cases, 42 constraint
definition/set cases, 12 constraint initialization cases, 12 DAE integration
cases, 8 Rosenbrock utility/regression cases, 10 process cases and 7 Jacobian
verification cases, all 117 dense/vector matrix cases, 14 species/phase/system
cases, 5 version cases, callback-rate override, OpenMP state-isolation, and all
190 sparse matrix cases, 27 grouped-view/view-category cases, 20 external-model
constraint cases, 18 aerosol integration cases and two Terminator integration
cases, 13 Chapman unit cases, one Chapman integration case, 12 Chapman
regression cases, and two flow-tube regression cases, plus all 32 analytical
integration cases, eight solver-builder cases and
eight complete tutorial programs. Tutorial comparisons preserve every printed
result and counter; measured durations and thread announcement interleaving are
validated and normalized as documented in DECISIONS.md. Additional
checks pass on both compilers: all three tutorial configurations versus
independent hand-built mechanisms (one and seven cells, ten integration
intervals), 49 configuration validation checks, independent formulas for three
additional generated rate types, all 379 TS1 initial fields by name, and a TS1
interval-clamp regression. Builds have no compiler warnings.

## TS1 final-state gate

MICM was rebuilt with `make micm-ref` using Clang++ and the same optimization
flags. All comparisons below use the required absolute tolerance `1e-30` and
relative tolerance `1e-10`, through the supplied `compare_final.py`.

| Input and run | C configuration | Checked values | Largest relative difference | Result |
| --- | --- | ---: | ---: | --- |
| 100 cells, seed 1, 100 steps, 30 s | Clang, 1 thread, groups of 32 | 21,100 | 9.7e-13 | PASS |
| Same | GCC, 1 thread, groups of 32 | 21,100 | 9.7e-13 | PASS |
| Same | GCC, 4 threads, groups of 7 | 21,100 | 9.7e-13 | PASS |
| 10,000 cells, seed 2026, 100 steps, 30 s | Clang, 10 threads, groups of 32 | 2,110,000 | 0 | PASS |
| Same | GCC, 10 threads, groups of 32 | 2,110,000 | 1.43e-13 | PASS |

Final extended runs also checked the following settings on both compilers,
all for 100 steps of 30 seconds after warm-up:

| Cells / seed | Group / threads / OpenMP schedule | Clang max relative difference | GCC max relative difference |
| --- | --- | ---: | ---: |
| 127 / 17 | 1 / 1 / static | 0 | 4e-14 |
| 127 / 17 | 7 / 4 / dynamic,1 | 0 | 4e-14 |
| 127 / 17 | 32 / 10 / guided,2 | 0 | 4e-14 |
| 127 / 17 | 1000 / 3 / static,2 | 0 | 4e-14 |
| 129 / 73 | 3 / 10 / static,1 | 0 | 3.92e-14 |
| 129 / 73 | 64 / 4 / dynamic,3 | 0 | 3.92e-14 |
| 129 / 73 | 129 / 1 / static | 0 | 3.92e-14 |

A group request larger than the cell count is capped to that count. The 127-cell
runs compare 26,797 fields each; the 129-cell runs compare 27,219 fields each.

The 10,000-cell MICM reference used ten threads, one group per thread. The C
benchmark used 313 groups, each with independent adaptive solver state. All
runs include an untimed warm-up interval before the timed intervals.

The retained one-cell TS1 regression checks the public solver's post-interval
nonnegative concentration clamp. It catches a difference that the 100-cell
sample alone did not expose.

## Diagnostic timings

Single runs on the supplied laptop; these are development measurements, not
idle-machine medians for ranking. No compiler optimization flags were tuned.

| Run | Clang C | GCC C | MICM Clang++ |
| --- | ---: | ---: | ---: |
| 100 cells, seed 1, 100 steps, 30 s, one thread | 516.17 ms | — | 1,316.91 ms |
| 10,000 cells, seed 2026, 100 steps, 30 s, ten threads | 9,625.25 ms | 10,602.78 ms | 27,421.51 ms |

The final GCC validation separately measured MICM at 26,035.22 ms. Clang C was
2.85 times faster than its paired MICM run. These measurements are single runs;
the evaluation harness determines the median on an idle machine.

For the full-size C run: 281,700 internal/accepted steps, zero rejected steps,
563,400 function calls, 281,700 Jacobian updates and decompositions, and 845,100
linear solves. MICM's larger groups report 9,000 internal/accepted steps,
18,000 function calls, 9,000 Jacobian updates and decompositions, and 27,000
solves. These counts sum over groups, not cells.

The final source audit found no forbidden extensions, target flags, extended
floating types or operating-system APIs in the C implementation. The case list
matches the supplied inventory exactly, and the TS1 translation reproduces the
committed mechanism. All C and generated C use the Makefile compiler and flags.
Linux execution was not available in this environment; the portability audit
and both compiler runs were performed on the supplied macOS machine.

Reproduce the extended numerical comparisons with `tools/validate_benchmarks.py`
as described in README.md. Logs and full comparison files are retained under
`build/validation/final-clang/` and `build/validation/final-gcc/`.
