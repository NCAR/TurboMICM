# Tier 2 verification and measurements

Apple M5 Pro, 18 CPU cores, macOS 26.7.1; Apple Clang 21 and GCC 15.3.0.
All C uses ISO C11, OpenMP, and `-O3 -DNDEBUG`, with no target or floating-point
tuning flags. Python 3.12.12 and PyYAML 6.0.3 run only during generation and tests.

## Gates

Both compiler builds pass:

- `make cpu`, `make bench`, and all 675 applicable `make test-cpu` cases;
  five C++ compile-time-only cases retain individual N/A explanations.
- `make test-specialized`: 36 generated-versus-generic comparisons across seven
  distinct configurations and two additional generator variants, at 1, 7, 32,
  and 33 cells. Concentrations agree to `1e-10` relative with no absolute floor.
- Tests cover all supported rate families, arbitrary runtime conditions, every
  custom parameter with distinct nonzero values, repeated reactants, source/loss
  terms, self-products, negative fixed rate coefficients, parameterized species,
  and an inert mechanism. They compare every sparse coordinate including fill and structural zeros, default
  and tight-tolerance trajectories, solver statistics, and all five method tables.
- Both TS1 benchmarks pass the 100-cell and 10,000-cell MICM final-state gates.
- All five disclosed day-night scenarios pass level 10, and two additional
  48-hour cases pass against independently run MICM references.
- All snapshot steps and cell indices are present; the final snapshots match
  the final-state CSV values exactly. A dynamic OpenMP schedule gives an
  identical benchmark final state.

| Largest relative TS1 difference from MICM | Clang | GCC |
|---|---:|---:|
| Generic, 100 cells | 4.25e-14 | 5.3e-14 |
| Generic, 10,000 cells | 2.03e-13 | 2.14e-13 |
| Specialized, 100 cells | 3.23e-14 | 3.66e-14 |
| Specialized, 10,000 cells | 2.96e-13 | 2.75e-13 |

The original generic API and its eight tutorial programs remain available.
The Tier 1 test coverage includes chemistry, all Rosenbrock methods, backward
Euler, constraints and initialization, arrays and sparse views, both LU algorithms,
external models, state ownership, metadata, and analytical/regression fixtures.

## Full benchmark

Input: 10,000 cells from seed 2026, 100 timed steps of 30 seconds after one
untimed warm-up; MICM default tolerances. Values below are medians of three
sequential local runs, measured by each executable's `elapsed_ms`. They are
observations on this machine, not the experimenter's ranking measurements.

| Program | 18 threads (ms) | 1 thread (ms) |
|---|---:|---:|
| Specialized C | 3014.059 | 35049.822 |
| Generic Tier 1 C | 4333.436 | 46475.159 |
| MICM | 12936.619 | 136001.636 |

The specialized solver is 1.44× faster than Tier 1 and
4.29× faster than MICM on 18 threads; the one-thread ratios are
1.33× and 3.88×.

Both C benchmarks use consecutive 32-cell groups, with a short final group.
MICM uses its benchmark default: one group per thread. All programs use the
same cells, step lengths, method, tolerances, and timing convention. Group-level
statistics differ accordingly: both C programs take 281,700 accepted group
steps and no rejected steps; MICM takes 16,200 on 18 threads and 900 on one.
Both C final states also pass comparison with MICM at each thread count.

## Clean builds and generated size

`make clean` followed by `make -j4 bench`, with the prefix environment already
created, measures generation and all generic/specialized benchmark builds.

| Compiler | Seconds | Generated C bytes |
|---|---:|---:|
| Clang | 13.77 | 990,603 |
| GCC 15 | 14.89 | 990,603 |

The generated file is identical under both compilers:

```text
fce4dbfb32c38705550a22ab2ab4357f6fa8840cb6da9ed9cd75235a5088a2ee
```

Default builds and both test targets complete without compiler warnings.
SPECIALIZATION.md records the measured emission alternatives, including slower variants.

## Day-night scenarios

Every run uses 32-cell groups, `--rtol 1e-6`, and `--atol 1e-15`, without a
warm-up. Timings include interpolation, photolysis updates, rate evaluation and
solves; snapshot copies and CSV output are excluded. Scenario times are single
observations. Digits are the supplied scorer's overall agreement and the lowest
per-species agreement over the snapshots.

| Scenario | Clang ms | GCC ms | MICM digits: Clang overall / worst | MICM digits: GCC overall / worst | High-accuracy reference: overall / worst |
|---|---:|---:|---:|---:|---:|
| diurnal-30min | 529.312 | 549.390 | 14.376 / 13.828 | 14.363 / 13.791 | 4.748 / 3.610 |
| diurnal-10min | 995.987 | 1007.047 | 14.285 / 13.655 | 14.292 / 13.686 | 5.052 / 3.824 |
| diurnal-5min | 1615.883 | 1607.240 | 14.230 / 13.688 | 14.233 / 13.537 | 5.214 / 4.328 |
| polluted-10min | 1174.724 | 1197.351 | 14.047 / 13.305 | 14.051 / 13.440 | 5.398 / 4.257 |
| clean-10min | 946.076 | 941.531 | 14.342 / 13.792 | 14.335 / 13.811 | 4.907 / 3.766 |

Every scenario has zero NaNs/infinities and zero concentrations below `-1e-20`
mol m-3. Their internal-step totals match the disclosed MICM totals exactly:
48,619; 90,980; 145,891; 106,393; and 84,273 in table order. Agreement with the
high-accuracy references is also the same as MICM's to the reported digits;
the much larger difference there comes from the prescribed method/tolerances.

The additional tests have 65 cells, seed 9127, 192 steps of 900 seconds, latitude
62, day 200, start time 23:00, and NOx ×3. They record steps 25, 50, 75, 100,
125, 150, 175, and 192, including a one-cell final group. The second test shifts
photolysis sample times by 31 seconds and ends table coverage at 36 hours to
check interpolation and extension at both endpoints.

| Additional case | Clang digits: overall / worst | GCC digits: overall / worst |
|---|---:|---:|
| extra | 14.261 / 13.852 | 14.279 / 13.873 |
| edge | 14.251 / 13.815 | 14.262 / 13.872 |

Raw samples, solver statistics, score summaries, build times, and the generated
source hash are retained in `tests/data/tier2_results.json`. README.md gives the
commands to regenerate the solvers and repeat the gates and timing runs.
