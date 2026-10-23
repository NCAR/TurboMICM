# Specialized portable chemistry solver

This repository contains Tier 1's generic ISO C11 chemistry library and a
Python build-time generator that emits fixed TS1 reaction, Jacobian, sparse LU,
and triangular-solve schedules. Runtime programs are C and use only the C
library, libm, and CPU OpenMP.

## Build and checks

The compiler is selected with `CC`; `OPTFLAGS` defaults to `-O3 -DNDEBUG`.
Apple clang and GCC 15 select their OpenMP flags in the Makefile.

```sh
make env
make clean
make bench CC=clang
make test-cpu CC=clang
make test-specialized CC=clang
make clean
make bench CC=gcc-15
make test-cpu CC=gcc-15
make test-specialized CC=gcc-15
```

`make env` creates the pinned build-time Python environment at `.conda/` and
keeps its package cache in `.conda-pkgs/`. Python and PyYAML are used only to
read configuration files and generate C.

`make test-cpu` runs the Tier 1 C checks and all cases in the supplied MICM
inventory. `make test-specialized` compares the TS1 specialized and generic
solvers and exercises generated forcing and Jacobian code on the two tutorial
configuration mechanisms.

## Configuration and regeneration

`tools/read_config.py` resolves CAMP `camp-files` paths and writes normalized
JSON. `tools/import_ts1.py` translates the permitted TS1 builder input to
`mechanisms/ts1.json`. `tools/generate_c.py mechanisms/ts1.json generated/ts1.c`
recreates the C build product. The generated source always includes the generic
mechanism tables. Defining `CHEM_SPECIALIZED` enables the generated TS1
evaluation and linear-algebra schedules; the Makefile sets it only for the two
specialized executables.

## TS1 benchmark

```sh
bin/ts1_bench_cpu [--threads N] CELLS STEPS DT CELLS_CSV
bin/ts1_bench_specialized [--threads N] CELLS STEPS DT CELLS_CSV
```

The CSV contains temperature, pressure, species, and custom rate-parameter
columns; omitted values use TS1 defaults. Both programs run one untimed warm-up
step, time `STEPS` steps of `DT` seconds with MICM's default tolerances, print
`elapsed_ms=` and solver statistics, and write `ts1_final.csv`.
`bin/ts1_bench_cpu` uses the generic implementation and
`bin/ts1_bench_specialized` opts into generated TS1 schedules.

## Day and night scenarios

```sh
bin/ts1_scenario_specialized [--threads N] [--group-size N] \
    [--rtol R] [--atol A] SCENARIO
```

The runner reads the scenario file and its relative `cells` and `photolysis`
paths, interpolates the photolysis factor at each step, and scales every
`j*` parameter from the original cell input. It writes `ts1_final.csv` and
`ts1_snapshots.csv`. The default group size creates one contiguous cell group
per worker; `--group-size 32` selects 32-cell groups. Default tolerances match
MICM; the scenario gate uses `--rtol 1e-6 --atol 1e-15`.

Compare benchmark output with the supplied MICM reference:

```sh
make micm-ref
python3 ../reference/compare_final.py ts1_final.csv ../reference/micm_ts1_final_100_seed1.csv
python3 ../scenarios/score.py ts1_snapshots.csv ../scenarios/diurnal-30min/micm_snapshots.csv --level 10
```

## Generic C API

`include/chemistry.h` defines the generic mechanism, rate, forcing, Jacobian,
matrix, constraint, and solver APIs. The generic Tier 1 path remains available
for arbitrary mechanisms and is built as `bin/ts1_bench_cpu`. See `DESIGN.md`
for the specialization boundary and `SPECIALIZATION.md` for measured effects.
`RESULTS.md` records validation and benchmark results for both compilers.
