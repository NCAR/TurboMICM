# Generated mechanism solver

This repository contains a portable C11 chemistry library, a generic mechanism
reader, and generated TS1 solvers. The library and executables use only the C
standard library, `libm`, and CPU OpenMP. Python runs during generation and
testing; generated programs do not need Python.

## Build

Create the pinned local Python environment, then build the Tier 1 library and
both TS1 benchmark programs:

```sh
make env
make cpu
make bench
```

`make cpu` builds `lib/libchem.a` and `bin/ts1_bench_cpu`. `make bench` also
builds `bin/ts1_bench_specialized` and `bin/ts1_scenario_specialized`. The
default compiler is Apple Clang; `OPTFLAGS` defaults to `-O3 -DNDEBUG`. Clean
between compiler changes so objects from different compilers are not reused:

```sh
make clean
make CC=clang bench
make clean
make CC=gcc-15 bench
```

`make test-cpu` runs the Tier 1 MICM case list and tutorial checks.
`make test-specialized` generates and compares TS1 and all three Tier 1 example
mechanisms against the generic solver at `1e-10` relative tolerance.

## Mechanism generation

`tools/generate_c.py` reads a normalized JSON mechanism or a MICM YAML/JSON
configuration with `camp-files`. The ordinary build emits immutable mechanism
tables for the generic library. The `--specialized` mode also emits fixed rate
formulas, reaction and Jacobian updates, sparse LU operations, and triangular
solves. All generated C is a build product. For example:

```sh
make mechanism CONFIG=examples/robertson/config.json
make test-specialized
```

The committed tutorial inputs cover YAML and JSON, built-in and user-defined
rates, and Robertson chemistry. `tools/translate_ts1.py` converts MICM's
read-only TS1 reaction table and initial state to `mechanisms/ts1.json`.

The public C API also supports linear, mass-action equilibrium, and
callback-defined DAE constraints. Generated kernels are used when their fixed
sparse pattern matches the solver. A mechanism with runtime constraints uses
the generic constraint-aware path.

## TS1 benchmark

Both benchmark commands accept a CELLS_CSV with a header and one row per cell.
Unknown cells use TS1's initial state. Temperature, pressure, species, and
custom rate parameter columns are matched by name.

```sh
bin/ts1_bench_cpu [--threads N] [--group-size G] CELLS STEPS DT CELLS_CSV
bin/ts1_bench_specialized [--threads N] [--group-size G] CELLS STEPS DT CELLS_CSV
```

Each runs one untimed warm-up step, followed by `STEPS` timed intervals of
`DT` seconds with MICM's default tolerances. The programs print elapsed time and
summed solver statistics, then write `ts1_final.csv`. By default, OpenMP uses
its maximum thread count and creates one contiguous cell group per thread. A
fixed run-time group size is selected with `--group-size`.

For example, using the committed 100-cell reference input:

```sh
bin/ts1_bench_specialized --threads 1 100 100 30 ../reference/ts1_cells_100_seed1.csv
python3 ../reference/compare_final.py ts1_final.csv ../reference/micm_ts1_final_100_seed1.csv
```

## Day-night scenarios

`bin/ts1_scenario_specialized` reads the scenario format described in
`../scenarios/README.md`:

```sh
bin/ts1_scenario_specialized --threads 18 --group-size 32 \
  --rtol 1e-6 --atol 1e-15 ../scenarios/diurnal-10min/scenario.txt
```

It updates every `j*` rate from the cell's noon value and the interpolated
photolysis factor at each step. It writes `ts1_final.csv` and
`ts1_snapshots.csv`, and times the photolysis update and integration work.
Scenario snapshots include every `snapshot_every` interval and the final
interval. See `RESULTS.md` for the measured benchmark, scenario scores, and
reference comparisons.
