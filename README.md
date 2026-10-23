# Portable C chemistry solvers

This repository contains the Tier 1 mechanism-generic ISO C11 chemistry library
and a Tier 2 Python generator that emits a solver specialized to a MICM-format
mechanism. Generated C is a build product. Python and PyYAML run only while
building; neither solver links Python.

## Build and test

Create the pinned local Python environment and build with either supported
compiler:

```sh
make env
make cpu CC=clang
make cpu CC=gcc-15
```

`OPTFLAGS` defaults to `-O3 -DNDEBUG`. The Makefile selects OpenMP flags for
Apple Clang and GCC. The Tier 1 regression gate is:

```sh
make test-cpu CC=clang
make test-cpu CC=gcc-15
```

It prints a PASS/FAIL/N/A line for each MICM case in the reference list and exits
nonzero on failure. The Tier 2 generator checks are:

```sh
make test-specialized CC=clang
make test-specialized CC=gcc-15
```

They generate TS1 and the two configuration mechanisms used by Tier 1 tests,
then compare each generated solver with the generic C library to `1e-10`
relative tolerance.

`make env` pins Python and PyYAML and puts conda's package cache in this
repository. The library and benchmarks have no Python runtime dependency.

## Tier 1 generic library

`include/chem.h` exposes the integration API, owning species/phase/system
definitions for code-built mechanisms, an owning `chem_state`, dense and sparse
linear algebra, reaction forcing, and Rosenbrock integration. The additional
headers provide constraints, strided group views, dense matrices, sparse
matrices, and Jacobian verification helpers. See [DESIGN.md](DESIGN.md) for
storage and ownership rules. `bin/ts1_bench_cpu` and `make test-cpu` continue to
use this generic library.

## Specialized TS1 benchmark and scenarios

`make bench` builds the generic benchmark and the specialized TS1 programs:

```sh
bin/ts1_bench_specialized [--threads N] CELLS STEPS DT CELLS_CSV
bin/ts1_scenario_specialized [--threads N] [--rtol R] [--atol A] [--group-size G] SCENARIO
```

The benchmark uses MICM's default tolerances, one untimed warm-up interval, and
the Tier 1 output format. The scenario runner reads a scenario file and its
relative CSV paths, sets every `j*` photolysis parameter at the beginning of
each step, and writes `ts1_final.csv` and `ts1_snapshots.csv`. It defaults to
MICM's tolerances and 32-cell groups. Use `--group-size 32 --rtol 1e-6 --atol
1e-15` for the disclosed-scenario gate settings. Scenario definitions and the
scoring tool are in the read-only `../scenarios` directory.

```sh
bin/ts1_scenario_specialized --threads 18 --group-size 32 \
  --rtol 1e-6 --atol 1e-15 ../scenarios/diurnal-30min/scenario.txt
python3 ../scenarios/score.py ts1_snapshots.csv \
  ../scenarios/diurnal-30min/micm_snapshots.csv --level 10
```

The TS1 descriptor is reproducibly translated by `tools/translate_ts1.py` from
the supplied MICM mechanism builder. `tools/generate.py` reads that descriptor
or a MICM-format YAML/JSON configuration and emits mechanism data and the
specialized C solver. Run `make bench` after changing the configuration or
generator to regenerate and compile the output.

See [DESIGN.md](DESIGN.md) for generator internals and parameters,
[SPECIALIZATION.md](SPECIALIZATION.md) for measured specialization effects,
[DECISIONS.md](DECISIONS.md) for implementation choices, and
[RESULTS.md](RESULTS.md) for verification and timing results.
