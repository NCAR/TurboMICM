# CPU chemistry port

This repository contains a C11 chemistry library, build-time configuration
reader, and TS1 benchmark. `make test-cpu` checks the reference cases with
either supported compiler; see [RESULTS.md](RESULTS.md) for coverage and
benchmark measurements.

## Build and test

Create the local, prefix conda environment once, then build with either compiler:

```sh
make env
make cpu CC=clang
make test-cpu CC=clang
make cpu CC=gcc-15
make test-cpu CC=gcc-15
```

`CONDA` can point to a different conda executable. The environment and its
package cache are kept in `.conda` and `.conda-pkgs` inside this directory.
`make env` leaves an existing environment alone. `OPTFLAGS` defaults to
`-O3 -DNDEBUG`; the Makefile uses the supplied `CC` and `OPTFLAGS` without
adding optimization or target flags. On macOS it selects the SDK using
`xcrun`; `LIBOMP_PREFIX` can override the location of clang's OpenMP runtime.

`make test-core` runs the independent analytic, configuration, constraint,
and linear-algebra checks. `make test-cpu` prints one result for each ID in
`tests/cases.txt` and exits unsuccessfully if any applicable case fails.
It also runs the core checks. Both commands use the local Python environment
only for build and test orchestration.

## Benchmark

`make bench` builds the same C benchmark as `make cpu`:

```sh
bin/ts1_bench_cpu [--threads N] [--group-size G] CELLS STEPS DT CELLS_CSV
```

The program reads a header and one row per cell, performs one untimed warm-up
step, times the requested steps, prints the elapsed time and solver statistics,
and writes `ts1_final.csv`. The defaults and tuning controls are described in
[DESIGN.md](DESIGN.md). To compare with the reference inputs in the experiment:

```sh
make micm-ref
.conda/bin/python ../reference/make_cells_csv.py 100 1 --output build/cells.csv
bin/ts1_bench_cpu --threads 1 100 100 30 build/cells.csv
.conda/bin/python ../reference/compare_final.py ts1_final.csv ../reference/micm_ts1_final_100_seed1.csv
```

`make micm-ref` is the only build target that depends on the neighboring MICM
source tree. The library and benchmark do not link MICM, Python, or a third-party
chemistry library.

## Mechanisms

The generator accepts a CAMP-style JSON or YAML entry file, resolves its
`camp-files`, and emits a C definition with a caller-selected symbol:

```sh
.conda/bin/python tools/generate.py path/to/config.yaml build/mechanism.c \
  --symbol my_mechanism --reorder
```

`include/chem.h` documents the C data structures and state API. A mechanism may
also be assembled directly in C. `tools/import_ts1.py` translates the supplied
TS1 source data to `mechanisms/ts1.json`; normal builds use that checked-in
configuration and do not need the source file. Format extensions and unit
choices are in [DECISIONS.md](DECISIONS.md).
