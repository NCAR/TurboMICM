# C CPU chemistry port

This repository contains an ISO C11 chemistry library, a TS1 benchmark driver,
and Python build tooling for generating C mechanism descriptors from MICM-style
YAML or JSON configuration files. Python is used only while building and testing;
the library and benchmark link no Python runtime.

## Build

Create the local, pinned Python environment once, then build the C library and
benchmark with either supported compiler:

```sh
make env
make cpu CC=clang
make cpu CC=gcc-15
```

`OPTFLAGS` defaults to `-O3 -DNDEBUG`. On macOS the Makefile selects the SDK and
OpenMP flags for Apple Clang; GCC uses `-fopenmp`. Build the test executables and
print one PASS/FAIL/N/A line for every case in the supplied MICM case list with:

```sh
make test-cpu CC=clang
make test-cpu CC=gcc-15
```

The runner exits nonzero while any case fails. The current coverage and known
unported cases are recorded in [RESULTS.md](RESULTS.md).

## Run TS1

The benchmark accepts a CSV with a header naming temperature, pressure, species,
and rate parameters. Unspecified values use TS1 defaults. It performs one
untimed warm-up interval, then the requested timed intervals, prints elapsed
milliseconds and summed solver statistics, and writes `ts1_final.csv` in the
current directory.

```sh
bin/ts1_bench_cpu [--threads N] CELLS STEPS DT CELLS_CSV
```

For example, from the repository root:

```sh
bin/ts1_bench_cpu --threads 4 100 100 30 ../reference/ts1_cells_100_seed1.csv
```

`make bench` builds the benchmark. `make micm-ref` builds the reference program
from the read-only MICM checkout for comparison. `tools/translate_ts1.py` keeps
the TS1 source-to-configuration translation reproducible; `tools/generate.py`
turns TS1 or either tutorial example configuration into a C descriptor.

See [DESIGN.md](DESIGN.md) for the architecture and extension points,
[DECISIONS.md](DECISIONS.md) for format choices, and [RESULTS.md](RESULTS.md)
for current verification data.
