# C chemistry solver

This repository contains a portable C11 CPU library, a build-time mechanism
reader/generator, and the TS1 benchmark. The executable and library link only
the C standard library, `libm`, and OpenMP. Python is used only to read YAML or
JSON mechanism files and emit C tables, or to drive the test cases.

## Build

Create the pinned Python environment once, then build the library and TS1
benchmark:

```sh
make env
make cpu
```

`make cpu` produces `lib/libchem.a` and `bin/ts1_bench_cpu`. `make bench` builds
the benchmark. The Makefile defaults to `CC=clang` and `OPTFLAGS=-O3 -DNDEBUG`.
To rebuild with another compiler, clean first so no objects from the prior
compiler are reused:

```sh
make clean
make CC=clang cpu
make clean
make CC=gcc-15 cpu
```

`OPTFLAGS` is passed through as given. On macOS, the Makefile sets the Xcode SDK
for linking and selects the OpenMP flags for Apple Clang or Homebrew GCC.

## Mechanisms

`tools/generate_c.py` reads JSON directly or follows `camp-files` entries in a
MICM YAML or JSON configuration. It turns the species, reactions, rate
parameters, and defaults into C tables. The runtime does not parse files and
does not depend on Python. Build a mechanism archive from any supported
configuration with:

```sh
make mechanism CONFIG=examples/robertson/config.json
```

The archive includes the generic solver and the generated `chem_model` table.
The committed tutorial inputs under `examples/` cover YAML and JSON, built-in
and user-defined rates, and the Robertson mechanism.

`tools/translate_ts1.py` translates the read-only TS1 reaction table and
initial-state values into `mechanisms/ts1.json`. The generated C model is built
as part of `make cpu`.

## Test

```sh
make test-cpu
```

This prints one `PASS`, `FAIL`, or `N/A` line for each of the 680 case IDs in
`../reference/micm_cases.txt`. It compiles the three example configurations,
compares both tutorial mechanisms and ten interval states with hand-built C
mechanisms, and runs the 100-cell TS1 reference comparison. `N/A` is reserved
for compile-time C++ type traits and CUDA-only behavior with no portable CPU C
equivalent.

## Run TS1

The command accepts a CSV with a header row and one row per cell. Unknown values
use TS1's initial state. Temperature, pressure, species, and custom rate
parameter columns are matched by name.

```sh
bin/ts1_bench_cpu [--threads N] [--group-size G] CELLS STEPS DT CELLS_CSV
```

The benchmark performs one untimed warm-up interval, then times `STEPS`
intervals of `DT` seconds. It prints elapsed time and summed solver statistics,
then writes the final concentrations to `ts1_final.csv`. By default it uses
OpenMP's maximum thread count and one contiguous cell group per thread.

For example, using the committed 100-cell reference input:

```sh
bin/ts1_bench_cpu --threads 1 100 100 30 ../reference/ts1_cells_100_seed1.csv
python ../reference/compare_final.py ts1_final.csv ../reference/micm_ts1_final_100_seed1.csv
```
