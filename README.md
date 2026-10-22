# CPU chemical kinetics in C

An ISO C11 port of MICM's CPU chemistry, matrix utilities, sparse LU solvers,
Rosenbrock methods, backward Euler, and algebraic constraints. The native library
and TS1 executable use only the C library, libm, and OpenMP. Python is used for
configuration generation and tests, never during a native solve.

`tools/specialize.py` also compiles any supported MICM configuration into fixed
chemistry and sparse linear algebra kernels. The specialized TS1 solver keeps
MICM's three-stage Rosenbrock method, tolerances, and adaptive control. See
[SPECIALIZATION.md](SPECIALIZATION.md) for the experiments and
[RESULTS.md](RESULTS.md) for the measured accuracy and performance.

## Build and test

Run these commands from this repository:

```sh
make env
make cpu CC=clang
make test-cpu CC=clang
make test-cpu CC=gcc-15
make bench CC=clang
make test-specialized CC=clang
make test-specialized CC=gcc-15
```

`make env` creates or updates the prefix environment `.conda` from the pinned
`environment.yml`; its package cache is `.conda-pkgs`. Set `CONDA=/path/to/conda`
if conda is not installed under `~/miniforge3`. All Python commands below use that
environment. No activation or changes to the base environment are needed.

`make cpu` produces `lib/libchem.so` and `bin/ts1_bench_cpu`. `make bench` also
builds `bin/ts1_bench_specialized` and `bin/ts1_scenario_specialized` from generated
C. The library filename is also `.so` on macOS. Compiler and
flag changes automatically rebuild the C objects. `OPTFLAGS` defaults to
`-O3 -DNDEBUG`; the Makefile adds language, warning, include, OpenMP, and link
flags only. On macOS it selects the Xcode SDK with `xcrun` and uses Homebrew's
libomp for Apple clang (`brew --prefix libomp`, overridable with `LIBOMP_PREFIX`).
On other systems the compiler must support `-fopenmp`.

`make test-specialized` generates and tests TS1, all three configuration-file
tutorials, a mechanism covering every supported rate family, and a mechanism
without reactions. It compares rates, forcing, Jacobians, LU, direct solves,
and evolving solutions against the generic library to `1e-10` relative, at both
default and tight tolerances, with multiple cell counts and changed runtime data.
Any discrepancy produces a nonzero exit status.

`make test-cpu` reports each of the 680 supplied case IDs individually and exits
nonzero on any failure. Cases cover native code through ctypes and standalone C
programs. No MICM source, executable, or Python package is needed by this target.
To rerun a subset after building:

```sh
.conda/bin/python tests/run_cases.py --filter 'rosenbrock:*' --traceback
```

See [RESULTS.md](RESULTS.md) for results from both compilers and the TS1 gate.

## Run TS1

```sh
make bench CC=clang
bin/ts1_bench_cpu --threads 1 100 100 30 ../reference/ts1_cells_100_seed1.csv
bin/ts1_bench_specialized --threads 18 100 100 30 ../reference/ts1_cells_100_seed1.csv
.conda/bin/python ../reference/compare_final.py ts1_final.csv ../reference/micm_ts1_final_100_seed1.csv
```

The command syntax is:

```text
bin/ts1_bench_cpu [--threads N] [--group-size G] CELLS STEPS DT CELLS_CSV
bin/ts1_bench_specialized [--threads N] [--group-size G] CELLS STEPS DT CELLS_CSV
```

The CSV header names `temperature`, `pressure`, species, or custom rate
parameters. The first `CELLS` rows are used; omitted columns take the translated
TS1 initial values. Air density is recomputed as pressure / (R temperature).
There is one untimed warm-up solve, followed by `STEPS` timed solves of `DT`
seconds. The executable prints elapsed milliseconds and summed solver counters,
then writes `ts1_final.csv` in its working directory. Setup and CSV writing are
outside the timer. Use separate working directories to retain multiple runs.

By default OpenMP uses all available cores. `--threads` overrides its thread
count, subject to the number of groups. `OMP_NUM_THREADS`, `OMP_SCHEDULE`, and
`CH_GROUP_SIZE` are also supported; see [DESIGN.md](DESIGN.md) for defaults and
the scope of each parameter.

Run a day-night scenario with:

```sh
bin/ts1_scenario_specialized --threads 18 --group-size 32 --rtol 1e-6 --atol 1e-15 ../scenarios/diurnal-10min/scenario.txt
.conda/bin/python ../scenarios/score.py ts1_snapshots.csv ../scenarios/diurnal-10min/micm_snapshots.csv --level 10
```

The scenario file supplies relative paths to the cells and photolysis CSVs,
`dt`, `steps`, and optional `snapshot_every`. Each step resets every parameter
whose name begins with `j` to its original cell value times the interpolated
factor. There is no warm-up. Outputs are `ts1_final.csv` and
`ts1_snapshots.csv`, with the final step always present. Defaults are MICM's
relative `1e-6` and absolute `1e-3` tolerances; the day-night gate uses `1e-15`
absolute tolerance as shown above. Snapshot copies and CSV output are untimed.

To reproduce a full comparison with a reference built using the same flags:

```sh
make cpu micm-ref CC=clang OPTFLAGS='-O3 -DNDEBUG'
mkdir -p runs/c runs/micm
.conda/bin/python ../reference/make_cells_csv.py 10000 2026 --output runs/cells.csv
(cd runs/c && ../../bin/ts1_bench_cpu 10000 100 30 ../cells.csv > timing.txt)
(cd runs/micm && ../../bin/micm_ts1_ref 10000 100 30 ../cells.csv > timing.txt)
.conda/bin/python ../reference/compare_final.py runs/c/ts1_final.csv runs/micm/ts1_final.csv
```

Only `make micm-ref` compiles C++ or uses MICM headers. It needs the supplied
read-only `../reference` and `../micm` inputs. Normal C builds are independent of
those directories. The committed configuration is enough to rebuild TS1.

## Configure and embed

The independent configuration reader accepts the YAML/JSON `camp-files` and
`camp-data` format from the tutorial examples. For example:

```sh
.conda/bin/python tools/configure.py configs/examples/robertson/config.json build/generated/robertson.c --name robertson
```

The generated functions `ch_build_robertson()` and `ch_init_robertson(state)`
construct and initialize ordinary library objects. Compile generated C with the
same `CC`, `OPTFLAGS`, and include flags as other C sources. All three supplied
configuration examples are built and checked against independent hand-built
mechanisms by `make test-cpu`.

For specialized numerical kernels:

```sh
.conda/bin/python tools/specialize.py configs/examples/robertson/config.json build/generated/custom.c --name custom
```

This writes `custom.c`, `custom_kernels.c`, `custom.h`, and `custom.json` (symbolic
metadata). Compile both C files with the same `CC` and `OPTFLAGS`, includes
`-Iinclude -Isrc -Itools`, and OpenMP flags used by the Makefile, and link the
generic library objects. Call `ch_build_custom`, `ch_state_create`,
`ch_init_custom`, `ch_rates_custom`, and `ch_solve_custom` with
`ch_rosenbrock_default(3)`. Set conditions, custom parameters, concentrations,
and tolerances through the ordinary state API. Recompute rates after changing
conditions or parameters. Use only states built from the matching generated
model; changes to fixed reaction data require regeneration.

Generated C is ignored by git and is never edited. Make regenerates TS1 and the
test mechanisms when their configuration, compiler generator, or generation
options change. `make clean bench` demonstrates complete regeneration. The
tuning controls and fallback behavior are listed in [DESIGN.md](DESIGN.md).

The retained translation command for TS1 is:

```sh
.conda/bin/python tools/translate_ts1.py ../micm/test/benchmark/ts1_mechanism.hpp configs/ts1.json
```

Public APIs are in `include/chem.h` and `include/chem_array.h`. A typical caller
creates a model, adds species/parameters/reactions, finalizes it, creates one
state per cell group, sets conditions and concentrations, calls
`ch_update_rates`, then calls `ch_solve`. Free states before their shared model.
`tools/tutorials.c` contains complete native examples. [DESIGN.md](DESIGN.md)
documents ownership, numerical conventions, and extension points.

The numerical methods and adapted test fixtures retain MICM's Apache-2.0
attribution. See `LICENSE` and `tests/reference/README.md`.
