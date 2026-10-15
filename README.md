# CPU chemistry library

An ISO C11 port in progress of MICM's CPU chemistry functionality. Numerical
formulas and test cases derive from MICM; attribution is in NOTICE and LICENSE.

```
make env
make cpu CC=clang
make test-cpu CC=clang
make cpu test-config CC=gcc-15
```

The pinned Python prefix environment and package cache live entirely in this
directory. `OPTFLAGS` defaults to `-O3 -DNDEBUG`. OpenMP flags are selected for the
compiler and operating system. Changes to the compiler or flags invalidate the
compiled C artifacts automatically. On Linux, use `CC=gcc` or `CC=clang` with an
installed OpenMP runtime. No third-party runtime library other than OpenMP is
required by the C implementation.

`make test-cpu` prints all 680 reference case identifiers and reports unported
cases as FAIL; it also runs configuration equivalence checks. The task is not
complete. See RESULTS.md and WORK_LOG.md for verified progress and remaining work.
`make test-config` runs the configuration checks separately.

Build and run the TS1 benchmark:

```
make bench
bin/ts1_bench_cpu [--threads N] CELLS STEPS DT CELLS_CSV
```

By default OpenMP uses all available cores. `--threads` limits the thread count.
Set `CHEM_GROUP_SIZE` to change the default group size of 32, for example:

```
CHEM_GROUP_SIZE=7 bin/ts1_bench_cpu --threads 4 100 100 30 ../reference/ts1_cells_100_seed1.csv
```

The driver performs one untimed warm-up step, times the requested steps, prints
elapsed milliseconds and summed group statistics, and writes `ts1_final.csv`.
The CSV accepts temperature, pressure, species and custom rate labels. Omitted
columns retain the configuration's initial values. Every input row is used for
its corresponding cell. Input files must use the unquoted CSV format supplied
by the reference tools.

Generate a new mechanism builder:

```
.conda/bin/python tools/configure.py path/to/config.yaml build/mechanism.c --symbol example
```

Compile the generated C with the same Makefile compiler/flags as the library,
then call `example_create` and optionally `example_initialize`. Both YAML and
JSON use the tutorial CAMP format; documented extensions are in DECISIONS.md.
`tools/translate_ts1.py` reproduces `mechanisms/ts1.json` from the provided TS1
header. No C++ source is executed during translation.

`make micm-ref` is the explicitly permitted reference-only C++ build, using
`clang++` and `OPTFLAGS`. The C library and benchmark do not link MICM or Python.
