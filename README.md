# CPU chemistry library

An ISO C11 port of MICM's CPU chemistry functionality. Numerical
formulas and test cases derive from MICM; attribution is in NOTICE and LICENSE.

```
make env
make cpu CC=clang
make test-cpu CC=clang
make cpu test-cpu bench CC=gcc-15
```

The pinned Python prefix environment and package cache live entirely in this
directory. `OPTFLAGS` defaults to `-O3 -DNDEBUG`. OpenMP flags are selected for the
compiler and operating system. Changes to the compiler or flags invalidate the
compiled C artifacts automatically. On Linux, use `CC=gcc` or `CC=clang` with an
installed OpenMP runtime. No third-party runtime library other than OpenMP is
required by the C implementation.

`make test-cpu` prints all 680 reference case identifiers and checks all ported cases, including eight tutorial stdout comparisons; it also
runs configuration equivalence checks. Both compilers currently pass all 680
cases. See RESULTS.md for numerical gates and benchmark measurements.
`make test-config` runs the configuration checks separately. Run an individual
tutorial with `bin/tutorial_cpu README_example` (or another program identifier
from `tests/cases.txt`).

Build and run the TS1 benchmark:

```
make bench
bin/ts1_bench_cpu [--threads N] CELLS STEPS DT CELLS_CSV
```

By default OpenMP uses all available cores. `--threads` limits the thread count.
Set `CHEM_GROUP_SIZE` to change the default group size of 32 and `OMP_SCHEDULE`
to choose static, dynamic or guided group scheduling (default: static), for example:

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

Reproduce the extended benchmark validation after building the desired compiler:

```
make micm-ref
.conda/bin/python tools/validate_benchmarks.py --compiler clang
```

Use `--compiler gcc` after a GCC build. The label selects the result directory;
it does not rebuild the executable. The driver checks the 10,000-cell gate and
smaller inputs across three seeds, varied group sizes and OpenMP schedules,
retaining results under `build/validation/final-<compiler>/`.

## Generated solvers

The Python generator now emits mechanism-specific C chemistry, sparse LU and
triangular solves, alongside the original generic library. Build and validate:

```
make env
make -j6 cpu bench test-cpu test-specialized CC=clang
make -j6 cpu bench test-cpu test-specialized CC=gcc-15
```

`make bench` builds all three drivers. The specialized benchmark has the same
arguments, warm-up, timed region, output file and statistics as the generic one:

```
bin/ts1_bench_specialized --threads 18 10000 100 30 cells.csv
bin/ts1_scenario_specialized --threads 18 --group-size 32 --rtol 1e-6 --atol 1e-15 ../scenarios/diurnal-10min/scenario.txt
```

The scenario driver reads the scenario's relative paths, interpolates the light
factor, scales every `j*` custom parameter from its original cell value, and runs
without warm-up. It writes `ts1_final.csv` and `ts1_snapshots.csv`. Its timer
includes parameter updates and integration, and excludes input and snapshot I/O.
Omitting tolerances selects MICM's defaults; the command above selects the
required day-night settings.

The default group size is 32. `make GROUP_SIZE=32 bench` selects the build-time
default and generated stride. `CHEM_GROUP_SIZE` overrides the driver's group size;
the scenario also accepts `--group-size`. A different runtime size remains valid
and uses the generic chemistry/algebra kernels for that group. Tail groups hold
only real cells, so padding never enters the error norm.

Generate another mechanism with the same CAMP format as the Tier 1 reader:

```
.conda/bin/python tools/specialize.py path/to/config.yaml build/specialized/example.c --symbol example --width 32 --depfile build/specialized/example.d
make build/specialized/example.o
```

Compile that file with the Makefile's C11 flags and link `build/libchem.a`.
Call `example_specialized_create`, `chem_state_create`, `example_initialize`,
`chem_solver_create`, and `example_specialized_attach`, in that order. Then call
`example_update_rates` and `example_solve` in place of the generic rate update and
solve. Their signatures match the corresponding generic functions; attach takes
one `ChemSolver *`. The generated mechanism is immutable and already reordered.
The original `example_create` remains available as a generic builder.

Generated C and dependency files live under `build/` and are never edited.
Changes to a configuration, any transitive CAMP include, generator, compiler,
flags or generation settings trigger the appropriate rebuild. See DESIGN.md for
all tuning parameters and SPECIALIZATION.md for the measured experiments.

`make test-specialized` checks TS1, all three tutorial configurations, the
additional-rate configuration, and a regenerated mutation fixture against the
generic library at `1e-10` relative tolerance. It checks rates, exact-input cache
invalidation, factors, triangular solves, integration and statistics, using
1, 7, 16, 32 and 33 cells and both default and tight absolute tolerances.

Reproduce the full benchmark and disclosed-scenario gate for an existing build:

```
make micm-ref
.conda/bin/python tools/validate_specialized.py --compiler clang
```

Use `--compiler gcc-15` after rebuilding with GCC. The argument labels the results;
it does not select or rebuild the compiler. Logs, CSVs, scores against both MICM
and the high-accuracy references, and a JSON summary are retained in `validation/`.
