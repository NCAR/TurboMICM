# CPU chemistry library

This repository translates the CPU chemistry mechanisms, rate constants, and numerical solvers to ISO C11. Python reads YAML/JSON configurations and generates C when building; the library and benchmark have no Python dependency.

Build and run the tests:

```sh
make env
make cpu CC=clang
make test-cpu CC=clang
make test-cpu CC=gcc-15
```

`make cpu` creates `bin/libchem.a`, the shared library, `bin/ts1_bench_cpu`, and the C tutorial runner. `make bench` builds the same binaries. Compilers and optimization flags are configurable, for example `make cpu CC=gcc-15 OPTFLAGS='-O3 -DNDEBUG'`. Compiler-specific intermediate directories prevent stale binaries when switching compilers. OpenMP is required. On macOS the Makefile sets the SDK and finds Homebrew libomp for Apple clang; other platforms use `-fopenmp`.

Tier 3 also builds `bin/ts1_scenario` from both `make cpu` and `make bench`. Run a disclosed scenario with the recommended numerical settings:

```sh
bin/ts1_scenario --method rodas3 --option rtol=3e-4 \
  --option group=16 --option schedule=dynamic \
  ../scenarios/diurnal-10min/scenario.txt
.conda/bin/python ../scenarios/score.py ts1_snapshots.csv \
  ../scenarios/diurnal-10min/reference_snapshots.csv --level 3
```

Omit `--method` for the exact inherited three-stage solver and its default tolerances. That default remains available for compatibility but is inaccurate or fails on some long-step scenarios. `OPERATING_POINTS.md` gives the five declared alternatives and their measured accuracy. New methods are opt-in additions; the original numerical API and benchmark remain available.

The scenario executable reads every cell, applies the piecewise linear photolysis factor at each step's start, and integrates with rates held fixed through that step. It has no warm-up. The printed `elapsed_ms` includes rate updates and integration but excludes reading, setup, validation, and writing. It writes the final state and requested snapshots with 17 significant digits. Use `--threads N` to limit workers; the default uses all OpenMP workers. All options, numerical policies, and extension points are documented in `DESIGN.md`.

Additional validation and reproducible scans:

```sh
make test-tier3 CC=clang
make test-tier3 CC=gcc-15
make bench CC=clang
.conda/bin/python tools/scan.py scans/operating-points-throughput.json --study local-points --repeats 5
make micm-scenario
.conda/bin/python tools/scan.py scans/micm-baselines.json --study local-reference \
  --backend micm --compiler clang++ --repeats 5
```

`test-tier3` compares independent Python prototypes with all new C formulas, including operation counts, and tests native scenario interpolation, partial groups, snapshots, rate caching, and scheduling. `ACCURACY.md` explains the accuracy policy, `SCANS.md` gives the research conclusions, and `scans/README.md` documents the raw data and reproduction. `RESULTS.md` includes work–precision figures and independent MICM comparisons. Python is used only by build, test, prototype, and analysis tools. NumPy and Matplotlib are pinned inside the private environment; the runtime library and benchmarks link neither.

The prefix environment `.conda` contains pinned build/test dependencies. All Python commands use `.conda/bin/python`; environment creation confines the package cache to `.conda-pkgs`. The library and executables run without that environment.

Run the benchmark with:

```
bin/ts1_bench_cpu --threads 1 100 100 30 ../reference/ts1_cells_100_seed1.csv
.conda/bin/python ../reference/compare_final.py ts1_final.csv ../reference/micm_ts1_final_100_seed1.csv
```

It runs one warm-up step, times the requested steps, prints elapsed time and summed solver statistics, and writes `ts1_final.csv`. Omit `--threads` to use all OpenMP workers. `--group N` changes the number of cells sharing an adaptive solver controller; the default is 32. `OMP_SCHEDULE` selects group scheduling, with static scheduling as the benchmark default.

For a larger input:

```sh
mkdir -p results
.conda/bin/python ../reference/make_cells_csv.py 10000 2026 --output results/cells.csv
bin/ts1_bench_cpu 10000 100 30 results/cells.csv
```

The CSV header accepts temperature, pressure, species names, and custom rate labels. Each row overrides that cell's defaults, and air density is recalculated from its temperature and pressure. Extra rows are allowed. Species and parameters absent from the CSV retain TS1's initial values.

`make micm-ref` builds the independent reference executable using the supplied read-only headers. It is separate from the C library and benchmark.

`make test-cpu` reports exactly one result for each of the 680 identifiers in `case-list.txt`, and exits unsuccessfully if any case fails. Python drives the C API through ctypes; the eight program cases run `bin/tutorials CASE`. Their reference stdout is retained in `tests/reference-output`. All three configuration examples are compared with independently code-built mechanisms, with generic and generated kernels. The suite also exercises generated TS1 kernels. Use `.conda/bin/python tests/run.py --filter TEXT --traceback` for focused diagnostics.

The C API is in `include/chem.h`, `include/chem_array.h`, `include/chem_system.h`, and `include/chem_verify.h`. `tests/tutorials.c` gives complete construction, state, solve, printing, and threading examples. See `DESIGN.md` for configuration generation and extension points, `DECISIONS.md` for translation choices, and `RESULTS.md` for validation and timings.
