# CPU chemistry library

This repository translates the CPU chemistry mechanisms, rate constants, and numerical solvers to ISO C11. It includes a generic CPU library and a Python generator for complete mechanism-specialized three-stage Rosenbrock solvers. Python reads YAML/JSON configurations and generates C when building; the library and executables have no Python dependency.

Build and run the tests:

```sh
make env
make cpu CC=clang
make test-cpu CC=clang
make test-cpu CC=gcc-15
make -j8 bench CC=clang
make -j8 test-specialized CC=clang
make -j8 test-specialized CC=gcc-15
```

`make cpu` creates `bin/libchem.a`, the shared library, `bin/ts1_bench_cpu`, and the C tutorial runner. `make bench` also builds `bin/ts1_bench_specialized` and `bin/ts1_scenario_specialized`. Compilers and optimization flags are configurable, for example `make bench CC=gcc-15 OPTFLAGS='-O3 -DNDEBUG'`. Compiler-specific intermediate directories prevent stale binaries when switching compilers. OpenMP is required. On macOS the Makefile sets the SDK and finds Homebrew libomp with `brew --prefix libomp` (override with `LIBOMP_PREFIX`); other platforms use `-fopenmp`.

The prefix environment `.conda` contains pinned build/test dependencies. All Python commands use `.conda/bin/python`; environment creation confines the package cache to `.conda-pkgs`. The library and executables run without that environment.

Run the benchmark with:

```
bin/ts1_bench_cpu --threads 1 100 100 30 ../reference/ts1_cells_100_seed1.csv
.conda/bin/python ../reference/compare_final.py ts1_final.csv ../reference/micm_ts1_final_100_seed1.csv
bin/ts1_bench_specialized --threads 18 100 100 30 ../reference/ts1_cells_100_seed1.csv
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

Run a day-night scenario with:

```sh
bin/ts1_scenario_specialized --threads 18 --group-size 32 --rtol 1e-6 --atol 1e-15 \
  ../scenarios/diurnal-10min/scenario.txt
.conda/bin/python ../scenarios/score.py ts1_snapshots.csv \
  ../scenarios/diurnal-10min/micm_snapshots.csv --level 10
```

Scenario paths resolve relative to the scenario file. Photolysis factors interpolate linearly and remain constant beyond the table. Every `j` parameter uses its initialized per-cell value times the current factor. This runner has no warm-up; it times rate updates and integration while excluding input, initialization, snapshots, and CSV output. It writes `ts1_final.csv` and, when `snapshot_every` is positive, `ts1_snapshots.csv` including the final step. Defaults are MICM's tolerances; use the tighter values above for day-night comparisons. `--group-size` (also `--group`) and `--threads` are runtime controls; `make GROUP_SIZE=N bench` sets the default group size.

`make test-specialized` generates, compiles, and tests TS1, all three Tier 1 configuration examples, the all-rate-family validation mechanism, a source/sink/duplicate-term fixture, and an automatically modified TS1 with a new reaction. It compares rate constants, forcing, Jacobians, LU, triangular solves, and three integrated intervals at `1e-30 + 1e-10 * abs(generic)` for 1, 7, 32, and 33 cells, changing conditions and custom parameters between intervals. Any failure exits nonzero.

To generate another solver:

```sh
.conda/bin/python tools/specialize.py mechanism.yaml generated/example_specialized.c \
  --name example_specialized
```

Compile that build product with the Makefile's `$(CC)` and `$(OPTFLAGS)`. Its entry points are `chem_load_example_specialized(1)`, `chem_rates_example_specialized(state)`, and `chem_solve_example_specialized(dt, state, stats)`. Create and populate `ChemState` with the generic C API, calculate rates after changing conditions or custom parameters, and call the generated solver. Its state must belong to the generated mechanism. Tolerances live in the state. Generated C is ignored by git and always regenerated from configuration and generator dependencies; never edit it.

The C API is in `include/chem.h`, `include/chem_array.h`, `include/chem_system.h`, and `include/chem_verify.h`. `tests/tutorials.c` gives complete construction, state, solve, printing, and threading examples. See `DESIGN.md` for configuration generation and tuning, `SPECIALIZATION.md` for specialization experiments, `DECISIONS.md` for translation choices, and `RESULTS.md` for validation and timings.
