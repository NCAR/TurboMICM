# Portable chemistry library

C11 CPU chemistry with generic mechanisms, Rosenbrock and backward Euler solvers,
algebraic constraints, external models, sparse LU, and grouped matrix utilities.
All 680 supplied reference cases are ported; no cases are marked N/A.

## Build and test

```sh
make env
make test-cpu CC=clang
make test-cpu CC=gcc-15
make cpu CC=clang
```

The runtime uses only the C standard library, libm, and OpenMP. Python is used
only for configuration reading, code generation, and development checks.
`make env` creates `.conda` from the pinned, platform-independent
`environment.yml`, with a local package cache. `CONDA` overrides the conda
executable path if needed.

`CC` selects the C compiler; `OPTFLAGS` defaults to `-O3 -DNDEBUG`. Changing either
invalidates compiled objects. On macOS, Make obtains SDKROOT with `xcrun` and
selects Homebrew libomp for Clang or `-fopenmp` for GCC. On other systems Clang
and GCC use `-fopenmp` and the installed OpenMP runtime. No architecture or extra
optimization flags are added.

`make test-cpu` reports each original case ID followed by PASS/FAIL and a reason.
Use `TEST_FILTER=substring` to select cases. It includes all numerical/layout
variants, generated-kernel regressions, configuration equivalence, and the eight
tutorial output comparisons. The elapsed-time fields and order of thread
announcements are validated separately; see DECISIONS.md.

## Run TS1

```sh
make bench CC=clang
bin/ts1_bench_cpu --threads 4 100 100 30 ../reference/ts1_cells_100_seed1.csv
```

Syntax: `bin/ts1_bench_cpu [OPTIONS] CELLS STEPS DT CELLS_CSV`.
Options are `--threads N`, `--group-size G`, `--schedule static|dynamic|guided`,
`--chunk K`, and `--generic`. By default it uses all processors reported by OpenMP
and groups of 32 cells with static scheduling. CSV columns override TS1's initial defaults by name.
The program performs one untimed warm-up step, times the requested subsequent
steps, prints elapsed milliseconds and accumulated solver statistics, and writes
`ts1_final.csv` in the working directory. `--generic` selects the table evaluator
instead of automatically generated chemistry kernels.

```sh
make micm-ref
.conda/bin/python tools/check_ts1.py --output build/ts1-checks
```

`make micm-ref` is the only target that includes the supplied reference headers;
it builds the C++ reference with the same `OPTFLAGS`. The check script tests the
C executable against the supplied 100-cell reference using nine thread/group, scheduling
and evaluator combinations. To generate another reference, run
`bin/micm_ts1_ref --threads 1 CELLS STEPS DT CELLS_CSV` in a separate directory
inside this checkout; it also writes `ts1_final.csv`. Compare saved files with
`.conda/bin/python ../reference/compare_final.py YOURS.csv MICM.csv`.
Reference paths in Make can be overridden with `MICM` and `REFERENCE`.

## Mechanisms and examples

`build/libchem.a` and the headers in `include/` are the library. Low-level
`chem_mechanism` descriptions borrow immutable arrays. The higher-level
`chem_solver.h` API owns copied systems and reactions and supports concurrent
states. See `examples/tutorials.c` for complete C examples; for instance:

```sh
bin/tutorial README_example
bin/tutorial solver_configuration
bin/tutorial openmp_tutorial
```

All eight plain-program IDs from `tests/cases.txt` are accepted by `bin/tutorial`.
To emit another mechanism from YAML or JSON:

```sh
.conda/bin/python tools/generate.py path/to/config.yaml generated/example.c --symbol example_mechanism
```

Compile the emitted source with the same C flags as the library and reference
`extern const chem_mechanism example_mechanism` from your application.
`--no-kernels` emits only the description, using the generic evaluator.
The three tutorial configuration directories in `configs/` are tested against
independent hand-built mechanisms.

The retained TS1 configuration can be regenerated without compiling its source:

```sh
.conda/bin/python tools/translate_ts1.py ../micm/test/benchmark/ts1_mechanism.hpp configs/ts1.json
```

DESIGN.md describes ownership, layouts, tuning parameters, and the extension
points for replacing a controller, solver, linear backend, or rate type.
RESULTS.md records validation and performance measurements.

## Specialized solvers and day-night runs

```sh
make bench test-specialized CC=clang
bin/ts1_bench_specialized --threads 18 10000 100 30 CELLS_CSV
bin/ts1_scenario_specialized --threads 18 --group-size 32 --rtol 1e-6 --atol 1e-15 \
  ../scenarios/diurnal-10min/scenario.txt
```

The specialized benchmark accepts the same options and has the same warm-up,
timing boundary, output and statistics as `ts1_bench_cpu`. `--generic` selects
the generic library for cross-checking. The scenario runner has no warm-up. It
reads cells and a piecewise-linear photolysis table relative to the scenario
file, updates all `j*` parameters at every step, and writes `ts1_final.csv` and
`ts1_snapshots.csv`. Snapshot I/O is excluded from its elapsed time. Relative and
absolute tolerances default to MICM's; the strict values above are appropriate
for the supplied day-night comparisons.

`make test-specialized` generates and compiles TS1, Chapman, all three tutorial
configuration mechanisms, an additional rate/duplicate-term/isolated-species
fixture, and a mechanism with no reactions. It compares rates, forcing,
Jacobians, factorization, triangular solutions and integration with the generic
library at relative tolerance 1e-10, using multiple group sizes and tolerances.
It fails if any comparison fails. Both Clang and GCC use the same target.

```sh
.conda/bin/python tools/check_scenarios.py --output build/day-night
.conda/bin/python tools/check_ts1.py --binary bin/ts1_bench_specialized \
  --output build/specialized-gate
```

To generate another complete solver:

```sh
.conda/bin/python tools/specialize.py path/to/config.yaml generated/example.c --symbol example
```

Compile the generated C with the library and reference
`extern const chem_specialized example` from `chem_specialized.h`. Create a
Markowitz/Mozart model from `example.mechanism`, create its state, set runtime
inputs, call `example.rates(&state)`, then `example.solve(&state, dt)`. See
DESIGN.md for ownership and build parameters, and SPECIALIZATION.md for the
specializations and measured alternatives. Python is never invoked by any
solver, benchmark or scenario executable.
