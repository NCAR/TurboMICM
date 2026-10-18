# Portable chemistry library

C11 CPU chemistry with generic mechanisms, Rosenbrock and backward Euler solvers,
algebraic constraints, external models, sparse LU, and grouped matrix utilities.
All 680 supplied reference cases are ported; no cases are marked N/A.

Tier 3 adds adaptive Rosenbrock, linearly implicit Euler extrapolation, and
SDIRK integration, with scenario-based accuracy and cost experiments. The
original integration paths and their roundoff compatibility remain available.
See [OPERATING_POINTS.md](OPERATING_POINTS.md) for evaluated configurations,
[ACCURACY.md](ACCURACY.md) for the accuracy policy, and [SCANS.md](SCANS.md) for
the controlled experiments. [RESULTS.md](RESULTS.md) contains measurements and
work–precision plots.

## Build and test

```sh
make env
make test-cpu CC=clang
make test-cpu CC=gcc-15
make cpu CC=clang
make bench CC=clang
make test-tier3 CC=clang
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

`make micm-ref` builds the original gate oracle from the supplied reference headers;
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

## Run a time-dependent TS1 scenario

```sh
make bench
bin/ts1_scenario --method extrap --option order=4 --option rtol=0.002 \
    --option atol=1e-18 ../scenarios/diurnal-30min/scenario.txt
.conda/bin/python ../scenarios/score.py ts1_snapshots.csv \
    ../scenarios/diurnal-30min/reference_snapshots.csv --level 3
```

`--threads N` limits OpenMP threads; the default uses all cores. `--option
group-size=32` selects the cell grouping. Omitting `--method` selects the
original MICM three-stage solver with its original tolerances. This default is
provided for compatibility; it fails or is inaccurate on some long scenarios.
All method names, options, formulas, cache rules, and counters are documented
in [DESIGN.md](DESIGN.md).

The runner reads the documented `key = value` scenario file, updates every
photolysis parameter at each external step, and freezes rates within that step.
There is no warm-up. It prints `elapsed_ms` and work counters and writes
`ts1_final.csv` and, when requested, `ts1_snapshots.csv` in the current directory.
Only rate updates, integration and validation are timed; CSV input/output and
setup are excluded. Rate caching, four-step blocks and dynamic group scheduling
are enabled by default; all remain runtime options. Caching recalculates every
affected photolysis rate, and blocking preserves every external forcing boundary. Use separate working directories to retain multiple runs.

## Reproduce the research

```sh
make env
make bench test-tier3 micm-scenario CC=clang
.conda/bin/python tools/scan.py scans/operating_points.json \
    --scan local-operating-points --threads 10 --repeats 5
.conda/bin/python tools/stress.py --scan local-stress-references
.conda/bin/python tools/scan.py scans/operating_points.json \
    --scan local-stress --scenario-root build/stress \
    --scenarios summer48 winter48 tropical48 --threads 4
```

Scan names must be new: scripts refuse to overwrite recorded CSVs. Source
configurations and every scored run are in `scans/`; raw outputs stay in
`build/scan-runs/`. Small-case prototypes, order checks and stiff-decay tests are
in `prototypes/`. Python dependencies are pinned in `environment.yml` and are
used only for development and analysis. The C library and both C benchmarks
never invoke Python. `make env` creates or updates only the local environment
and uses the local package cache. The optional `micm-scenario` target compiles
the supplied read-only C++ oracle; it is not part of the C library.

For the final compiler gates and larger timing domain:

```sh
.conda/bin/python tools/validate.py --compiler gcc-15
.conda/bin/python tools/validate.py --compiler clang
.conda/bin/python tools/scaling.py
.conda/bin/python tools/scan.py scans/operating_points.json \
    --scan local-10000 --scenario-root build/scaling --threads 10 --repeats 5 \
    --discard-repeats --reuse-identical-scores
.conda/bin/python tools/report.py
```

The 10,000-cell timing inputs have independently supplied references only for
their first 1,000 cells. The scorer checks finite/nonnegative values on all
10,000 cells, but the accuracy claim covers the reference prefix. The report
script reads the retained final scan names and regenerates tables and plots.
Run builds and timing scans sequentially; changing compilers during a scan would
invalidate its performance comparison.
