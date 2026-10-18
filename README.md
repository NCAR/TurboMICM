# CPU chemical kinetics in C

An ISO C11 chemistry library with MICM-compatible CPU solvers and seven added
adaptive ODE formulas for atmospheric chemistry. The retained implementation
includes matrix utilities, sparse LU, Rosenbrock methods, backward Euler, and
algebraic constraints. Native code uses only the C library, libm, and OpenMP.
Python prepares configurations, prototypes, tests and research reports; it never
runs during a native solve.

The recommended atmospheric-model configuration is **balanced RODAS3** with a
200-second internal step cap. See [OPERATING_POINTS.md](OPERATING_POINTS.md) for
the four declared configurations, [ACCURACY.md](ACCURACY.md) for the accuracy
argument, and [RESULTS.md](RESULTS.md) for measured work–precision results.

## Build and test

Run these commands from this repository:

```sh
make env
make cpu CC=clang
make test-cpu CC=clang
make test-cpu CC=gcc-15
```

`make env` creates or updates the prefix environment `.conda` from the pinned
`environment.yml`; its package cache is `.conda-pkgs`. Set `CONDA=/path/to/conda`
if conda is not installed under `~/miniforge3`. All Python commands below use that
environment. Environment registration is disabled and plotting caches stay
under `build/`. No activation or changes to the base environment are needed.

`make cpu` produces `lib/libchem.so` and `bin/ts1_bench_cpu`; `make bench` also
builds `bin/ts1_scenario`. The library filename is also `.so` on macOS. Compiler and
flag changes automatically rebuild the C objects. `OPTFLAGS` defaults to
`-O3 -DNDEBUG`; the Makefile adds language, warning, include, OpenMP, and link
flags only. On macOS it selects the Xcode SDK with `xcrun` and uses Homebrew's
libomp for Apple clang. On other systems the compiler must support `-fopenmp`.

`make test-cpu` reports each of the 680 supplied case IDs individually and exits
nonzero on any failure. Cases cover native code through ctypes and standalone C
programs. It also runs 147 independent C/prototype comparisons, stiff analytical
and formal-order tests, failure tests, and scenario I/O/forcing/thread checks.
No MICM source, executable, or Python package is needed by this target.
To rerun a subset after building:

```sh
.conda/bin/python tests/run_cases.py --filter 'rosenbrock:*' --traceback
```

See [RESULTS.md](RESULTS.md) for results from both compilers and the TS1 gate.

## Run a chemistry scenario

```sh
make bench CC=clang
bin/ts1_scenario --method rodas3 --option rtol=0.003 --option atol=1e-18 --option h_max=200 --option group_size=32 ../scenarios/diurnal-10min/scenario.txt
```

The runner accepts `--threads N`, `--method NAME`, and repeated
`--option KEY=VALUE`. Without a method it uses the original three-stage
Rosenbrock and unchanged MICM defaults. New methods and every option are
documented in [DESIGN.md](DESIGN.md). All cores are used by default.

A scenario specifies its cells CSV, step duration, number of steps, a
piecewise-linear photolysis table, and a snapshot interval. Each external step
updates photolysis from the original cell parameters and holds it constant
through that solve. There is no warmup. The runner prints native `elapsed_ms`
and solver counters, writes `ts1_final.csv` and `ts1_snapshots.csv` in its working
directory, and returns nonzero on failure. Timing includes rate updates and
integration only. Use separate working directories to retain multiple runs.

## Reproduce the research

[SCANS.md](SCANS.md) describes the controlled experiments; `scans/*.json` stores
their configurations and `scans/*.csv` stores every native observation, including
failures and the supplied scorer's JSON output. The dense prototype and its
scope are described in [prototypes/README.md](prototypes/README.md).

```sh
# Generate local 10,000-cell inputs and independent 48-hour cases.
.conda/bin/python tools/prepare_experiments.py --large --stress

# Rebuild and check both complete compatibility gates, sequentially.
.conda/bin/python tools/validate.py --compiler clang
.conda/bin/python tools/validate.py --compiler gcc-15
make bench CC=clang

# Score every declared point on every disclosed 1,000-cell scenario.
.conda/bin/python tools/scan.py scans/operating-points.json --log runs/recheck-1000.csv --run-root runs/recheck-1000

# Primary-size native timings: median of five, with prefix accuracy explicitly labeled.
.conda/bin/python tools/scan.py scans/operating-points.json --log runs/recheck-10000.csv --run-root runs/recheck-10000 --scenario-root runs/large --reference-root ../scenarios --prefix-cells 1000 --repeat 5
```

The larger inputs reproduce the disclosed 1,000-cell prefix exactly. The
remaining cells have no supplied reference, so their accuracy is not claimed
from that prefix score. The native integrator checks every accepted state for
finiteness and nonnegativity. Do not run builds or other compute jobs alongside
timing experiments. Full reference/stress reproduction and the final measurement
sequence are in [SCANS.md](SCANS.md). Regenerate the committed report and figures
from the final CSVs with:

```sh
.conda/bin/python tools/report_results.py
.conda/bin/python tools/plot_scans.py
```

## Run the retained Tier 1 benchmark

```sh
make bench CC=clang
bin/ts1_bench_cpu --threads 1 100 100 30 ../reference/ts1_cells_100_seed1.csv
.conda/bin/python ../reference/compare_final.py ts1_final.csv ../reference/micm_ts1_final_100_seed1.csv
```

The command syntax is:

```text
bin/ts1_bench_cpu [--threads N] [--group-size G] CELLS STEPS DT CELLS_CSV
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

To reproduce a full comparison with a reference built using the same flags:

```sh
make cpu micm-ref CC=clang OPTFLAGS='-O3 -DNDEBUG'
mkdir -p runs/c runs/micm
.conda/bin/python ../reference/make_cells_csv.py 10000 2026 --output runs/cells.csv
(cd runs/c && ../../bin/ts1_bench_cpu 10000 100 30 ../cells.csv > timing.txt)
(cd runs/micm && ../../bin/micm_ts1_ref 10000 100 30 ../cells.csv > timing.txt)
.conda/bin/python ../reference/compare_final.py runs/c/ts1_final.csv runs/micm/ts1_final.csv
```

`make micm-ref` and `make bin/micm_scenario` compile the unchanged supplied C++
reference runners using read-only inputs under `../reference`, `../scenarios`
and `../micm`. Normal C builds are independent of those directories. The
committed configuration is enough to rebuild TS1.

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

The retained translation command for TS1 is:

```sh
.conda/bin/python tools/translate_ts1.py ../micm/test/benchmark/ts1_mechanism.hpp configs/ts1.json
```

Public APIs are in `include/chem.h` and `include/chem_array.h`. A typical caller
creates a model, adds species/parameters/reactions, finalizes it, creates one
state per cell group, sets conditions and concentrations, calls
`ch_update_rates`, then calls `ch_solve` or creates a `ch_integrator` for the new
ODE methods. Free an integrator before its state, and states before their model.
`tools/tutorials.c` contains complete native examples. [DESIGN.md](DESIGN.md)
documents ownership, numerical conventions, and extension points.

The numerical methods and adapted test fixtures retain MICM's Apache-2.0
attribution. See `LICENSE` and `tests/reference/README.md`.
