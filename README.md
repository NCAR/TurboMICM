# Portable chemical kinetics in C

ISO C11 chemistry library with Rosenbrock and backward Euler integration,
analytic sparse Jacobians, differential-algebraic constraints, external models,
and dense/sparse matrix utilities. Runtime dependencies are the C standard
library, libm and OpenMP. Python is used only for configuration generation and
tooling. See [DESIGN.md](DESIGN.md) for the public interfaces and extension guide,
[RESULTS.md](RESULTS.md) for validation, and [DECISIONS.md](DECISIONS.md) for choices.

## Build and test

```sh
make env
make cpu
make test-cpu
make CC=gcc-15 test-cpu cpu
```

`make env` creates `.conda` from the pinned, cross-platform `environment.yml` and
keeps its package cache in `.conda-pkgs`. Set `CONDA=/path/to/conda` if necessary.
Generation targets create the environment automatically when needed. Nothing
is installed into the base environment.

`CC` defaults to clang and `OPTFLAGS` to `-O3 -DNDEBUG`. Switching either rebuilds
all C objects. On macOS Make sets `SDKROOT` and selects the requested compiler's
OpenMP flags. Apple clang uses `OMP_PREFIX=/Users/dfillmor/homebrew/opt/libomp`; GCC uses `-fopenmp`.
On Linux the compiler uses `-fopenmp`. `OMPFLAGS` and `OMPLIBS` can be overridden
for another OpenMP installation. Make adds language/warning/include/OpenMP/link
flags only. `make clean` removes build products and generated C.

The build produces `build/libchem.a`, `bin/test_cpu`, and
`bin/ts1_bench_cpu`. `make test-cpu` first validates the configuration reader,
then prints one result for each of the exact 680 IDs in `tests/cases.txt`.
Assertions remain active under `-DNDEBUG`. The five N/A cases contain exclusively
C++ compile-time concept/category assertions; numerical and runtime cases pass.
The three tutorial configurations are also checked against independently built
C descriptors, including rates, forcing, Jacobians and trajectories.

## TS1 benchmark

```sh
make bench
bin/ts1_bench_cpu --threads 4 100 100 30 ../reference/ts1_cells_100_seed1.csv
.conda/bin/python ../reference/compare_final.py ts1_final.csv ../reference/micm_ts1_final_100_seed1.csv
```

The command accepts `[--threads N] [--group-size N] CELLS STEPS DT CELLS_CSV`.
It runs one untimed warm-up, times the requested steps, prints `elapsed_ms` and
summed solver statistics, and writes `ts1_final.csv` in the working directory.
The CSV begins with `temperature,pressure` followed by named species. Input
columns can name temperature, pressure, species or custom rate parameters.
Omitted columns retain TS1's initial defaults. The first `CELLS` rows are used.
Species output order follows the model's permutation; comparison is by name.

By default OpenMP uses all available threads and groups of up to 32 cells.
The default groups shrink for small inputs to keep available threads occupied.
`--threads` limits threads, `--group-size` or `CK_GROUP_SIZE` changes group size,
and `OMP_SCHEDULE` changes scheduling (default static). A final partial group
contains exactly its remaining cells. `make GROUP_SIZE=256 cpu` changes the
compiled default. All sizes must be positive; `STEPS=0` performs only warm-up.

The six photolysis defaults called out in the task remain 1/s in the translated
configuration. The supplied perturbed inputs set them to zero. To reproduce a
full-size reference comparison, use the local Python environment:

```sh
mkdir -p results
.conda/bin/python ../reference/make_cells_csv.py 10000 2026 --output results/cells-10000-seed2026.csv
make micm-ref
# Run the two programs in separate working directories to preserve both CSVs.
```

`make micm-ref` alone builds the C++ reference from `../reference` and `../micm`,
with clang++ and the same `OPTFLAGS`. Library, benchmark and ordinary test builds
use only this repository's C/configuration files. Reference tests must be run
from a copy of their build directory inside this repository.

## Mechanism input and tutorials

Read `include/chem.h` for descriptors, state management and integration;
`include/chem_matrix.h` for matrix/view utilities; and `include/chem_verify.h`
for finite-difference validation. Construct `ck_species`, `ck_reaction`, and
`ck_mechanism` arrays in C, or generate them from CAMP YAML/JSON:

```sh
.conda/bin/python tools/configure.py mechanisms/tutorial/robertson/config.json generated/example.c --symbol example
make bin/tutorial_cpu
bin/tutorial_cpu README_example
```

`tools/import_ts1.py` retains the one-time C++-to-configuration translation.
Regenerate with:

```sh
.conda/bin/python tools/import_ts1.py ../micm/test/benchmark/ts1_mechanism.hpp mechanisms/ts1.json
```

Normal builds read the committed JSON and run the same configuration reader as
for the tutorials. No C++ translation is needed at build or runtime. Generated
C contains immutable descriptors and is compiled by Make with `CC`/`OPTFLAGS`.
Numerical formulas and test data derived from MICM retain Apache-2.0 attribution
in [NOTICE](NOTICE) and [LICENSE](LICENSE).

## Generated mechanism solvers

```sh
make -j6 bench test-specialized
bin/ts1_bench_specialized --threads 18 100 100 30 ../reference/ts1_cells_100_seed1.csv
bin/ts1_scenario_specialized --threads 18 --group-size 32 --rtol 1e-6 --atol 1e-15 ../scenarios/diurnal-10min/scenario.txt
make CC=gcc-15 -j6 bench test-cpu test-specialized
```

`tools/specialize.py` reads the same CAMP JSON/YAML configurations as the generic
reader and emits C11 chemistry, rate and sparse linear algebra kernels. The
three-stage Rosenbrock method retains MICM's coefficients, error norm and step
controller. `make bench` builds both Tier 1 and specialized benchmarks and the
specialized scenario driver. Their runtime requires C, libm and OpenMP only.

`make test-specialized` generates TS1, the three Tier 1 configuration-file
mechanisms and fixtures covering every supported configurable rate and a zero-order mechanism. An additional build tests combined forcing/Jacobian
emission, row-oriented LU updates, and odd preferred/tile widths together. It compares
rates, forcing, Jacobians, factorizations, direct solves, trajectories and solver
statistics against the generic library, using a relative bound of `1e-10`.
Widths 1, 3, 32, 37 and 64 exercise full and partial groups; tolerances,
environmental inputs and custom parameters change during testing. Python tests
also exercise regeneration after a reaction/parameter edit.

The benchmark has the same arguments, warm-up, timed region, output and default
tolerances as `ts1_bench_cpu`. Scenario files follow `../scenarios/README.md`:
paths resolve relative to the scenario file, each step scales every `j*` custom
parameter from its original CSV value, and there is no warm-up. Scenario timing
includes interpolation, rate updates and solves, excluding input/output and
snapshot copies. `ts1_final.csv` and requested `ts1_snapshots.csv` use 17 digits.
Scenario group size defaults to 32 and never shrinks according to thread count.
`--rtol` and `--atol` override relative and uniform absolute tolerances; omitted
options retain MICM defaults. Day-night verification uses `1e-6` and `1e-15`.

To generate a different mechanism:

```sh
.conda/bin/python tools/specialize.py path/to/config.yaml generated/custom.c --symbol custom
```

Compile the output with the same include paths as Make, and link
`build/libchem.a`, OpenMP and libm. Include `chem_specialized.h`, declare
`extern const ck_specialized custom`, construct a model from `custom.mechanism`,
and use `custom.update_rates(state)` and `custom.solve(state, dt, &stats)`.
These kernels require states from that descriptor's default ordered model.
All generated C is disposable and is regenerated from configuration dependencies.
See [SPECIALIZATION.md](SPECIALIZATION.md) for measured experiments and
[DESIGN.md](DESIGN.md) for tuning parameters.

The disclosed scenario gate can be reproduced with
`tools/check_scenarios.sh clang 18` after a clang build, or
`tools/check_scenarios.sh gcc 18` after a GCC build. Each run stays under
`results/<label>/`, checks final states and snapshots against MICM, and reports
snapshot agreement with the high-accuracy reference.
