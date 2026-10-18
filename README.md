# Portable chemical kinetics in C

ISO C11 chemistry library with Rosenbrock, backward Euler and Euler extrapolation,
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
OpenMP flags. Apple clang uses `/opt/homebrew/opt/libomp`; GCC uses `-fopenmp`.
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

## Tier 3 numerical research

`make bench` also builds the scenario driver. For example:

```sh
bin/ts1_scenario --method ros4 --option rtol=1e-4 --option atol=1e-18 ../scenarios/diurnal-10min/scenario.txt
.conda/bin/python ../scenarios/score.py ts1_snapshots.csv ../scenarios/diurnal-10min/reference_snapshots.csv --level 3
make test-methods
```

The default method is `micm3`, retaining MICM's original parameters and
tolerances. New methods are additions selected at runtime. See
[DESIGN.md](DESIGN.md) for method names, all options and how to add one;
[ACCURACY.md](ACCURACY.md) for the scientific error budget;
[REFERENCES.md](REFERENCES.md) for the literature used. C/prototype agreement
and scenario parser checks are in `make test-methods`; the original Tier 1 gate
remains `make test-cpu`. `bin/micm_scenario`, if built from the supplied read-only
source, enables an additional independent scenario-driver parity check.

Research tooling runs only under `.conda/bin/python`. `tools/scan.py` invokes
the C driver for timing/counters and the supplied scorer for accuracy; it writes
one CSV row per attempted configuration/scenario. Example:

```sh
.conda/bin/python tools/scan.py --suite family --output scans/family.csv
```

The five evaluated configurations and the model recommendation are in
[OPERATING_POINTS.md](OPERATING_POINTS.md), backed by
[operating_points.json](operating_points.json). The candidate methods include
RODAS3.1, PI/H211b control and Euler extrapolation; the selected points use
current Jacobians and exact retry-forcing reuse. [SCANS.md](SCANS.md) summarizes
the controlled experiments and links the raw CSV records. Scientific prototype
and controller-roundoff details are in [prototypes/README.md](prototypes/README.md).

To repeat the declared accuracy checks on the independent 48-hour cases:

```sh
# First build the supplied reference driver using the command in scenarios/README.md.
.conda/bin/python tools/stress.py --references --point recommended
.conda/bin/python tools/stress.py --point L2
.conda/bin/python tools/stress.py --point L2_margin
.conda/bin/python tools/stress.py --point L3
.conda/bin/python tools/stress.py --point L4
```

`tools/scale_scan.py` generates the 10000-cell version of diurnal-10min and
measures grouping/scheduling. Its accuracy score uses only the disclosed
1000-cell prefix; its finite/nonnegative gate checks all 10000 cells.
`tools/operating_points.py --repeat 5` reproduces the final C timing protocol.
`tools/micm_scan.py` measures the independently compiled original C++ baseline.
All generated large CSVs stay in ignored `results/` directories; concise
accuracy/work/time records stay in versioned `scans/`.

The report and work–precision figures are reproducible with
`tools/report_results.py` and `tools/plot_scans.py`. Use `.conda/bin/python` for
both. Run timing studies into a fresh `--output` CSV when retaining the committed
reference records; the scan tools append every run instead of discarding history.
`tools/scale_scan.py --points` checks all five declared points at 10000 cells.
