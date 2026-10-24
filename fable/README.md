# A C port of MICM's CPU chemistry solver

`libtmc.a` (header `include/tmc.h`) reproduces MICM's CPU solver in ISO C11
with OpenMP: species, reactions, and all of MICM's rate-constant types; a
solver builder; Rosenbrock and backward Euler integrators; sparse LU; and
constraints and external models. Python tooling reads MICM's configuration
format and generates, for any mechanism, the C tables and kernels the library
builds from. See DESIGN.md for the architecture, DECISIONS.md for the choices,
and RESULTS.md for the measurements.

## Requirements

- A C11 compiler with OpenMP: Apple clang with libomp, or GCC. On macOS the
  Makefile sets `SDKROOT` and uses libomp from `/opt/homebrew/opt/libomp`.
- GNU Make, and conda for the build-time Python environment.
- MICM's source tree at `../micm` and the experiment's `../reference`
  directory (or `make MICM=... REFERENCE=...`). The benchmark and the library
  do not use MICM; the TS1 translation reads MICM's `ts1_mechanism.hpp`, and
  the tests read MICM's test data.

## Build

```
make env                       # the conda environment in .conda (once)
make cpu                       # library, generated TS1 kernels, bin/ts1_bench_cpu
make cpu CC=gcc-15             # the same with GCC
make bench                     # alias of cpu
make micm-ref                  # MICM's reference driver, bin/micm_ts1_ref, with clang++
```

`OPTFLAGS` (default `-O3 -DNDEBUG`) and `CC` are the only compiler choices;
the Makefile adds language, warning, include, OpenMP, and link flags.
Build-time tuning parameters: `WIDTH` (lanes per block, default 2),
`GROUP_SIZE` (cells per solver state in the benchmark, default 8),
`KERNEL_MODE` (`fast` or `exact`), `TS1_ORDERING` (`minimum-fill`, or
`markowitz` for MICM's species order, the default in exact mode),
`KERNEL_CHUNK`, and `GEN_FLAGS` passed to the generator. The generated files live in `generated/` and are not
committed.

## Test

```
make test-cpu                  # every case of ../reference/micm_cases.txt: <id> PASS|FAIL|N/A <reason>
make test-cpu TEST_FILTER=rosenbrock
```

The suite is C (`tests/`), one function per MICM case, run in the order of
the reference list. tests/README.md describes the conventions.

## Run the TS1 benchmark

```
bin/ts1_bench_cpu [--threads N] [--group-size G] [--schedule static|dynamic|guided] [--chunk K] [--generic] CELLS STEPS DT CELLS_CSV
```

For example, the ranking run with the experiment's input:

```
.conda/bin/python ../reference/make_cells_csv.py 10000 2026 --output cells.csv
bin/ts1_bench_cpu 10000 100 30 cells.csv
.conda/bin/python ../reference/compare_final.py ts1_final.csv micm_ts1_final.csv
```

The program prints `elapsed_ms=` and the solver statistics summed over the
timed steps, and writes `ts1_final.csv` in the working directory.
`--generic` runs the library's generic kernels instead of the generated TS1
kernels; `--threads 1` and `--group-size` change the threading and the
grouping (32 cells per group by default),
and all of them must agree with MICM within the gate.

## Use the library

```c
#include "tmc.h"
#include "ts1.h"   /* generated */

tmc_system system = ts1_system();
tmc_reaction *reactions;
int n = ts1_reactions(&reactions);
tmc_builder *b = tmc_builder_create_rosenbrock(tmc_rosenbrock_three_stage());
tmc_builder_set_system(b, &system);
tmc_builder_set_reactions(b, reactions, n);
tmc_builder_set_kernels(b, ts1_kernels());      /* optional: the generated kernels */
tmc_error error;
tmc_solver *solver = tmc_builder_build(b, &error);
tmc_state *state = tmc_solver_create_state(solver, 8);   /* 8 cells, one step size */
ts1_initial_state(state);                               /* the configuration's values in every cell */
/* set conditions (tmc_state_conditions), concentrations (tmc_state_set, or
   tmc_state_set_concentrations per species), parameters ... */
tmc_solver_update_state_parameters(solver, state);       /* rate constants */
tmc_result result = tmc_solver_solve(solver, 30.0, state);
```

A mechanism can also be built in code with `tmc_species_create`,
`tmc_phase_add_species`, `tmc_reaction_create`, `tmc_reaction_add_reactant`,
`tmc_rate_arrhenius` and the other constructors, as the tests do.

## Other mechanisms

```
.conda/bin/python tools/generate.py path/to/config.json generated/mine --name mine --width 2
```

reads a configuration (YAML or JSON in the format of `configs/examples/`)
and writes `generated/mine.c` and `mine.h`; compile `mine.c` with the library
and call `mine_system()`, `mine_reactions()`, and `mine_kernels()`.
`tools/mechanism.py` prints a mechanism's size, sparsity, and fingerprint.

## Citation

To cite TurboMICM, use
[10.5281/zenodo.23618864](https://doi.org/10.5281/zenodo.23618864).
See the [project citation](../README.md#citation) for release details.
