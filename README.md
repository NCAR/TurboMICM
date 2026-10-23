# Portable chemistry library

C11 chemistry integration with build-time YAML/JSON configuration reading,
OpenMP across independent groups of cells, and SIMD across contiguous cells.

All 675 applicable reference cases pass with Clang and GCC; five cases concern
C++ compile-time mechanics and have documented N/A reasons. See RESULTS.md for
the TS1 accuracy checks and local timing observations.

## Build and test

```sh
make env
make cpu CC=clang
make test-cpu CC=clang
make test-cpu CC=gcc-15
```

`make env` creates the prefix environment `.conda` using `environment.yml`, with
its package cache in `.conda-pkgs`. Python is used for configuration generation
and test orchestration. The library and benchmark run without Python.

`CC` selects the compiler; `OPTFLAGS` defaults to `-O3 -DNDEBUG`. Changing either
rebuilds the C outputs. macOS builds set SDKROOT automatically. Clang on macOS
uses Homebrew libomp; GCC and Linux Clang use `-fopenmp`.

`build/libchem.so` exposes the public API in `include/chem.h`. Descriptors can be
assembled in C or generated from configurations. Mechanisms borrow descriptor
arrays, which must remain alive until the mechanism and its states are freed.
A mechanism can be shared by independent states on different threads.
`chem_builder.h` assembles gas species, external-model variables and parameters,
and constraints before constructing a mechanism. Builder-produced descriptor
tables are owned by the mechanism; callback contexts and equation arrays remain
borrowed.

`chem_system.h` provides owned species properties, phase metadata, and surface
rate construction. `chem_jacobian.h` provides finite-difference and sparsity
diagnostics. Named state updates, independent state copies, explicit ownership
transfer, and separate or combined LU factors are available through `chem.h`.

`include/chem_array.h` supplies typed dense/grouped arrays, strided views,
reusable functions, and reductions. `include/chem_sparse.h` supplies CSR/CSC
matrices with the same view operations across blocks. Sparse numeric storage is
a borrowed array of blocks by structural nonzeros, so callbacks can combine
sparse elements, dense columns, and vectors through one interface.

`make test-cpu` prints every supplied case ID and its status, including explicit
reasons for five N/A cases. Supplemental checks verify configuration
reading, configuration-versus-code equivalence, and a stiff TS1 regression.

## Benchmark

```sh
make bench CC=clang
bin/ts1_bench_cpu [--threads N] [--group-size N] CELLS STEPS DT CELLS_CSV
```

The executable performs one untimed warm-up, times STEPS solves of DT seconds,
prints elapsed milliseconds and aggregate solver statistics, and writes
`ts1_final.csv` in the current directory. The CSV input accepts temperature,
pressure, species names, and custom parameter labels. Omitted columns retain
the configuration's InitState defaults. Unknown or duplicate columns and invalid
numbers produce an error. The supplied reference inputs zero the six photolysis
rates whose InitState defaults are unphysical at surface conditions.

For a local reference comparison:

```sh
make micm-ref
mkdir -p build/reference-run
# Run bin/micm_ts1_ref from build/reference-run using absolute input paths.
.conda/bin/python ../reference/compare_final.py ts1_final.csv REFERENCE_FINAL.csv
```

Only `make micm-ref` compiles against the reference C++ source. To run the original
reference tests, copy `../micm/build` locally and run executables from that copy.
Do not invoke ctest in the copy: its paths still point to the original build.

## Configuration

```sh
.conda/bin/python tools/configure.py CONFIG OUTPUT.c --symbol mechanism_data
.conda/bin/python tools/translate_ts1.py \
  ../micm/test/benchmark/ts1_mechanism.hpp configs/ts1.json
```

The reader follows `camp-files`, collects `CHEM_SPEC` and `MECHANISM` entries,
expands branched reactions, converts rate units, and emits immutable C data.
`--no-reorder` disables Markowitz species ordering. See DECISIONS.md for format
extensions and units. The committed TS1 JSON is translated data, not a mechanism
implemented specially in the C solver.

## Specialized build and day-night scenarios

```sh
make bench CC=clang
make test-specialized CC=clang
bin/ts1_bench_specialized --threads 18 10000 100 30 CELLS_CSV
bin/ts1_scenario_specialized --threads 18 --group-size 32 \
  --rtol 1e-6 --atol 1e-15 ../scenarios/diurnal-10min/scenario.txt
```

The specialized benchmark has the generic benchmark's warm-up, timed region,
statistics, and final CSV format. The scenario runner starts directly from the
input, applies each step's photolysis factor, and writes `ts1_final.csv` and
`ts1_snapshots.csv` when a snapshot interval is requested. It accepts arbitrary
relative scenario paths and nonuniform photolysis tables with constant endpoints.
Run each executable in a separate output directory when retaining comparisons.

```sh
.conda/bin/python tools/specialize.py CONFIG OUTPUT.c --prefix mechanism
```

The output defines `mechanism_data`, `mechanism_update_rates`, and
`mechanism_solve`. Compile it with the generic library and `-Iinclude -Isrc`, using
the same C11/OpenMP flags as the Makefile. Create states with `mechanism_data`.
`make` automatically regenerates TS1 and the configuration-file test mechanisms;
generated C and objects live under ignored `build/` and are never hand-edited.
See DESIGN.md for tuning controls and SPECIALIZATION.md for measured alternatives.

`make test-scenarios` runs all five disclosed scenarios and scores their snapshots
against MICM at level 10, then against the high-accuracy references for scale.
For retained accuracy and timing evidence after building:

```sh
make micm-ref
.conda/bin/python tools/evaluate.py --scenarios --benchmarks \
  --repeats 3 --directory validation/local
```

This runs the programs sequentially, retains logs and CSVs, checks both benchmark
outputs against MICM, and records results in `results.json`. Python orchestrates
separate C executables; it is never involved in a solver or timed region.
