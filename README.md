# Portable chemistry library

C11 chemistry integration with build-time YAML/JSON configuration reading,
OpenMP across independent groups of cells, and SIMD across contiguous cells.

The port is in progress. The full reference case gate currently reports pending
cases as FAIL; see RESULTS.md for verified coverage and outstanding work.

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

`make test-cpu` prints every supplied case ID and its status. A pending case is a
failure, never an automatic pass. Supplemental checks verify configuration
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
