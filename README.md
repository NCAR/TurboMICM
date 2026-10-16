# Portable chemistry library

Work in progress: the TS1 integration gate and the initial rate-constant tests
pass, but the complete CPU test port is not yet finished. See WORK.md and
RESULTS.md for the current evidence and remaining scope.

The runtime is C11 with the standard C library, libm, and OpenMP. Configuration
reading and C generation use a private Python environment at build time.

```
make env
make cpu CC=clang
make test-cpu CC=clang
make cpu CC=gcc-15
make test-cpu CC=gcc-15
make bench
bin/ts1_bench_cpu --threads 4 100 100 30 ../reference/ts1_cells_100_seed1.csv
```

`make test-cpu` currently fails for cases that have not yet been ported. To run
only the implemented rate tests, use `make test-cpu TEST_FILTER=_rate_constant:`.
Every selected case prints its original ID, PASS/FAIL/N/A, and a reason.

`CC` selects the C compiler and `OPTFLAGS` defaults to `-O3 -DNDEBUG`. Compiler
or flag changes invalidate the build objects. The Makefile chooses the OpenMP
flags for Clang on macOS, Clang elsewhere, or GCC; on macOS it sets SDKROOT using
`xcrun`. The conda environment and package cache live inside this checkout.
`CONDA` may override the conda executable path. Dependencies are version-pinned
without platform-specific build strings.

The benchmark accepts `--threads N`, `--group-size G`, and `--generic` before
`CELLS STEPS DT CELLS_CSV`. It performs one untimed warm-up step, times the
requested subsequent steps, reports their accumulated statistics, and writes
`ts1_final.csv`. Omitting `--threads` uses all processors OpenMP reports.
`--generic` selects the table-driven chemistry evaluator for verification.

```
make micm-ref
bin/micm_ts1_ref --threads 1 100 100 30 ../reference/ts1_cells_100_seed1.csv
.conda/bin/python ../reference/compare_final.py YOUR_FINAL.csv REFERENCE_FINAL.csv
```

The reference binary also writes `ts1_final.csv`; run it in a separate directory
or preserve the C output first. Reference paths can be overridden with `MICM`
and `REFERENCE` Make variables. Normal library/benchmark builds include no
reference headers or code and have no Python runtime dependency.

To regenerate the retained TS1 configuration from the supplied data:

```
.conda/bin/python tools/translate_ts1.py ../micm/test/benchmark/ts1_mechanism.hpp configs/ts1.json
```

See DESIGN.md for the library interfaces and DECISIONS.md for format extensions.
