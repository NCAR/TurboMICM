# CPU chemistry library

This repository translates the CPU chemistry mechanisms, rate constants, and numerical solvers to ISO C11. Python reads YAML/JSON configurations and generates C when building; the library and benchmark have no Python dependency.

Build the private tooling environment with `make env`, then build with `make cpu`. Compilers and optimization flags are configurable, for example `make cpu CC=gcc-15 OPTFLAGS='-O3 -DNDEBUG'`. OpenMP is required. On macOS the Makefile sets the SDK and finds Homebrew libomp for Apple clang.

Run the benchmark with:

```
bin/ts1_bench_cpu --threads 10 100 100 30 ../reference/ts1_cells_100_seed1.csv
```

It runs one warm-up step, times the requested steps, prints elapsed time and summed solver statistics, and writes `ts1_final.csv`. Omit `--threads` to use the OpenMP default. `--group N` changes the number of cells sharing an adaptive solver controller.

`make micm-ref` builds the independent reference executable using the supplied read-only headers. It is separate from the C library and benchmark.

The test suite and broader solver functionality are under construction; progress and measurements are recorded in `RESULTS.md`.
