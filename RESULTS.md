# Results

Measured on the experiment's machine: Apple M5 (4 performance and 6
efficiency cores), 16 GB, macOS 27.0.1, Apple clang 21.0.0 and GCC 15.2.0,
`OPTFLAGS=-O3 -DNDEBUG` for the port and for MICM. Times are from the
programs' own `elapsed_ms`; the experimenter's timing session is the ranking.

## Test cases

(to be filled from `make test-cpu` with each compiler)

| Compiler | PASS | FAIL | N/A |
|---|---:|---:|---:|
| Apple clang 21 | | | |
| GCC 15 | | | |

## TS1 gate

(to be filled from `tools/check_benchmark.py`)

## TS1 timing

(to be filled: 10 threads and 1 thread, MICM's default and fastest
configurations, generic kernels, exact mode, group sizes)
