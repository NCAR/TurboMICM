# CPU chemistry port (in progress)

Create the local tooling environment with `make env`. `make cpu` builds the C11
mass-action library. `make check CC=clang` or `make check CC=gcc-15` builds and
runs the current focused tests; the tests read tutorial configuration input from
the adjacent `micm` checkout. `OPTFLAGS` can be supplied to make.

A preliminary single-cell three-stage Rosenbrock solver is available through
`include/rosenbrock.h`. The required full test suite, remaining rate types, additional
solver features, and TS1 benchmark remain to be implemented. Run `make clean`
before switching compilers.

Regenerate the TS1 configuration from the reference source with
`.conda/bin/python tools/import_ts1.py ../micm/test/benchmark/ts1_mechanism.hpp configs/ts1`.
The current checks compare the translation's initial values with the reference
defaults CSV and exercise generated TS1 rate and forcing evaluation.

Build the preliminary TS1 driver with `make bench`, then run
`bin/ts1_bench_cpu [--threads N] CELLS STEPS DT CELLS_CSV`. It writes
`ts1_final.csv`. `make check-ts1` compares 100 cells and 100 steps with the
supplied reference, and `make micm-ref` builds the independent timing reference.
The current solver uses independent dense factorizations per cell; it is not yet
optimized for the requested benchmark scale.
