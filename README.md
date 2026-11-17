# CPU chemistry port (in progress)

Create the local tooling environment with `make env`. `make cpu` builds the C11
mass-action library. `make check CC=clang` or `make check CC=gcc-15` builds and
runs the current focused tests; the tests read tutorial configuration input from
the adjacent `micm` checkout. `OPTFLAGS` can be supplied to make.

A preliminary single-cell three-stage Rosenbrock solver is available through
`include/rosenbrock.h`. The required full test suite, rate constants, additional
solver features, and TS1 benchmark remain to be implemented. Run `make clean`
before switching compilers.
