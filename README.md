# CPU chemistry port (work in progress)

`make env` creates the local conda prefix from `environment.yml`. `make cpu`
compiles the C library and generated TS1 tables; `make test-reader` runs the
configuration and forcing/Jacobian checks. Choose the compiler using `CC` and
optimization flags using `OPTFLAGS`. Run `make clean` when switching compilers.

`tools/import_ts1.py` translates the reference TS1 builder into the checked-in
configuration. `tools/generate.py CONFIG OUTPUT` reads a tutorial-style JSON or
YAML configuration and writes a C table. The runtime is C11 and does not load
Python or MICM.

The numerical solver, benchmark, and full case-by-case compatibility tests are
not yet implemented. The existing checks are not a substitute for those gates.
