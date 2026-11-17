# CPU chemistry port

`make env` creates a local conda prefix from `environment.yml`. `make cpu`
builds the C library and TS1 benchmark, and `make test-cpu` reports the status
of every reference case (676 pass and four C++ template-only cases are N/A).
Select a compiler with `CC` and optimization flags
with `OPTFLAGS` (default `-O3 -DNDEBUG`). Run `make clean` when switching
compilers. `make test-reader` runs additional reader, Jacobian, and solver
checks. `make micm-ref` builds the reference benchmark with the same `OPTFLAGS`.

Run `bin/ts1_bench_cpu [--threads N] CELLS STEPS DT CELLS_CSV`. It performs a
warm-up step and the requested timed steps and writes `ts1_final.csv`.

`tools/import_ts1.py` translates the reference TS1 builder into the checked-in
configuration. `tools/generate.py CONFIG OUTPUT` reads a tutorial-style JSON or
YAML configuration and writes a C table. The runtime is C11 and does not load
Python or MICM.

The C interfaces for chemistry, matrix storage, state, constraints, and solver
construction are in `include/`. Tutorial configurations and reference outputs
used by the tests are included, so ordinary builds and tests do not require a
neighboring source tree. `make micm-ref` and the TS1 import script are optional
comparison and translation tools that do require the supplied reference tree.
