# Design

`include/chem.h` is the public C API. `src/rates.c` implements the rate formulas
and preparation of derived parameters. `src/mechanism.c` owns species and
reaction metadata and performs stable rate-type sorting. `src/state.c` manages
cell-contiguous concentrations, parameters and rates, and evaluates mass-action
forcing and Jacobians. The current dense Jacobian evaluator is also a reference
path for validating the planned sparse evaluator.

Configuration files pass through `tools/configure.py`: included CAMP records
are collected, names and numeric values are validated, branched reactions are
expanded, units converted, and a C builder emitted. Generated C is compiled by
the Makefile with the selected compiler and optimization flags. The runtime
requires no configuration parser or Python. TS1's data is retained in
`mechanisms/ts1.json`; `tools/translate_ts1.py` reproduces it from the supplied
read-only input.

To add a rate type, add a tagged parameter member in `chem.h`, supply defaults,
preparation and evaluation in `rates.c`, and add a configuration type mapping
and independent reference tests. Callbacks have an explicit context pointer
whose lifetime remains the caller's responsibility.

Solver/controller/linear-solver interfaces and tuning parameters are under
construction. They must be completed and documented before the full task is
considered finished.
