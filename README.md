# Portable chemistry C port

This repository contains a generic ISO C11 chemistry mechanism library, Python
configuration readers and code generation, a CPU Rosenbrock integrator, and a
TS1 benchmark driver. Python is used at build time only; the executable is C and
links only the C library, libm, and OpenMP.

## Build

The compiler defaults to the system C compiler. Select a compiler with `CC` and
optimization flags with `OPTFLAGS` (default `-O3 -DNDEBUG`):

```sh
make cpu CC=clang
make clean
make cpu CC=gcc-15
make test-cpu CC=clang
```

`make test-cpu` runs the C unit checks and reports each case in the supplied
MICM inventory. Implemented cases print `PASS`; cases without a meaningful C
equivalent print `N/A`; unported cases print `FAIL`. The target exits nonzero
until every non-`N/A` case is implemented.

`make env` creates the pinned Python tooling environment under `.conda/`, with
its package cache under `.conda-pkgs/`. Use `.conda/bin/python` for YAML input.
JSON reading uses Python's standard library.

## Configuration and code generation

`tools/read_config.py` resolves CAMP `camp-files` entries relative to the input
configuration and writes normalized JSON. For example:

```sh
.conda/bin/python tools/read_config.py path/to/config.yaml build/mechanism.json
```

`tools/import_ts1.py` translates the permitted TS1 builder source to
`mechanisms/ts1.json`; `tools/generate_c.py` emits immutable C arrays at
`generated/ts1.c`. The generated representation is mechanism-generic and
contains species metadata, stoichiometry, rate parameters, and initial defaults.

## Run the TS1 benchmark

```sh
bin/ts1_bench_cpu [--threads N] CELLS STEPS DT CELLS_CSV
```

The CSV has a header of `temperature`, `pressure`, species names, and custom
rate parameter labels. Missing columns use TS1 initial defaults. The executable
uses the first `CELLS` rows, takes one untimed warm-up step, then runs `STEPS`
timed steps of length `DT` seconds. It prints elapsed time and solver counters
and writes the final state to `ts1_final.csv`.

To build the permitted comparison executable, run `make micm-ref`. The supplied
comparison utility checks every species and cell against the reference:

```sh
python3 ../reference/compare_final.py ts1_final.csv path/to/reference.csv
```

## C API

`include/chemistry.h` defines the generic mechanism records, non-owning state view,
strided matrix views, rate evaluation, forcing and Jacobian routines, linear and
equilibrium constraint evaluation helpers, version API, and solver interface.
Matrix views support checked element access, row-major vector conversion,
row-to-vector copies, element traversal, fill, AXPY, and min/max/sum reductions.
The constraint helpers
build algebraic masks and Jacobian patterns and add residual and negative-Jacobian
terms to caller-provided arrays. Masked forcing and Jacobian helpers accumulate
reaction contributions while leaving algebraic rows unchanged. The optional
constrained Rosenbrock entry point
projects algebraic variables with Newton iterations before integration and after
each accepted step, and reports the initial projection iteration count in
`chem_solver_stats.constraint_init_iterations`. The mechanism evaluator covers
Arrhenius, photolysis, user-defined, Troe, surface, ternary activation (including
JPL19), tunneling, branched, Taylor-series, reversible, and callback rates.
Constrained calls support scalar and sequential per-cell grouped operation; each
grouped cell has its own adaptive controller.
`chem_rosenbrock_builtin_parameters` supplies MICM's two-, three-, four-, and
DAE parameter sets; pass a parameter structure to
`chem_rosenbrock_solve_with_parameters` to use a built-in or custom method. The
existing `chem_rosenbrock3_*` functions use the default three-stage set.

See `DESIGN.md` for the architecture and `RESULTS.md` for validation and current
limitations.
