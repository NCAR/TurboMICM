# Portable chemistry C port

This repository is being developed as an ISO C11 CPU port of MICM chemistry
functionality. The current commit contains a generic species/reaction data model,
mass-action forcing for the implemented rate types, configuration readers, and a
TS1 mechanism translation. The ODE integrator, command-line benchmark, and full
MICM case port are not implemented yet; `bin/ts1_bench_cpu` currently reports
that it is unavailable.

## Build

The compiler defaults to the system C compiler. Set `CC` to `clang` or `gcc-15`
to select a compiler, and set `OPTFLAGS` to change optimization flags. The
Makefile only adds the C11 language, warning, include, OpenMP, and link flags.

```sh
make cpu CC=clang
make clean
make cpu CC=gcc-15
make test-cpu CC=clang
```

`make test-cpu` checks the implemented C routines and emits one status line for
each case in the supplied MICM case inventory. At this stage most ports are
reported as failures, and the target exits nonzero until they are implemented.

`make env` creates the pinned Python configuration-tooling environment at
`.conda/`; its package cache also stays in this checkout. Use `.conda/bin/python`
to read YAML inputs. JSON reading needs only Python's standard library.

## Configuration tooling

`tools/read_config.py` resolves CAMP `camp-files` entries relative to the
configuration and emits a normalized JSON document. For example:

```sh
.conda/bin/python tools/read_config.py path/to/config.yaml build/mechanism.json
```

`tools/import_ts1.py` translates the allowed TS1 MICM builder source into
`mechanisms/ts1.json`. It checks the expected 210 species and 547 reactions.
`tools/generate_c.py` converts that normalized JSON to `generated/ts1.c`, which
is compiled into the CPU executable. The runtime benchmark path remains to be
completed.

## C API

`include/chemistry.h` declares mechanism, species, reaction, and stoichiometric
term records. `src/mechanism.c` implements Arrhenius, Troe, ternary activation,
tunneling, branched, surface, and user-supplied rates, along with mass-action
forcing. This API is an early foundation and does not yet implement the
requested solver behavior.

## Validation status

Both requested compiler commands currently build the C foundation. Configuration
reader checks have been run on the JSON and YAML tutorial mechanisms and the TS1
translation counts have been checked. This is not a MICM equivalence result;
there is no integrator or final-state benchmark yet. See `RESULTS.md` for the
current status.
