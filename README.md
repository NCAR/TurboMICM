# CPU chemistry library

An ISO C11 port in progress of MICM's CPU chemistry functionality. Numerical
formulas and test cases derive from MICM; attribution is in NOTICE and LICENSE.

```
make env
make cpu CC=clang
make test-cpu CC=clang
```

The Python prefix environment and package cache live entirely in this directory.
`OPTFLAGS` defaults to `-O3 -DNDEBUG`. OpenMP flags are selected for the compiler
and operating system. To change compilers during development, run `make clean`
first, then build with `CC=gcc-15` (or the compiler available on your platform).

`make test-cpu` prints all 680 reference case identifiers and reports unported
cases as FAIL. The task is not yet complete. See RESULTS.md for verified progress.

Generate a mechanism builder:

```
.conda/bin/python tools/configure.py mechanisms/ts1.json build/ts1.c --symbol ts1
```

`make micm-ref` is the explicitly permitted reference-only C++ build. The C
library does not link MICM. The benchmark and remaining solvers are in progress.
