# Results

Validated on 2026-10-16 on the supplied Apple M5 machine (10 CPU cores, 16 GB
memory, macOS 27.0.1). C compilers were Apple clang 21.0.0 and Homebrew GCC
15.2.0 (`gcc-15`). Every C build used `OPTFLAGS='-O3 -DNDEBUG'`, with the
Makefile's C11, warning, include, PIC, and OpenMP flags. The reference was built
by `make micm-ref` with clang++, C++20, and the same `-O3 -DNDEBUG` flags.

## CPU cases

| Build | PASS | FAIL | N/A |
| --- | ---: | ---: | ---: |
| `make test-cpu CC=clang` | 680 | 0 | 0 |
| `make test-cpu CC=gcc-15` | 680 | 0 | 0 |

These are all 672 GoogleTest IDs and eight plain-program IDs in the supplied
case list. The committed per-case records are
[`validation/clang-cases.txt`](validation/clang-cases.txt) and
[`validation/gcc-cases.txt`](validation/gcc-cases.txt). Their IDs were checked
against the supplied list for missing entries and duplicates. Final builds had
no compiler warnings.

The runs include all three generated tutorial configurations compared with their
hand-built equivalents, both sparse LU algorithms, all five Rosenbrock tables,
backward Euler, differential/algebraic constraints, external processes, matrix
and view policies, metadata, ownership/error handling, analytical and regression
solutions, and native OpenMP execution. The native controller/linear-solver
replacement check also passed with both compilers. The complete 3,601-row
wall-loss trajectories pass the original absolute tolerance of 1e-10.

All eight plain-program comparisons match their deterministic reference output.
Timing values and concurrent startup ordering are handled as specified in
[`DECISIONS.md`](DECISIONS.md); numerical output and printed solver counters are
compared exactly. There are no N/A exceptions.

## TS1 numerical gate

All runs use the default three-stage Rosenbrock parameters and tolerances,
one warm-up step, then 100 steps of 30 seconds. The supplied comparison script
checks every column by name, using
`abs(y - reference) <= 1e-30 + 1e-10 * abs(reference)`.

| C compiler | Input | Checked values | Largest relative difference | Result |
| --- | --- | ---: | ---: | --- |
| clang | 100 cells, seed 1 | 21,100 | 3.23e-14 | PASS |
| GCC | 100 cells, seed 1 | 21,100 | 3.66e-14 | PASS |
| clang | 10,000 cells, seed 2026 | 2,110,000 | 2.96e-13 | PASS |
| GCC | 10,000 cells, seed 2026 | 2,110,000 | 2.75e-13 | PASS |

The 211 columns are temperature, pressure, and 209 evolving species. TS1 has
210 declared species and 547 reactions; the remaining species is the
parameterized third body. The small run uses the supplied MICM final state;
the full run uses a newly executed reference built with the flags above.
Every repeated native timing run also passed the full comparison. The largest
relative difference in the full runs occurred for BCARY in cell 4179.

Additional clang gate runs on the 100-cell input passed with group sizes
1, 3, 7, 31, 32, 64, and 257 (the last capped at 100), thread counts 1, 2, 3,
and 10, and static/dynamic/guided scheduling with different chunk sizes. The
largest relative difference stayed 3.23e-14. The full 10,000-cell tuning sweep
also passed for group sizes 16, 32, 64, 128, and 256. Group assignment does not
change any group's arithmetic.

The retained translator was rerun from the supplied TS1 source, and its output
was byte-for-byte identical to the committed `configs/ts1.json`.

## Timing

Workload: **10,000 distinct cells, seed 2026, 100 timed steps, DT = 30 s**.
Times are the programs' `elapsed_ms`, excluding initialization, the single
warm-up solve, and output. Runs were executed sequentially, with no concurrent
build or test job. The C implementation uses groups of 32 and dynamic scheduling;
MICM uses its reference program's default of one group per worker. Both use
all 10 cores unless explicitly restricted.

| Program | Threads | Samples | Median elapsed_ms | Speedup over corresponding MICM run |
| --- | ---: | ---: | ---: | ---: |
| C, clang | 10 | 5 | **8,236.661** | **3.25x** |
| C, GCC | 10 | 3 | 9,556.184 | 2.80x |
| MICM, clang++ | 10 | 3 | 26,803.373 | 1.00x |
| C, clang | 1 | 1 | 45,282.753 | 3.23x |
| MICM, clang++ | 1 | 1 | 146,247.295 | 1.00x |

Raw all-core times in milliseconds:

- C/clang: 7,772.185; 8,436.673; 8,308.654; 8,236.661; 8,227.012.
- C/GCC: 9,556.184; 9,585.853; 9,465.783.
- MICM/clang++: 26,749.312; 27,120.543; 26,803.373.

The single-core figures are single observations for comparison, not medians
from repeated runs. These local measurements are separate from the experiment's
idle-machine, undisclosed-seed ranking procedure. Linux/x86-64 was not available
for execution; portability was checked through the source restrictions and the
two supplied compilers on this host.

Each cell took nine internal steps per timed solve, with no rejections. Counters
are summed over independent states, so their totals differ with group count:

| Counter | C, 313 groups | MICM, 10 groups | MICM, one group |
| --- | ---: | ---: | ---: |
| internal_steps | 281,700 | 9,000 | 900 |
| accepted | 281,700 | 9,000 | 900 |
| rejected | 0 | 0 | 0 |
| function_calls | 563,400 | 18,000 | 1,800 |
| jacobian_updates | 281,700 | 9,000 | 900 |
| decompositions | 281,700 | 9,000 | 900 |
| solves | 845,100 | 27,000 | 2,700 |

[`validation/benchmarks.json`](validation/benchmarks.json) retains the raw
elapsed times, counters, and native gate results. Reproduction commands are in
[`README.md`](README.md). Environment recreation/update was checked with
`make env`; the final clang executable was rebuilt and matched the binary used
for the recorded measurements byte for byte. Dynamic dependency inspection
showed only the compiler's OpenMP runtime and the system C runtime for both C
executables, with no Python, MICM, or MUSICA dependency.

Input SHA-256 checksums:

```text
configs/ts1.json:
73843ae0472aa6e69f44a75967e5059d1a592276c58c8c761e4f5ec755a1e65e
tests/cases.txt:
95615128d26687fbf371e99551062160d9d0c0aff4e5c8c36ca38e4dec840da4
10,000-cell seed-2026 CSV generated by the supplied script:
43febf31a0a1f8e0acfe02d93ac96d95ecc024df9572547dd1000926a85ad8b5
```
