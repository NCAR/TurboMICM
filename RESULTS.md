# Results

## Tier 2: mechanism compiler

Validated on 2026-10-22 on the supplied Apple M5 Pro (18 cores, 24 GB),
macOS 26.7.1. Compilers: Apple clang 21.0.0 and Homebrew GCC 15.3.0.
All C builds used `OPTFLAGS='-O3 -DNDEBUG'` plus the Makefile's language,
warning, include, PIC, OpenMP, and link flags. The MICM baseline used clang++,
C++20, and the same optimization flags. No compiler optimization flags were tuned.

The selected generator settings are 32-cell groups, 32 operation units per
function, four units per SIMD loop, a full-group cell tile, constant-stride
variants enabled, and Mozart LU. OpenMP uses dynamic group scheduling by default.

### Tests and gates

| Check | clang | GCC |
| --- | --- | --- |
| `make cpu` and `make bench` | PASS | PASS |
| `make test-cpu` | 680 PASS, 0 FAIL, 0 N/A | 680 PASS, 0 FAIL, 0 N/A |
| `make test-specialized` | All six mechanisms PASS | All six mechanisms PASS |
| Generic and specialized 100-cell MICM gate | PASS | PASS |
| Generic and specialized 10,000-cell MICM gate, every timing run | PASS | PASS |
| All five disclosed day-night scenarios, final states and snapshots | PASS | PASS |
| Additional 48-hour and interpolation scenarios | PASS | PASS |

Generator tests cover TS1, the three Tier 1 configuration-file mechanisms, all
supported rate families together, and a reaction-free mechanism. Each checks
rates, forcing, negative Jacobian, LU, direct solves, and five evolving steps
against the generic implementation at `1e-10` relative plus a `1e-30` absolute
floor. They exercise 1, 3, 7, 32, and 33 cells, default and tight tolerances,
forced initial rejections, and changed conditions, independent air density, and
custom parameters. The two compiler builds generated identical source hashes.
Additional parameter checks pass with seven-cell fixed groups, Doolittle LU,
five-unit functions and three-unit loops on both compilers, and with runtime
strides only, 64-unit functions and nine-unit loops on clang. These results are
recorded in `validation/parameter-checks.json`.

The full benchmark gate checks all 2,110,000 values against a freshly run MICM
reference. The largest differences remain below `3e-13` relative. Complete
per-case results are in `validation/clang-cases.txt` and `validation/gcc-cases.txt`;
timings, exact gate output, source hashes, and both snapshot/final scenario scores
are in [`validation/clang-tier2.json`](validation/clang-tier2.json) and
[`validation/gcc-tier2.json`](validation/gcc-tier2.json).

### Full TS1 benchmark

Workload: 10,000 cells, seed 2026, 100 timed steps of 30 seconds, after one
untimed warm-up. All 18-thread figures are medians of three sequential runs;
one-thread figures are single observations. No build or test ran concurrently
with a timing run. C solvers use groups of 32; MICM uses its benchmark default
of one group per worker. The supplied cell generator's output is identical to
the Tier 1 seed-2026 input (checksum below).

| Program | Compiler | 18 threads, seconds | 1 thread, seconds |
| --- | --- | ---: | ---: |
| MICM | clang++ | 11.906 | 136.212 |
| Generic C | clang | 3.872 | 45.126 |
| Specialized C | clang | **3.063** | **37.398** |
| Generic C | GCC | 3.982 | 52.858 |
| Specialized C | GCC | 3.379 | 39.581 |

The specialized clang solver is **1.26× faster than the generic library and
3.89× faster than MICM on 18 threads**. On one thread the corresponding ratios
are 1.21× and 3.64×. These are local measurements on a disclosed seed, separate
from the experimenter's idle-machine, undisclosed-seed ranking.

Every C run takes 281,700 internal steps across 313 groups, with 281,700 accepted,
zero counted rejections, 563,400 function calls, 281,700 Jacobian updates and
factorizations, and 845,100 triangular solves. Each cell follows the same nine
internal steps per timed interval. MICM's different number of groups changes
its summed counters, as in Tier 1.

### Build and generated size

| Clean `make -j6 bench` | Seconds |
| --- | ---: |
| clang | 12.589 |
| GCC | 43.713 |

These wall times start after `make clean`, with the conda environment already
present. They include configuration generation, the generic library/benchmark,
and both specialized executables; reference C++ and test builds are separate.
The two TS1 generated C files total **2,863,137 bytes** (2.73 MiB), including
the constructor, dispatch, and numerical kernels. Headers and JSON metadata
are not included in that byte count. Generated C is a reproducible build product.

### Day-night results

All runs use consecutive 32-cell groups, relative tolerance `1e-6`, and absolute
tolerance `1e-15`. The following scores compare all four disclosed snapshots
(4,000 rows per scenario) with MICM. Both snapshot and final-state files pass
`score.py --level 10`, with no nonfinite value and none below `-1e-20`.

| Scenario | clang overall digits | clang worst species (digits) | GCC overall digits | GCC worst species (digits) |
| --- | ---: | --- | ---: | --- |
| diurnal-30min | 14.365 | ISOP (13.780) | 14.363 | ISOP (13.791) |
| diurnal-10min | 14.288 | ISOP (13.717) | 14.292 | ISOP (13.687) |
| diurnal-5min | 14.229 | ISOP (13.643) | 14.232 | ISOPNO3 (13.537) |
| polluted-10min | 14.071 | ISOPNO3 (13.428) | 14.050 | ISOP (13.432) |
| clean-10min | 14.334 | ISOP (13.813) | 14.335 | ISOP (13.811) |

The reference-solution scores below are the same to the displayed precision
with either compiler. They reflect the accuracy of MICM's three-stage method
at the requested tolerances, rather than a change in numerical method.

| Scenario | Internal steps (both compilers and MICM) | clang ms | GCC ms | Reference overall digits | Reference worst species (digits) |
| --- | ---: | ---: | ---: | ---: | --- |
| diurnal-30min | 48,619 | 539.574 | 570.636 | 4.748 | ISOP (3.610) |
| diurnal-10min | 90,980 | 1,020.266 | 1,060.500 | 5.052 | ISOP (3.824) |
| diurnal-5min | 145,891 | 1,626.636 | 1,734.086 | 5.214 | ISOP (4.328) |
| polluted-10min | 106,393 | 1,189.592 | 1,268.854 | 5.398 | ISOP (4.257) |
| clean-10min | 84,273 | 949.237 | 980.386 | 4.907 | ISOP (3.766) |

Additional checks use independently run MICM three-stage results:

- A 48-hour summer case: 65 cells, seed 4242, 55° latitude, day 172, start at
  19:00, NOx multiplied by three, 192 steps of 900 s, snapshots every six hours.
  Clang agrees to 14.246 digits overall, worst species 13.871; GCC to 14.270,
  worst 13.903. The final group contains one cell.
- An interpolation case: 33 cells, eight 450-second steps, photolysis samples
  at 600, 1,234, and 2,200 s, snapshots at steps 3, 6, and 8. It tests both
  endpoint clamps, interpolation between rows, and an irregular final snapshot.
  Clang agrees to 14.866 digits overall, worst 14.464; GCC to 14.856,
  worst 14.492.

Omitting the snapshot interval or setting it longer than the run also produces
the correct final-only snapshot and passes the MICM final-state score; see
`validation/scenario-format-checks.json`. Source inspection found no prohibited
extensions, intrinsics, or precision types. The default native executables link
only the system C runtime and OpenMP; the generator environment contains the
specified Python 3.12.12, NumPy 2.2.6, and PyYAML 6.0.3.

To reproduce all clean builds, gates, timings, and scenario scores:

```sh
make env
.conda/bin/python tools/validate_tier2.py
```

This writes only in this checkout. The supplied `../reference`, `../scenarios`,
and `../micm` inputs are read-only. Individual commands and generator tuning
controls are documented in [README.md](README.md) and [DESIGN.md](DESIGN.md);
all retained and discarded tuning experiments are in
[SPECIALIZATION.md](SPECIALIZATION.md). Linux/x86-64 execution was unavailable;
the source restrictions and both supplied compilers were checked on this host.

## Tier 1 historical results

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
