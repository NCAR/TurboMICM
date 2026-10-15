# Verification of the Sol entries, 2026-10-15

The three Sol entries were rebuilt from their branches, retested, audited, and
timed alongside MICM on the run machine (MacBook Pro, Apple M5 with 4 Super and 6
Efficiency cores, 16 GB, macOS 27.0.1, on AC power and otherwise idle).

| File | Contents |
|---|---|
| `tests.md` | `make test-cpu` results for each entry and compiler |
| `timings.csv` | one row per timed run: elapsed_ms, solver statistics, process wall time, gate |
| `summary.md` | median and range of elapsed_ms per program and thread count |
| `environment.txt` | the machine and the exact commands timed |

## Builds and tests

Each entry was checked out from `runs/sol/<effort>` into a directory whose
`micm/` and `reference/` link to the run root's. It was built from its commit
with `make env`, then `make cpu bench`, then `make test-cpu`, once with
`CC=clang` and once with `CC=gcc-15`, both with `OPTFLAGS="-O3 -DNDEBUG"`.

- **Tests** (`tests.md`): every build reports 675 PASS, 0 FAIL, and 5 N/A, with
  each of the 680 case ids exactly once. The same five C++-only cases are N/A
  everywhere. This confirms the agents' own reports.
- **Flags:** the compile and link commands contain only `-std=c11`, warnings,
  `-O3 -DNDEBUG`, include paths, OpenMP flags, a group-size define, and (sol-max,
  for its shared library) `-fPIC`. None of them is a target, fast-math, LTO, or
  profile-guided flag.
- **Libraries:** each benchmark links only `libomp` and `libSystem`.
- **Sources:** no intrinsics headers, assembly, vector extensions,
  `long double`, or Apple frameworks.

## Is the same work timed the same way?

The timed region is the same in MICM's driver (`reference/micm_ts1_ref.cpp`) and
in all three entries' `src/bench.c`:

- a wall clock (`steady_clock` or `omp_get_wtime`) around the STEPS solves of
  every cell group, stopped when the last thread finishes;
- one untimed warm-up step per group before it;
- rate constants, building the mechanism and states, reading the input, and
  writing `ts1_final.csv` all outside it;
- static assignment of groups to threads.

The same input was used everywhere: 10,000 cells from `make_cells_csv.py` with
seed 2026, SHA-256 `43febf31…ad8b5`, the file the agents report timing. MICM's
final state for it matches, byte for byte, the MICM output sol-max generated
independently (SHA-256 `13329009…6ba7d`).

The one real difference was how cells are grouped. The entries solve 32-cell
groups (313 states, each with its own step-size control). MICM's driver
originally used one state per thread: all 10,000 cells on one thread, or 1,000
per thread on ten. So MICM was also timed with `--group-size 32` (`micm-g32`).
With identical grouping the work is identical: MICM and all three entries take
281,700 internal steps (9 per group per step, none rejected), with 563,400
function calls, 281,700 Jacobian updates, and 845,100 linear solves.

## Timing

`admin/time_entries.py`: 5 rounds, each running all five programs at 10 threads
and then at 1 thread, in an order that rotates from round to round, with a 5 s
pause between runs. Every one of the 50 runs passed the gate against MICM's final
state. The largest relative differences were 2.03e-13 (sol-high), 0 (sol-xhigh,
bit-identical to MICM), and 1.26e-13 (sol-max).

| Program | 10 threads, median ms (range) | vs MICM | vs MICM, 32-cell groups | 1 thread, median ms (range) | vs MICM | vs MICM, 32-cell groups |
|---|---:|---:|---:|---:|---:|---:|
| MICM, one group per thread | 27,561 (26,574–27,911) | 1.00× | | 144,630 (143,788–145,659) | 1.00× | |
| MICM, 32-cell groups | 23,906 (22,228–24,646) | 1.15× | 1.00× | 135,586 (133,792–138,680) | 1.07× | 1.00× |
| sol-high | 8,717 (8,409–8,836) | 3.16× | 2.74× | 46,748 (46,697–46,771) | 3.09× | 2.90× |
| sol-xhigh | 8,160 (7,949–8,420) | 3.38× | 2.93× | 44,172 (44,065–44,192) | 3.27× | 3.07× |
| sol-max | 7,159 (6,857–7,266) | 3.85× | 3.34× | 37,623 (37,560–37,647) | 3.84× | 3.60× |

All five programs scale alike from 1 to 10 threads: 5.25× (MICM), 5.67× (MICM,
32-cell groups), 5.36× (sol-high), 5.41× (sol-xhigh), and 5.26× (sol-max). One
thread is the steadier measurement: its runs spread by 0.2–0.3 % for the
entries and 1.3–3.6 % for MICM, against 5–10 % at 10 threads on the M5's mixed
cores.

## Findings

- **Correctness and portability claims hold.** Every test and gate result the
  agents reported was reproduced from their committed code, with both compilers.
- **The timing comparison was fair to MICM, apart from grouping.** Giving MICM
  the entries' 32-cell grouping makes it 7 % faster on one thread and 15 % on ten.
  Most of the entries' advantage, 2.7–3.6× over MICM with identical grouping and
  identical solver work, is in their own code.
- **The order is sol-max, sol-xhigh, sol-high,** at both thread counts. On one
  thread the three are clearly separated. At ten threads the ranges of sol-high
  and sol-xhigh overlap slightly, and their medians differ by 7 %.
- **Self-reported times were optimistic by 2–7 %.** The agents measured 8,084,
  7,714, and 6,935 ms (sol-max: median of 5) where these runs give medians of
  8,717, 8,160, and 7,159 ms. Their MICM measurements were also 0.5–6 % faster
  than today's, so the machine itself ran faster during their sessions. Their
  speedup ratios agree with these to within 8 %: 3.39×, 3.43×, and 3.75× against
  3.16×, 3.38×, and 3.85×.
- **MICM's context time has drifted.** The prompt cites 137,383 ms on one core,
  measured on 2026-10-14; today's median is 144,630 ms, 5 % slower. Comparisons
  made within one interleaved session are unaffected.

These use seed 2026, which the agents knew. The ranking itself will use an
undisclosed seed.
