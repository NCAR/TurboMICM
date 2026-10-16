# Verification of the Astra entries, 2026-10-15 and 16

astra-xhigh was rebuilt from its branch, retested, and audited, and all five
finished entries were timed in one interleaved session alongside MICM's default,
MICM with the entries' 32-cell groups, and MICM's fastest configuration, at 10
and 1 threads. This completes astra-high's timing, which had been measured only
at 10 threads, with the 24 configurations, and adds the 1-thread comparison with
MICM's fastest configuration. The machine is the run machine (MacBook Pro, Apple
M5 with 4 Super and 6 Efficiency cores, 16 GB, macOS 27.0.1), on AC power.

| File | Contents |
|---|---|
| `tests.md` | `make test-cpu` results for astra-xhigh with each compiler |
| `timings.csv` | one row per timed run: elapsed_ms, solver statistics, process wall time, gate |
| `summary.md` | median and range of elapsed_ms per program and thread count |
| `environment.txt` | the machine and the exact commands timed |
| `leaderboard.md` | `admin/leaderboard.py` at 10 threads over all three timing sessions of 2026-10-15 |
| `leaderboard-1-thread.md` | the same at 1 thread, from this session alone |

## Builds, tests, and audit

astra-xhigh (`fbd365c`, 11 commits, Codex CLI 0.163.0) was checked out from
`runs/astra/xhigh` into a directory whose `micm/` and `reference/` link to the run
root's, and built with `make env`, `make cpu bench OPTFLAGS="-O3 -DNDEBUG"`, and
`make test-cpu`, once with `CC=clang` and once with `CC=gcc-15`, with `SDKROOT`
left to the Makefile.

- **Tests** (`tests.md`): 680 PASS, 0 FAIL, 0 N/A with each compiler, every id
  exactly once. This confirms the agent's own report.
- **Flags:** the compile and link commands contain only `-std=c11`, warnings,
  `-O3 -DNDEBUG`, include paths, and OpenMP flags. None is a target, fast-math,
  LTO, or profile-guided flag.
- **Libraries:** the benchmark links only the OpenMP runtime (`libomp` with clang,
  `libgomp` with gcc-15) and `libSystem`.
- **Sources:** no intrinsics headers, assembly, vector extensions, `long double`,
  target-specific attributes or pragmas, or Apple frameworks, in the tracked or
  the generated C.
- **Same work, same clock:** `src/bench.c` runs one untimed warm-up solve per
  group, then times the STEPS solves of every group with `omp_get_wtime()`
  around one `omp parallel for`, with static scheduling by default; input, setup,
  and output are outside the timer. That matches MICM's driver and the other
  entries. Its default group is 32 cells, and it takes exactly the same solver
  work as MICM with 32-cell groups: 281,700 internal steps, none rejected, 563,400
  function calls, 281,700 Jacobian updates, and 845,100 linear solves.

The four earlier entries were rebuilt with clang from their branches, and MICM's
two configurations from the same generated copy of `reference/micm_ts1_ref.cpp`
as the configuration sessions (`admin/micm_configs.sh`), which is byte-identical.
The input was regenerated with `make_cells_csv.py 10000 2026` (SHA-256
`43febf31…ad8b5`), and MICM's default result on it matches the recorded one
(SHA-256 `13329009…6ba7d`).

## Timing

`admin/time_entries.py` from 23:49 to 00:55: 5 rounds, each running all eight
programs at 10 threads and then at 1 thread, in an order that rotates from round
to round, with a 5 s pause between runs. Every one of the 80 runs passed the gate.
The largest relative differences from MICM's default were 2.03e-13 (astra-xhigh
and sol-high), 1.26e-13 (sol-max), and 0 (astra-high, sol-xhigh, and all three
MICM programs).

| Program | 10 threads, median ms (range) | vs MICM | vs MICM, 32-cell groups | vs MICM's fastest | 1 thread, median ms (range) | vs MICM | vs MICM, 32-cell groups | vs MICM's fastest |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| MICM, one group per thread | 27,782 (26,717–27,836) | 1.00× | | | 144,668 (143,832–145,461) | 1.00× | | |
| MICM, 32-cell groups | 23,407 (22,855–23,735) | 1.19× | 1.00× | | 136,210 (134,031–136,927) | 1.06× | 1.00× | |
| MICM's fastest: vector8, in-place Mozart, 32-cell groups | 9,797 (9,762–9,956) | 2.84× | 2.39× | 1.00× | 50,659 (50,087–50,959) | 2.86× | 2.69× | 1.00× |
| sol-max | 7,175 (6,945–7,250) | 3.87× | 3.26× | 1.37× | 37,684 (37,610–37,698) | 3.84× | 3.61× | 1.34× |
| astra-xhigh | 7,277 (7,116–7,361) | 3.82× | 3.22× | 1.35× | 38,328 (38,095–38,365) | 3.77× | 3.55× | 1.32× |
| sol-xhigh | 8,238 (7,886–8,359) | 3.37× | 2.84× | 1.19× | 44,105 (44,079–44,252) | 3.28× | 3.09× | 1.15× |
| sol-high | 8,729 (8,331–8,747) | 3.18× | 2.68× | 1.12× | 46,702 (46,663–46,828) | 3.10× | 2.92× | 1.08× |
| astra-high | 9,049 (8,857–9,211) | 3.07× | 2.59× | 1.08× | 49,247 (49,180–49,561) | 2.94× | 2.77× | 1.03× |

From 1 to 10 threads the programs scale by 5.2–5.8×: 5.21× (MICM), 5.82× (MICM,
32-cell groups), 5.17× (MICM's fastest), and 5.25–5.44× for the entries. On one
thread the entries' runs spread by 0.2–0.8 % and MICM's by 1.1–2.1 %; on ten, all
spread by 2.0–5.7 %.

## Findings

- **astra-xhigh's correctness and portability claims hold,** with both compilers.
- **The order is sol-max, astra-xhigh, sol-xhigh, sol-high, astra-high, then
  MICM's fastest configuration,** at both thread counts. astra-xhigh is 1.4 %
  slower than sol-max at ten threads, where their ranges overlap, and 1.7 % slower
  on one thread, where they do not.
- **Effort helped Astra:** xhigh is 1.24 times faster than high at ten threads
  and 1.28 times on one, and it used 31 % more tokens (1.88 M against 1.43 M).
- **MICM's fastest configuration keeps its lead over the default on one thread**
  (2.86×, against 2.84× at ten), and every entry still beats it: by 1.03× (astra-high)
  to 1.34× (sol-max) on one thread, and 1.08× to 1.37× at ten.
- **The measurements reproduce.** The Sol entries' medians are within 1.0 % of
  their first verification session's, at both thread counts, and MICM's 1-thread
  median (144,668 ms) is within 0.03 % of that session's (144,630 ms).
- **Self-reports:** astra-xhigh reported 7.31 s against MICM's 28.08 s; these runs
  give 7.28 s and 27.78 s. astra-high had reported 9.63 s, 6 % slower than its
  9.05 s here.

Other load: the astra-xhigh Codex session stayed open and idle. Claude Code's own
session made one brief progress check; nothing else from the experiment ran.

These use seed 2026, which the agents knew. The ranking itself will use an
undisclosed seed.
