# Verification of luna-xhigh, 2026-10-20

luna-xhigh was rebuilt from its branch, retested, and audited, and all eight
finished entries were timed in one interleaved session alongside MICM's default,
MICM with 32-cell groups, and MICM's fastest configuration, at 10 and 1 threads.
The method, machine, input, reference, and builds are those of
[`2026-10-19-luna-max-verification`](../2026-10-19-luna-max-verification/README.md):
the same clang builds of the seven earlier entries and MICM's three programs in
`~/NSF-NCAR/tier1-work`, and the same input and MICM final state, whose SHA-256
sums match the earlier sessions'.

| File | Contents |
|---|---|
| `tests.md` | `make test-cpu` results for luna-xhigh with each compiler |
| `timings.csv` | one row per timed run: elapsed_ms, solver statistics, process wall time, gate |
| `summary.md` | median and range of elapsed_ms per program and thread count |
| `environment.txt` | the machine and the exact commands timed |
| `load-audit.log` | one sample a minute: load, AC power, every process above 5 % of a core |
| `leaderboard.md` | `admin/leaderboard.py` at 10 threads over all six timing sessions |
| `leaderboard-1-thread.md` | the same at 1 thread, over the four sessions that ran it |

## Builds, tests, and audit

luna-xhigh (`ada55d5`, 28 commits, Codex CLI 0.163.0) was built with `make env`,
`make cpu bench OPTFLAGS="-O3 -DNDEBUG"`, and `make test-cpu`, once with
`CC=clang` and once with `CC=gcc-15`, by `admin/verify_tier3.py --stages
build,tests,ts1,audit`.

- **Tests** (`tests.md`): 662 PASS, 0 FAIL, 18 N/A with each compiler, every id
  exactly once, as the agent reported. Eleven of the 18 N/A cases are N/A in
  no other entry: the five `Version` cases, a lambda rate entry,
  `SolverBuilder.ThrowsUnusedSpecies`, and four external-model constraint
  cases.
- **TS1 gate:** passes with both compilers; largest relative difference from
  MICM 1.77e-13 with clang (BCARY, cell 6848) and 2.33e-13 with gcc-15 (BCARY,
  cell 4179).
- **Flags, libraries, sources:** only `-std=c11`, warnings, `-O3 -DNDEBUG`,
  include paths, and OpenMP flags; every binary links only the OpenMP runtime
  and `libSystem`. The scan of the 41 C sources for banned constructs matched
  three lines, all a local variable named `vector_size` in
  `tests/test_grouped_views.c:164-170`, not GCC's vector attribute.
- **Warnings:** none, from the library, the benchmark, or the tests, with
  either compiler.
- **Same work, same clock:** `src/ts1_bench.c` runs one untimed warm-up
  interval, then times the STEPS intervals with C11's `timespec_get`
  (`src/ts1_bench.c:250-254`, `355-369`); file input and output are outside
  the timer. One `omp parallel for schedule(static)` over cell groups runs each
  group through all the timed intervals in one library call
  (`src/ts1_bench.c:357-365`), and that call, inside the timer, first builds
  the group's sparse ordering and fill, allocates its work arrays, copies the
  group into its species-major layout, and evaluates the rate constants, and
  at the end copies the group back (`src/chem.c:1080-1113`, `1249-1251`). MICM's
  driver does that kind of setup before its timer. It happens once per group
  per timed run, not per step; its cost was not measured. There is one
  group per thread (`src/ts1_bench.c:328`), as in MICM's default driver and
  luna-max, and it does MICM's default solver work: 9,000 internal steps at 10
  threads and 900 on one, with the same function calls, Jacobian updates,
  decompositions, and solves.

## Timing

The session waited for an idle machine (`time-when-idle.sh`): no Codex session
record written for 10 minutes, total CPU below one core, and AC power, for 10
checks in a row a minute apart. It waited from 12:34 to 12:45; luna-xhigh's
record last wrote at 12:20.

`admin/time_entries.py` from 12:45 to 14:24: 5 rounds, each running all eleven
programs at 10 threads and then at 1 thread, in rotating order, with a 5 s pause
between runs. Every one of the 110 runs passed the gate. The largest relative
differences from MICM's default were 2.96e-13 (astra-max), 2.78e-13
(luna-max), 2.03e-13 (astra-xhigh and sol-high), 1.77e-13 (luna-xhigh),
1.26e-13 (sol-max), and 0 (astra-high, sol-xhigh, and the three MICM programs).

| Program | 10 threads, median ms (range) | vs MICM | vs MICM's fastest | 1 thread, median ms (range) | vs MICM | vs MICM's fastest |
|---|---:|---:|---:|---:|---:|---:|
| MICM, one group per thread | 28,198 (26,177–28,474) | 1.00× | | 144,435 (143,680–145,585) | 1.00× | |
| MICM, 32-cell groups | 23,305 (22,774–23,421) | 1.21× | | 135,977 (134,396–136,905) | 1.06× | |
| MICM's fastest: vector8, in-place Mozart, 32-cell groups | 9,743 (9,627–9,808) | 2.89× | 1.00× | 50,987 (50,064–51,075) | 2.83× | 1.00× |
| sol-max | 7,097 (7,081–7,266) | 3.97× | 1.37× | 37,661 (37,597–37,679) | 3.84× | 1.35× |
| astra-xhigh | 7,369 (7,264–7,478) | 3.83× | 1.32× | 38,275 (38,247–38,360) | 3.77× | 1.33× |
| sol-xhigh | 8,188 (8,006–8,358) | 3.44× | 1.19× | 44,211 (44,096–44,214) | 3.27× | 1.15× |
| astra-max | 8,545 (8,185–8,663) | 3.30× | 1.14× | 45,308 (45,271–45,362) | 3.19× | 1.13× |
| sol-high | 8,710 (8,603–8,735) | 3.24× | 1.12× | 46,802 (46,759–46,839) | 3.09× | 1.09× |
| astra-high | 9,232 (8,923–9,449) | 3.05× | 1.06× | 49,313 (49,219–49,349) | 2.93× | 1.03× |
| luna-xhigh | 21,443 (21,374–21,506) | 1.32× | 0.45× | 59,916 (59,903–59,993) | 2.41× | 0.85× |
| luna-max | 31,551 (30,857–31,756) | 0.89× | 0.31× | 183,302 (182,864–183,504) | 0.79× | 0.28× |

## Findings

- **luna-xhigh's correctness and portability claims hold,** with both compilers.
- **luna-xhigh is faster than MICM's default,** which does the same solver work
  with the same grouping: 1.32 times at ten threads and 2.41 times on one. It
  is 1.47 times faster than luna-max at ten threads and 3.06 times on one.
  MICM's fastest configuration is 2.20 times faster than luna-xhigh at ten
  threads and 1.18 times on one; sol-max is 3.02 and 1.59 times faster.
- **Scaling:** from 1 to 10 threads luna-xhigh speeds up 2.79 times, the least
  of any program in the session (the other entries 5.19–5.81×, MICM's three
  programs 5.12–5.83×). Its solver work at each thread count matches MICM's
  default, so the difference lies in how its time per step grows with threads;
  the analysis pass is to measure why.
- **The other seven keep their order,** sol-max, astra-xhigh, sol-xhigh,
  astra-max, sol-high, astra-high, luna-max, at both thread counts.
- **The measurements reproduce.** Every program's medians are within 2.2 % of
  the 2026-10-19 session at ten threads and within 0.6 % on one thread.
- **Self-report:** luna-xhigh reported 21,397 ms against MICM's 27,851 ms in one
  ten-thread run, against 21,443 and 28,198 ms here.

Other load (`load-audit.log`, 99 samples): AC attached throughout, the battery
held at 80 % by the system. Processes above 5 % of a core other than the
benchmark: macOS's media analysis daemon in 7 samples (mean 0.7 core, at most
1.0), Spotlight's indexer in 3 (0.1–0.2 core), the timing script's gate check, in
Python, in 1 (0.4 core), Claude Code in 1 (0.1 core), and the biome sync daemon in 1
(0.1 core). No agent ran: luna-xhigh's Codex session had been quit at 12:31;
this Claude Code session waited, apart from one read-only look at the logs at
13:03; and the agent the SSE had working on the wiki's slides in another
session showed no process on this machine above 5 % of a core. The rotation
spreads the daemons' load over every program.

These use seed 2026, which the agents knew. The ranking itself will use an
undisclosed seed.
