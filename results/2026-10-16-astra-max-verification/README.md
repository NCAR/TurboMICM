# Verification of astra-max, 2026-10-16

astra-max was rebuilt from its branch, retested, and audited, and all six finished
entries were timed in one interleaved session alongside MICM's default, MICM with
the entries' 32-cell groups, and MICM's fastest configuration, at 10 and 1
threads. The method, machine, input, and reference are those of
[`2026-10-15-astra-verification`](../2026-10-15-astra-verification/README.md),
whose builds of the other five entries and of MICM were reused.

| File | Contents |
|---|---|
| `tests.md` | `make test-cpu` results for astra-max with each compiler |
| `timings.csv` | one row per timed run: elapsed_ms, solver statistics, process wall time, gate |
| `summary.md` | median and range of elapsed_ms per program and thread count |
| `environment.txt` | the machine and the exact commands timed |
| `leaderboard.md` | `admin/leaderboard.py` at 10 threads over all four timing sessions |
| `leaderboard-1-thread.md` | the same at 1 thread, over the two sessions that ran it |

## Builds, tests, and audit

astra-max (`0dc7f67`, 7 commits, Codex CLI 0.163.0) was built with `make env`,
`make cpu bench OPTFLAGS="-O3 -DNDEBUG"`, and `make test-cpu`, once with
`CC=clang` and once with `CC=gcc-15`.

- **Tests** (`tests.md`): 680 PASS, 0 FAIL, 0 N/A with each compiler, every id
  exactly once, as the agent reported.
- **Flags, libraries, sources:** only `-std=c11`, warnings, `-O3 -DNDEBUG`,
  include paths, and OpenMP flags; the benchmark links only the OpenMP runtime
  and `libSystem`; no banned construct in the tracked or the generated C.
- **Same work, same clock:** `tools/ts1_bench.c` runs one untimed warm-up solve
  per group, then times every group's STEPS solves with `omp_get_wtime()` around
  one `omp parallel for`; input, setup, and output are outside the timer. Its
  default group is 32 cells, and it takes exactly MICM's 281,700 internal steps
  with 32-cell groups. Unlike MICM's driver and the other entries, it assigns
  groups to threads dynamically, one at a time (`dynamic,1`, overridable with
  `OMP_SCHEDULE`); the prompt leaves scheduling to the agent.

## Timing

`admin/time_entries.py` from 02:58 to 04:09: 5 rounds, each running all nine
programs at 10 threads and then at 1 thread, in rotating order, with a 5 s pause
between runs. Every one of the 90 runs passed the gate. The largest relative
differences from MICM's default were 2.96e-13 (astra-max), 2.03e-13 (astra-xhigh
and sol-high), 1.26e-13 (sol-max), and 0 (astra-high, sol-xhigh, and the three
MICM programs).

| Program | 10 threads, median ms (range) | vs MICM | vs MICM's fastest | 1 thread, median ms (range) | vs MICM | vs MICM's fastest |
|---|---:|---:|---:|---:|---:|---:|
| MICM, one group per thread | 27,549 (26,521–27,795) | 1.00× | | 144,459 (143,714–144,866) | 1.00× | |
| MICM, 32-cell groups | 23,280 (22,575–23,676) | 1.18× | | 135,465 (133,920–135,970) | 1.07× | |
| MICM's fastest: vector8, in-place Mozart, 32-cell groups | 9,767 (9,561–9,822) | 2.82× | 1.00× | 50,950 (50,661–51,070) | 2.84× | 1.00× |
| sol-max | 7,136 (6,942–7,203) | 3.86× | 1.37× | 37,625 (37,598–37,663) | 3.84× | 1.35× |
| astra-xhigh | 7,295 (7,095–7,349) | 3.78× | 1.34× | 38,479 (38,383–38,514) | 3.75× | 1.32× |
| sol-xhigh | 8,267 (8,060–8,281) | 3.33× | 1.18× | 44,222 (44,026–44,298) | 3.27× | 1.15× |
| astra-max | 8,307 (8,257–8,430) | 3.32× | 1.18× | 45,388 (45,296–45,457) | 3.18× | 1.12× |
| sol-high | 8,625 (8,469–8,753) | 3.19× | 1.13× | 46,790 (46,737–46,911) | 3.09× | 1.09× |
| astra-high | 9,084 (8,945–9,182) | 3.03× | 1.08× | 49,256 (49,215–49,392) | 2.93× | 1.03× |

## Findings

- **astra-max's correctness and portability claims hold,** with both compilers.
- **The order is sol-max, astra-xhigh, sol-xhigh, astra-max, sol-high,
  astra-high, then MICM's fastest configuration,** at both thread counts.
  sol-xhigh and astra-max are 0.5 % apart at ten threads, with overlapping
  ranges, and 2.6 % apart on one thread, where they do not overlap.
- **For Astra, max effort did not give the fastest code.** astra-max is 1.14
  times slower than astra-xhigh at ten threads (1.18 on one) and 1.09 times
  faster than astra-high (1.09 on one). It used the fewest tokens of the three
  (1.02 M, against 1.43 M at high and 1.88 M at xhigh) and the shortest session
  (1 h 41 min).
- **Scaling:** from 1 to 10 threads astra-max speeds up 5.46 times, the most of
  the entries (5.27–5.43×); it is the only one that schedules groups dynamically.
- **The measurements reproduce.** Every program's medians are within 1.2 % of
  the previous night's session, at both thread counts.
- **Self-report:** astra-max reported 8.24 s, against 8.31 s here.

Other load: none from the experiment. astra-max's Codex session, its work done,
stayed open and idle until it was closed during the session.

These use seed 2026, which the agents knew. The ranking itself will use an
undisclosed seed.
