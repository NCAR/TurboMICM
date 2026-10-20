# Verification of luna-max, 2026-10-19

luna-max was rebuilt from its branch, retested, and audited, and all seven
finished entries were timed in one interleaved session alongside MICM's default,
MICM with 32-cell groups, and MICM's fastest configuration, at 10 and 1 threads.
The method, machine, input, and reference are those of
[`2026-10-16-astra-max-verification`](../2026-10-16-astra-max-verification/README.md).
The entries and MICM were rebuilt for this session in `~/NSF-NCAR/tier1-work`
(`prepare.sh` there): the six earlier entries with clang, from their branches,
by `admin/verify_tier3.py --stages build,ts1`, each passing the gate before
the session; MICM's programs from `admin/micm_configs_source.py`'s copy of
`reference/micm_ts1_ref.cpp`, as `admin/micm_configs.sh` builds them. The input
and MICM's final state match the earlier sessions' SHA-256 sums.

| File | Contents |
|---|---|
| `tests.md` | `make test-cpu` results for luna-max with each compiler |
| `timings.csv` | one row per timed run: elapsed_ms, solver statistics, process wall time, gate |
| `summary.md` | median and range of elapsed_ms per program and thread count |
| `environment.txt` | the machine and the exact commands timed |
| `load-audit.log` | one sample a minute: load, AC power, every process above 5 % of a core |
| `leaderboard.md` | `admin/leaderboard.py` at 10 threads over all five timing sessions |
| `leaderboard-1-thread.md` | the same at 1 thread, over the three sessions that ran it |

## Builds, tests, and audit

luna-max (`da2d722`, 5 commits, Codex CLI 0.163.0) was built with `make env`,
`make cpu bench OPTFLAGS="-O3 -DNDEBUG"`, and `make test-cpu`, once with
`CC=clang` and once with `CC=gcc-15`, by `admin/verify_tier3.py --stages
build,tests,ts1,audit`.

- **Tests** (`tests.md`): 671 PASS, 0 FAIL, 9 N/A with each compiler, every id
  exactly once, as the agent reported.
- **TS1 gate:** passes with both compilers; largest relative difference from
  MICM 2.78e-13 with clang and 2.2e-13 with gcc-15 (BCARY, cell 4179).
- **Flags, libraries, sources:** only `-std=c11`, warnings, `-O3 -DNDEBUG`,
  include paths, and OpenMP flags; every binary links only the OpenMP runtime
  and `libSystem`; no banned construct in the 22 C sources.
- **Warnings:** the library and benchmark build without warnings under clang;
  gcc-15 gives three `-Wint-in-bool-context` warnings, on allocation-size
  guards of the form `n * m ? n * m : 1` (`src/solver.c:837-838`,
  `src/benchmark.c:307`). The tests add 38 `-Wmissing-field-initializers`
  warnings with each compiler (and two `-Wmissing-braces` with gcc-15), for
  partly initialized structs in `tests/test_cpu.c` and
  `tests/test_example_runtime.c`.
- **Same work, same clock:** `src/benchmark.c` runs one untimed warm-up step,
  then times the STEPS calls of `chem_solver_step` with `omp_get_wtime()`; input,
  setup, and output are outside the timer. Each step is one `omp parallel for
  schedule(static)` over cell groups (`src/solver.c:980`). With no group size
  given, there is one group per thread (`src/solver.c:830-832`), as in MICM's
  default driver, and it takes MICM's default solver work: 9,000 internal steps
  at 10 threads and 900 on one, with the same function calls, Jacobian
  updates, decompositions, and solves. It is the only entry not to use 32-cell
  groups by default.

## Timing

The session waited for an idle machine (`time-when-idle.sh`): no Codex session
record written for 10 minutes, total CPU below one core, and AC power, for 10
checks in a row a minute apart. It waited from 21:25 to 23:06: the SSE's
separate Codex session last wrote at 21:34, and then Spotlight's knowledge
indexer used a full core from 21:47 until 22:57.

`admin/time_entries.py` from 23:06 to 00:36: 5 rounds, each running all ten
programs at 10 threads and then at 1 thread, in rotating order, with a 5 s pause
between runs. Every one of the 100 runs passed the gate. The largest relative
differences from MICM's default were 2.96e-13 (astra-max), 2.78e-13
(luna-max), 2.03e-13 (astra-xhigh and sol-high), 1.26e-13 (sol-max), and 0
(astra-high, sol-xhigh, and the three MICM programs).

| Program | 10 threads, median ms (range) | vs MICM | vs MICM's fastest | 1 thread, median ms (range) | vs MICM | vs MICM's fastest |
|---|---:|---:|---:|---:|---:|---:|
| MICM, one group per thread | 27,876 (25,847–28,582) | 1.00× | | 144,522 (143,739–145,253) | 1.00× | |
| MICM, 32-cell groups | 23,497 (22,284–23,905) | 1.19× | | 135,382 (133,882–136,753) | 1.07× | |
| MICM's fastest: vector8, in-place Mozart, 32-cell groups | 9,808 (9,195–10,023) | 2.84× | 1.00× | 50,673 (50,599–51,211) | 2.85× | 1.00× |
| sol-max | 7,115 (6,947–7,224) | 3.92× | 1.38× | 37,613 (37,573–37,694) | 3.84× | 1.35× |
| astra-xhigh | 7,297 (7,191–7,418) | 3.82× | 1.34× | 38,312 (38,157–38,375) | 3.77× | 1.32× |
| sol-xhigh | 8,171 (7,938–8,272) | 3.41× | 1.20× | 44,191 (44,156–44,258) | 3.27× | 1.15× |
| astra-max | 8,408 (8,090–8,446) | 3.32× | 1.17× | 45,392 (45,294–45,488) | 3.18× | 1.12× |
| sol-high | 8,621 (8,459–8,736) | 3.23× | 1.14× | 46,757 (46,640–46,854) | 3.09× | 1.08× |
| astra-high | 9,030 (9,000–9,270) | 3.09× | 1.09× | 49,532 (49,227–49,611) | 2.92× | 1.02× |
| luna-max | 30,914 (30,484–31,633) | 0.90× | 0.32× | 183,244 (182,292–184,055) | 0.79× | 0.28× |

## Findings

- **luna-max's correctness and portability claims hold,** with both compilers.
- **luna-max is slower than MICM's default,** which does the same solver work
  with the same grouping: 1.11 times slower at ten threads and 1.27 times
  slower on one. MICM's fastest configuration is 3.15 times faster at ten
  threads (3.62 on one), and sol-max 4.35 times (4.87 on one).
- **Scaling:** from 1 to 10 threads luna-max speeds up 5.93 times, the most of
  any program in the session (the other entries 5.25–5.49×, MICM's default
  5.18×).
- **The other six keep their order,** sol-max, astra-xhigh, sol-xhigh,
  astra-max, sol-high, astra-high, at both thread counts. sol-xhigh and
  astra-max are 2.9 % apart at ten threads (0.5 % in the previous session) and
  2.7 % apart on one.
- **The measurements reproduce.** Every program's medians are within 1.3 % of
  the 2026-10-16 session at ten threads and within 0.6 % on one thread.
- **Self-report:** luna-max reported 185,107 ms on one thread, against 183,244
  here, and that MICM was faster than its code, as measured here.

Other load (`load-audit.log`, 91 samples): AC attached throughout, the battery
held at 80 % by the system. Processes above 5 % of a core other than the
benchmark: macOS's media analysis daemon in 22 samples (mean 0.9 core, at most
1.2), the PDF service in 5 (0.5 core), the Wi-Fi daemon in 6 and the system
policy daemon in 4 (0.1 core each), and seven others in one sample each. No
agent ran: the SSE's Codex session stayed
open but idle, and this Claude Code session waited without running commands.
The rotation spreads the daemons' load over every program.

These use seed 2026, which the agents knew. The ranking itself will use an
undisclosed seed.
