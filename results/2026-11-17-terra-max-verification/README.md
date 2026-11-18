# Verification of terra-max, 2026-11-17

terra-max, the first Terra entry (gpt-6.1-terra at max effort, Codex CLI
0.169.0), was rebuilt from its branch, retested, and audited with both
compilers on 2026-11-17 at 01:37, twelve minutes after the run ended. It was
timed with every other entry that evening, in one session (the SSE's choice,
2026-11-16; see [Timing](#timing)).

| File | Contents |
|---|---|
| `tests.md` | `make test-cpu` results for terra-max with each compiler |

## Builds, tests, and audit

terra-max (`9071941`, 22 commits) was built with `make env`, `make cpu bench
OPTFLAGS="-O3 -DNDEBUG"`, and `make test-cpu`, once with `CC=clang` (Apple clang
21.0.0, clang-2100.0.123.102, the compiler that built the Luna entries and
MICM's timed programs) and once with `CC=gcc-15` (Homebrew GCC 15.2.0), by
`admin/verify_tier3.py --stages build,tests,ts1,audit`, the method of
[`2026-10-20-luna-high-verification`](../2026-10-20-luna-high-verification/README.md).

- **Tests** (`tests.md`): 676 PASS, 0 FAIL, 4 N/A with each compiler, every id
  exactly once, as the agent reported. The four N/A cases are MICM's checks of
  C++ template concepts, which no C port can express. `tests/run_cases.py`
  takes each PASS or FAIL line from one of 27 C test programs (and one Python
  ordering test), each printing one line per case from a check run for that
  case: `tests/reference_lu.c:139-163`, for example, runs a factor or solve
  reconstruction for each of its 32 LU and linear-solver ids, parameterized by
  the case's matrix family, ordering, and variant. Eight tutorial cases compare
  printed output with MICM's expected text, with timings and thread numbers
  masked (`tests/run_cases.py:34-58`). A case it does not report fails as
  "case has not yet been ported".
- **TS1 gate:** passes with both compilers; largest relative difference from
  MICM 2.51e-13 with clang and 1.57e-13 with gcc-15 (BCARY, cell 4179, in
  both).
- **Flags, libraries, sources:** only `-std=c11`, warnings, `-O3 -DNDEBUG`,
  include paths, and OpenMP flags; all 34 programs link only the OpenMP
  runtime and `libSystem`. The scan of 62 C sources (tracked and generated)
  for banned constructs matched nothing.
- **Warnings:** none, from the library, the benchmark, or the tests, with
  either compiler.
- **Same work, same clock:** `src/ts1_bench.c` solves 64-cell groups by
  default (`src/ts1_bench.c:161`; `--group-size` or `CHEM_GROUP_SIZE` changes
  it), 157 groups for 10,000 cells, each with its own solver state, allocated
  before the timer (`:186-191`). It computes the rate constants once, in one
  untimed warm-up interval per group (`:196-202`), then times the 100
  intervals with `omp_get_wtime` (`:206-215`), one `omp parallel for` over the
  groups, static schedule unless `OMP_SCHEDULE` is set (`:195`). Its solver
  work at ten threads: 141,300 internal steps, all accepted, 282,600 function
  calls, 141,300 Jacobian updates and decompositions, and 423,900 solves, 9
  steps per group and interval, as MICM takes in each of its groups (9,000
  steps for its default 10 groups); the agent's own run of MICM with 64-cell
  groups reported the same counts. MICM's timed programs use one group per
  thread and 32-cell groups, so the session compares terra-max's 64-cell
  grouping against both.
- **Self-report** (one line, not a timing): its RESULTS.md gives 9,632 ms
  against 23,835 ms for MICM with 64-cell groups, in one run each at 10
  threads.

## Timing

Timed on 2026-11-17, 20:49 to 23:39, with every other entry and MICM's three
programs, in
[`2026-11-17-tier1-all-twelve-timing`](../2026-11-17-tier1-all-twelve-timing/README.md):
5 rounds at 10 and 1 threads, every run passing the gate.

| | 10 threads, median ms (range) | vs MICM's default | vs MICM's fastest | 1 thread, median ms (range) | vs MICM's default | vs MICM's fastest |
|---|---:|---:|---:|---:|---:|---:|
| terra-max | 8,668 (8,591–8,707) | 3.15× | 1.11× | 43,691 (43,676–43,762) | 3.35× | 1.17× |

It is sixth of the twelve entries at ten threads and third on one, and gains
5.04× from one thread to ten. With its 64-cell groups it is 2.63× faster than
MICM with 32-cell groups at ten threads and 3.11× on one; MICM with 64-cell
groups was not timed in the session. Its self-report, 9,632 ms against 23,835 ms
for MICM with 64-cell groups in one ten-thread run, compares with 8,668 ms here.

## Compliance

`admin/scan_commands.py` on the session record: 600 commands, none reaching
outside the run directory and its inputs, none naming another run or the
experiment's private folders, no installs or downloads, no connector, MCP, or
web-search calls, and no subagents (no subagent record names the session).
All 167 file edits went to 55 files inside the run directory.
