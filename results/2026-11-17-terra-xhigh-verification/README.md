# Verification of terra-xhigh, 2026-11-17

terra-xhigh, the second Terra entry (gpt-6.1-terra at xhigh effort, Codex CLI
0.169.0), was rebuilt from its branch, retested, and audited with both
compilers on 2026-11-17 at 04:40, ten minutes after the run ended. It is timed
later, with the other two Terra entries, luna-high, and every other entry, in
one session; its times will be added here then.

| File | Contents |
|---|---|
| `tests.md` | `make test-cpu` results for terra-xhigh with each compiler |

## Builds, tests, and audit

terra-xhigh (`288f274`, 37 commits) was built with `make env`, `make cpu bench
OPTFLAGS="-O3 -DNDEBUG"`, and `make test-cpu`, once with `CC=clang` (Apple
clang 21.0.0) and once with `CC=gcc-15` (GCC 15.2.0), by
`admin/verify_tier3.py --stages build,tests,ts1,audit`, the method of
[`2026-11-17-terra-max-verification`](../2026-11-17-terra-max-verification/README.md).

- **Tests** (`tests.md`): 676 PASS, 0 FAIL, 4 N/A with each compiler, every id
  exactly once, as the agent reported; the N/A cases are terra-max's four, MICM's
  checks of C++ template concepts. `tools/test_cpu.py` runs 25 C test programs,
  each printing `<id> 1|0` for its cases (a repeated id is an error), and
  compares the eight tutorial programs' printed output with MICM's, timings and
  interleaved thread messages masked. The checks are per case:
  `tests/test_lu_cases.c:131-158`, for example, factors each LU case's matrix
  family at every ordering and grouped width.
- **TS1 gate:** passes with both compilers; largest relative difference from
  MICM 1.77e-13 with clang (BCARY, cell 6848) and 2.33e-13 with gcc-15
  (BCARY, cell 4179).
- **Flags, libraries, sources:** only `-std=c11`, warnings, `-O3 -DNDEBUG`,
  include paths, and OpenMP flags. The benchmark, the one program in `bin/`,
  links only the OpenMP runtime and `libSystem`; the 29 test programs in
  `build/`, outside the verifier's audit of `bin/`, were checked by hand and
  link only `libSystem` (the tutorial program also the OpenMP runtime). The
  scan of 56 C sources for banned constructs matched nothing.
- **Warnings:** none, with either compiler.
- **Same work, same clock:** `src/ts1_bench.c` solves each cell on its own,
  with its own step-size control, scalar per cell, one workspace per thread
  (`src/ts1_bench.c:176-198`), an `omp for` over the 10,000 cells, static
  schedule unless `OMP_SCHEDULE` is set (`:174`). It computes the rate
  constants while reading the input, before the timer (`:117`), runs one
  untimed warm-up interval per cell (`:180-184`), and times the 100 intervals
  with `omp_get_wtime` (`:187-198`). Its solver work at ten threads: 9,000,000
  internal steps, all accepted, 18,000,000 function calls, 9,000,000 Jacobian
  updates and decompositions, and 27,000,000 solves: MICM's 9 steps per
  interval, for each cell rather than each group.
- **Self-report** (one line, not a timing): its RESULTS.md gives 30,496 ms
  against 29,617 ms for MICM's default, one run each at 10 threads.

## Compliance

`admin/scan_commands.py` on the session record: 647 commands, none reaching
outside the run directory and its inputs, none naming another run or the
experiment's private folders, no installs or downloads, no connector, MCP, or
web-search calls, and no subagents (no subagent record names the session).
All 137 file edits went to 48 files inside the run directory.
