# Fable's reference entry: timing sessions and gate matrix, 2026-10-17

The entry in `fable/` on `main` is Claude Fable's own solution to the Tier 1
prompt, built outside the contest at the SSE's request (wiki: Fable's
reference entry). These are the session's measurements, made with
`admin/time_entries.py` and `fable/tools/check_benchmark.py` on the run
machine (Apple M5, 4 performance and 6 efficiency cores, macOS 27.0.1, Apple
clang 21, `-O3 -DNDEBUG`), on AC power with the machine otherwise idle after
a cool-down.

- `session-1/`: 05:01 to 05:37, the build before the minimum-fill ordering
  (commit `842e99d`), whose default grouping was 8 cells, timed as `fable`,
  `fable-g16`, and `fable-g32` (its one-thread default was already 32 cells),
  with MICM's default configuration and its fastest (vector ordering of 8,
  in-place Mozart, 32-cell groups). 5 rounds at 10 and 1 threads; all 50 runs
  within 1.1e-13 of MICM.
- `session-2/`: 07:02 to 07:37, the final build (commit `ec2b170`, minimum
  fill, 32-cell groups by default), the same programs; `fable` and
  `fable-g32` are the same configuration timed twice. 5 rounds; all 50 runs
  within 1.04e-13 of MICM.
- `session-3/`: 08:02 to 08:37, the final commit (`5e770a9`, which adds the
  second stage's forcing formed block by block; the kernels are otherwise
  session 2's), the same programs. 5 rounds; all 50 runs within 1.04e-13.
- `gate-matrix/`: the build of session 2 in every benchmark setting (groupings,
  one thread, the generic kernels, the schedules), each checked against
  MICM's final state: every setting passes, 1.04e-13 at most, the generic
  kernels exact.

The input is 10,000 cells from `reference/make_cells_csv.py 10000 2026`,
100 steps of 30 s after one warm-up step, checked with
`reference/compare_final.py` against MICM's default configuration.

| Program | 10 threads, median ms | vs MICM's fastest | 1 thread, median ms | vs MICM's fastest |
|---|---:|---:|---:|---:|
| MICM default (session 3) | 26,824 | | 145,280 | |
| MICM fastest (session 3) | 9,892 | 1.00× | 50,900 | 1.00× |
| fable, final commit (session 3) | 4,594 | 2.15× | 25,558 | 1.99× |
| the same, timed as `fable-g32` (session 3) | 4,558 | 2.17× | 25,536 | 1.99× |
| fable, minimum fill (session 2) | 4,757 | 2.05× (MICM's fastest: 9,731) | 25,786 | 1.97× (50,735) |
| fable, before minimum fill, 32-cell groups (session 1) | 4,850 | 2.02× (MICM's fastest: 9,782) | 26,364 | 1.93× (50,792) |
