# Results

Measured on the experiment's machine: Apple M5 (4 performance and 6
efficiency cores), 16 GB, macOS 27.0.1, Apple clang 21.0.0 and GCC 15.2.0,
`OPTFLAGS=-O3 -DNDEBUG` for the port and for MICM. Times are the programs'
own `elapsed_ms`; the experimenter's timing session is the ranking.

## Test cases

`make test-cpu` prints one line per id of `reference/micm_cases.txt` (680
ids, each once; the sets were diffed), then the port's two extra checks on
stderr. Neither compiler warns.

| Compiler | PASS | FAIL | N/A |
|---|---:|---:|---:|
| Apple clang 21 | 673 | 0 | 7 |
| GCC 15 | 673 | 0 | 7 |

The seven N/A cases are MICM tests of C++ types alone, with nothing that
runs: static assertions on view-category tags and template policies
(`view_category_sfinae`, four cases; `grouped_views`, one) and which C++
overload or template argument a builder passes to a finalize callback
(`solver_builder`, two). The extra checks: the generated TS1 kernels against
the generic kernels on the 100-cell reference input (largest relative
difference 3.5e-15), and the experiment's gate on that input (3.8e-14 in
fast mode with clang, 3.1e-14 with GCC; 0 in exact mode).

## TS1 gate

`tools/check_benchmark.py` on 10,000 cells (seed 2026), 100 steps of 30 s,
against MICM's final state; the gate is |y - y_MICM| <= 1e-30 + 1e-10 |y_MICM|.

| Setting | Arguments | Gate | Largest relative difference | Internal steps |
|---|---|---|---:|---:|
| default (all threads, groups of 32) | | PASS | 1.1e-13 | 281,700 |
| 1 thread | `--threads 1` | PASS | 1.1e-13 | 281,700 |
| groups of 4 | `--group-size 4` | PASS | 1.1e-13 | 2,250,000 |
| groups of 8 | `--group-size 8` | PASS | 1.1e-13 | 1,125,000 |
| groups of 16 | `--group-size 16` | PASS | 1.1e-13 | 562,500 |
| groups of 64 | `--group-size 64` | PASS | 1.1e-13 | 141,300 |
| generic kernels | `--generic` | PASS | 0 | 281,700 |
| static schedule | `--schedule static` | PASS | 1.1e-13 | 281,700 |
| guided schedule | `--schedule guided` | PASS | 1.1e-13 | 281,700 |
| 2 threads, groups of 3 | `--threads 2 --group-size 3` | PASS | 1.1e-13 | 3,000,600 |

Every grouping takes 9 internal steps per group per 30 s step on this input,
as MICM does (MICM's 32-cell groups: 281,700 internal steps). The generic
kernels reproduce MICM exactly; the generated kernels in fast mode differ by
the reassociated sums and the reciprocal multiply, 1.1e-13 at most. Exact
mode (`make KERNEL_MODE=exact`) differs by 0 on the 100-cell input.

## TS1 timing

`admin/time_entries.py`, 5 rounds in rotating order at 10 and then 1
threads, 10,000 cells (seed 2026), 100 steps of 30 s, every run checked
against MICM's final state (all 50 passed, largest relative difference
1.1e-13). MICM's programs are `reference/micm_ts1_ref.cpp` as built by
`admin/micm_configs.sh`: the default (standard ordering, in-place Mozart,
one cell group per thread) and the fastest configuration found on
2026-10-16 (vector ordering of 8, in-place Mozart, 32-cell groups).

| Program | 10 threads, median ms (range) | vs MICM | vs MICM's fastest | 1 thread, median ms (range) | vs MICM | vs MICM's fastest |
|---|---:|---:|---:|---:|---:|---:|
| MICM, default | 26,979 (26,652–27,303) | 1.00× | | 145,166 (144,632–146,272) | 1.00× | |
| MICM's fastest: vector8, in-place Mozart, 32-cell groups | 9,782 (9,233–9,866) | 2.76× | 1.00× | 50,792 (50,354–51,024) | 2.86× | 1.00× |
| this port, 32-cell groups (the default) | 4,850 (4,816–4,916) | 5.56× | 2.02× | 26,364 (26,228–26,465) | 5.51× | 1.93× |
| this port, 16-cell groups | 4,917 (4,898–5,055) | 5.49× | 1.99× | 33,324 (30,559–34,126) | 4.36× | 1.52× |
| this port, 8-cell groups | 5,030 (4,957–5,129) | 5.36× | 1.94× | 26,316 (26,146–26,404) * | | |

\* The session timed the build before 32 became the default, whose 1-thread
default was already 32-cell groups, so its "8-cell" row at 1 thread is a
second 32-cell measurement; back to back on the same machine, 8-cell groups
take 6.9 s per 20 steps at 1 thread against 5.2 s for 32-cell groups.

Against MICM's fastest configuration the port is 2.0× at 10 threads, the
ranking configuration, and 1.9× at 1 thread. The 10-thread speedups of the
two codes over their own 1-thread times are the same (5.2×), so both are
bound by the same thing there: the mixed cores' throughput.

### The one-thread target

The one-thread ratio is 1.93× in the session above and 1.97× on a cooler
machine; 2× would need 3.5 %. Everything identified was measured, back to
back, on the final code (DECISIONS.md, "Measured rejects"):

- the LU panels finished target by target and 32-cell groups were the two
  gains, 9 % and 24 % at one thread, and are in;
- forming each solve's right-hand side inside the kernel instead of a
  separate pass and a copy: bit-identical, 1.5 % slower (76 KB more code);
- the solution, error, and norm in one pass: bit-identical, 3 % more
  instructions, no faster;
- removing the compiler's loop prefetches: 2.8 % fewer instructions, 3 %
  slower;
- the error norm's division costs nothing measurable (a build without it
  runs in the same time);
- a greedy minimum-fill ordering would save 5 % of the elimination
  multiply-adds and 3 % of the solve's, about 2 % of the run, and give up
  bit-for-bit comparison with MICM.

What remains is the formulation's floor: 47,700 instructions per cell and
internal step at 3.75 per cycle, with the LU, solves, Jacobian, and forcing
each at one to two instructions per multiply-add, the same operation
counts as MICM's method, and code that no longer fits the instruction cache
however it is grouped. The remaining ideas (two-row register tiles in the
trailing update, the pivot reciprocals stored in the diagonal) are each
worth about 1 % of the instructions, less than the machine's run-to-run
spread, so the one-thread target was not reached and is judged out of reach
for this design.

## Where the time goes

Instructions per cell and internal step (`/usr/bin/time -l`, one thread,
32-cell groups, differenced between 20- and 40-step runs so setup cancels):
47,700, at 3.75 instructions per cycle. Static counts of the generated
kernels per block of two cells: LU 24,500 (11,300 loads, 2,900 stores),
Jacobian 15,600, forcing 8,400 (run twice per step), solve 9,400 (run three
times): 85,200 per internal step, 42,600 per cell; the rest is the
library's stage combinations, error norm, and clamp. A sampling profile of
the same run: LU 27 %, the three solves 30 %, the Jacobian 18 %, the two
forcings 14 %, the library's vector passes and step control 12 %.

The kernels' code (496 KB) exceeds the 192 KB instruction cache, and a
group's blocks run in lock step through each loop body, so the group size
sets how often that code is fetched: at one thread, groups of 4 cells run
at 2.4 instructions per cycle, 8 at 2.8, and 32 at 3.7. After hours of
heavy load this laptop runs the kernels 30-50 % slower at the same clock
(MICM's loops and the generic kernels are unaffected), so the session above
was run on a cool, idle machine.

Generator variants, instructions per cell and internal step at 8-cell
groups (the gate holds for each; the panel size changes no arithmetic, so
its results are byte for byte identical):

| Variant | Instructions | Change | Kernel code (KB) |
|---|---:|---:|---:|
| panel 4 (rows and columns of the panel updated statement by statement) | 51,532 | | 513 |
| `--panel 16`, same scheme | 50,682 | -1.6 % | 503 |
| `--panel 64`, same scheme | 50,154 | -2.7 % | 500 |
| `--panel 211` (one panel), same scheme | 52,403 | +1.7 % | 515 |
| panel 64, every element finished target by target (the default) | 48,387 | -6.1 % | 484 |
| `--interleave 2` (with 4-pivot panels) | 51,768 | +0.5 % | 513 |
| `--interleave 8` (with 4-pivot panels) | 52,628 | +2.1 % | 516 |
| `--row-budget 32` (with 4-pivot panels) | 52,017 | +0.9 % | 515 |
| `--row-budget 128` (with 4-pivot panels) | 51,818 | +0.6 % | 513 |

Back to back at one thread with 32-cell groups, the target-by-target panels
took 5.23 s per 20 steps against 5.46 s for the statement-by-statement
64-pivot panels and 5.74 s for 4-pivot panels. Block offsets the compiler
cannot analyze remove its loop prefetches and 2.8 % of the instructions but
cost 3 % of the time, so the prefetches stay.
