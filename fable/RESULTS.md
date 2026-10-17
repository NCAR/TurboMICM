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
difference 5.9e-15), and the experiment's gate on that input (4.5e-14 in
fast mode with clang, 3.3e-14 with GCC; 0 in exact mode, which keeps MICM's
species ordering).

## TS1 gate

`tools/check_benchmark.py` on 10,000 cells (seed 2026), 100 steps of 30 s,
against MICM's final state; the gate is |y - y_MICM| <= 1e-30 + 1e-10 |y_MICM|.

| Setting | Arguments | Gate | Largest relative difference | Internal steps |
|---|---|---|---:|---:|
| default (all threads, groups of 32) | | PASS | 1.0e-13 | 281,700 |
| 1 thread | `--threads 1` | PASS | 1.0e-13 | 281,700 |
| groups of 4 | `--group-size 4` | PASS | 1.0e-13 | 2,250,000 |
| groups of 8 | `--group-size 8` | PASS | 1.0e-13 | 1,125,000 |
| groups of 16 | `--group-size 16` | PASS | 1.0e-13 | 562,500 |
| groups of 64 | `--group-size 64` | PASS | 1.0e-13 | 141,300 |
| generic kernels | `--generic` | PASS | 0 | 281,700 |
| static schedule | `--schedule static` | PASS | 1.0e-13 | 281,700 |
| guided schedule | `--schedule guided` | PASS | 1.0e-13 | 281,700 |
| 2 threads, groups of 3 | `--threads 2 --group-size 3` | PASS | 1.0e-13 | 3,000,600 |

Every grouping takes 9 internal steps per group per 30 s step on this input,
as MICM does (MICM's 32-cell groups: 281,700 internal steps). The generic
kernels reproduce MICM exactly; the generated kernels in fast mode differ by
the reassociated sums, the reciprocal multiply, and the minimum-fill pivot
order, 1.0e-13 at most. Exact mode (`make KERNEL_MODE=exact`) differs by 0
on the 100-cell input.

## TS1 timing

`admin/time_entries.py`, 5 rounds in rotating order at 10 and then 1
threads, 10,000 cells (seed 2026), 100 steps of 30 s, every run checked
against MICM's final state, on AC power with the machine otherwise idle after
a cool-down. MICM's programs are `reference/micm_ts1_ref.cpp` as built by
`admin/micm_configs.sh`: the default (standard ordering, in-place Mozart,
one cell group per thread) and the fastest configuration found on
2026-10-16 (vector ordering of 8, in-place Mozart, 32-cell groups). The data
are in `results/2026-10-17-fable-reference-entry`.

**The final build** (minimum-fill ordering, 32-cell groups), 2026-10-17
07:02 to 07:37; all 50 runs within 1.04e-13 of MICM:

| Program | 10 threads, median ms (range) | vs MICM | vs MICM's fastest | 1 thread, median ms (range) | vs MICM | vs MICM's fastest |
|---|---:|---:|---:|---:|---:|---:|
| MICM, default | 26,789 (25,628–27,801) | 1.00× | | 145,266 (144,453–145,969) | 1.00× | |
| MICM's fastest: vector8, in-place Mozart, 32-cell groups | 9,731 (9,410–9,987) | 2.75× | 1.00× | 50,735 (50,320–51,009) | 2.86× | 1.00× |
| this port, default (32-cell groups) | 4,757 (4,583–4,862) | 5.63× | 2.05× | 25,786 (25,663–25,850) | 5.63× | 1.97× |
| the same configuration, timed again as `--group-size 32` | 4,608 (4,592–4,757) | 5.81× | 2.11× | 25,719 (25,626–25,804) | 5.65× | 1.97× |
| this port, 16-cell groups | 4,788 (4,636–4,875) | 5.60× | 2.03× | 36,115 (32,728–37,613) | 4.02× | 1.40× |

**The build before the minimum-fill ordering** (MICM's ordering, 8-cell
groups by default), 05:01 to 05:37; all 50 runs within 1.1e-13:

| Program | 10 threads, median ms (range) | vs MICM's fastest | 1 thread, median ms (range) | vs MICM's fastest |
|---|---:|---:|---:|---:|
| MICM, default | 26,979 (26,652–27,303) | | 145,166 (144,632–146,272) | |
| MICM's fastest | 9,782 (9,233–9,866) | 1.00× | 50,792 (50,354–51,024) | 1.00× |
| this port, 32-cell groups | 4,850 (4,816–4,916) | 2.02× | 26,364 (26,228–26,465) | 1.93× |
| this port, 16-cell groups | 4,917 (4,898–5,055) | 1.99× | 33,324 (30,559–34,126) | 1.52× |
| this port, 8-cell groups | 5,030 (4,957–5,129) | 1.94× | (its 1-thread default was 32 cells) | |

Against MICM's fastest configuration the port is 2.05× at 10 threads, the
ranking configuration, and 1.97× at 1 thread. The 10-thread speedups of the
two codes over their own 1-thread times are the same (5.2–5.4×), so both are
bound by the same thing there: the mixed cores' throughput.

### The one-thread target

The one-thread ratio is 1.97× in the final session (1.93× before the
minimum-fill ordering); 2× would need another 1.5 %. Everything identified
was measured, back to back, on the final code (DECISIONS.md, "Measured
rejects"):

- the LU panels finished target by target, 32-cell groups, and the
  minimum-fill ordering were the gains, 9 %, 24 %, and 2.7 % at one thread,
  and are in;
- forming each solve's right-hand side inside the kernel instead of a
  separate pass and a copy: bit-identical, 1.5 % slower (76 KB more code);
- the solution, error, and norm in one pass: bit-identical, 3 % more
  instructions, no faster;
- removing the compiler's loop prefetches: 2.8 % fewer instructions, 3 %
  slower;
- skipping the Jacobian assembly's zero stores into the LU's fill positions
  (the LU starting them from zero): bit-identical, fewer instructions, 25 %
  slower; the pivot reciprocals kept in the diagonal: bit-identical, 5 %
  slower;
- the error norm's division costs nothing measurable (a build without it
  runs in the same time);
- minimum degree and weighted fill-plus-cost orderings leave more fill than
  greedy minimum fill;
- the LU's trailing update in row pairs sharing their loads of the panel's
  row values: bit-identical, 0.6 % more instructions, the same time;
- the second stage's forcing with its combination y + a k formed block by
  block in a thread-local buffer: bit-identical, 0.5 % faster at one thread,
  the one small fusion kept.

What remains is the formulation's floor: 45,900 instructions per cell and
internal step at 3.75 per cycle, with the LU, solves, Jacobian, and forcing
each at one to two instructions per multiply-add, the same operation counts
as MICM's method with a better pivot order, and code that no longer fits the
instruction cache however it is grouped. Every idea identified was
implemented and measured, the last being the row-pair tiles above; the
rejects show that fewer instructions do not reliably mean less time on this
core, and the machine's thermal state moves the one-thread ratio by about
2 %, more than the shortfall. The one-thread target was not reached and is
judged out of reach for this design.

## Where the time goes

Instructions per cell and internal step (`/usr/bin/time -l`, one thread,
32-cell groups, differenced between 20- and 40-step runs so setup cancels):
45,900 with the minimum-fill ordering, 47,700 with MICM's, at 3.75
instructions per cycle. Static counts of the generated kernels per block of
two cells with MICM's ordering: LU 24,500 (11,300 loads, 2,900 stores),
Jacobian 15,600, forcing 8,400 (run twice per step), solve 9,400 (run three
times): 85,200 per internal step, 42,600 per cell; the rest is the
library's stage combinations, error norm, and clamp. A sampling profile of
the same run: LU 27 %, the three solves 30 %, the Jacobian 18 %, the two
forcings 14 %, the library's vector passes and step control 12 %. The
minimum-fill ordering takes the LU from 2,381 to 2,293 elements, its
multiply-adds from 6,986 to 6,314 per cell and the solves' from 4,553 to
4,377.

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
