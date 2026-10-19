# Results

Tier 3's results come first; the Tier 1 measurements follow under "Tier 1".

# Tier 3

The operating points of OPERATING_POINTS.md on the six disclosed scenarios,
against MICM's three-stage Rosenbrock, all at 1,000 cells on the Apple M5
laptop (10 cores) with Apple clang 21, built with `-O3 -DNDEBUG`. Accuracy
is `scenarios/score.py`'s: digits overall and for the worst species, with
the gate (no NaN or infinity, nothing below -1e-20 mol m-3). Times are
medians of three runs at 10 threads in one session (scan S13, 06:27 to 06:39
on 2026-10-19, idle machine, AC power; `scans/S13-confirmation.csv`), and
one run at one thread (S14). The experimenter's ranking times 10,000-cell
scenarios, median of five, in its own session; those numbers are the
ranking's, not these.

## Operating points

Every run passes the gate, and every point's output is byte-identical
across the three rounds, so the digits are exact. Per scenario: digits
overall (worst species and its digits), median elapsed ms at 10 threads.

| Point | clean-10min | diurnal-10min | diurnal-30min | diurnal-5min | polluted-10min | tier1-30s | Lowest (worst) | Level | Time, s |
|---|---|---|---|---|---|---|---|---|---:|
| L2 | 2.58 (ISOP 1.52) 175 | 2.56 (ISOP 1.39) 205 | 2.52 (ISOP 1.37) 162 | 2.73 (NO 1.87) 281 | 2.42 (ISOP 1.29) 247 | 5.09 (C6H5O2 4.40) 89 | 2.42 (1.29) | 2 | 1.160 |
| balanced | 3.26 (ISOP 2.26) 365 | 3.25 (ISOP 2.17) 410 | 3.30 (ISOP 2.32) 335 | 3.75 (ISOP 2.89) 488 | 3.17 (ISOP 2.17) 482 | 5.80 (ISOP 5.02) 128 | 3.17 (2.17) | 3 (declared 2) | 2.208 |
| L3 | 3.61 (ISOP 2.49) 399 | 3.84 (ISOP 2.58) 449 | 3.65 (ISOP 2.52) 387 | 4.12 (ISOP 3.00) 644 | 3.95 (ISOP 2.85) 550 | 5.87 (ISOPO2 5.32) 217 | 3.61 (2.49) | 3 | 2.647 |
| recommended | 3.62 (MALO2 2.87) 467 | 3.65 (OClO 2.69) 471 | 4.24 (ISOP 3.10) 372 | 3.87 (BrO 3.00) 830 | 4.52 (ISOP 3.40) 599 | 6.53 (MEKO2 5.86) 303 | 3.62 (2.69) | 3 | 3.041 |
| accurate | 4.81 (ISOP 3.66) 530 | 4.72 (ISOP 3.49) 556 | 4.57 (ISOP 3.54) 488 | 3.87 (BrO 3.00) 870 | 4.54 (ISOP 3.41) 701 | 6.84 (MEKO2 6.30) 303 | 3.87 (3.00) | 3 | 3.449 |

Operation counts per interval and 32-cell group, the mean over the six
scenarios (the per-scenario values are in the CSV): internal steps,
rejections, forcing evaluations, Jacobians, factorizations, solves.

| Point | Steps | Rejected | Forcings | Jacobians | Factorizations | Solves |
|---|---:|---:|---:|---:|---:|---:|
| L2 | 2.12 | 0.21 | 6.2 | 2.1 | 2.1 | 8.5 |
| balanced | 4.16 | 0.35 | 12.1 | 4.2 | 4.2 | 16.6 |
| L3 | 2.11 | 0.10 | 14.7 | 2.0 | 6.3 | 19.0 |
| recommended | 1.62 | 0.06 | 17.7 | 1.6 | 6.5 | 22.6 |
| accurate | 1.91 | 0.11 | 21.0 | 1.8 | 7.7 | 26.8 |

Photolysis rate constants are recomputed in 0.5 of the intervals (none by
night, when the rates stay zero). The 30-minute scenario takes two to three
times the steps per interval of the 10-minute ones; the 5-minute scenario
and Tier 1's 30-second steps take about one step per interval for every
point, which is why the extrapolation points, at 9 and 14 substeps per step,
are dearest there (644 to 870 ms on diurnal-5min against 281 for L2).

## Against MICM's three-stage Rosenbrock

*(MICM at the ranking's baseline tolerances on the same six scenarios in the
same session, `scans/MICM-baselines.csv`, and each point's speedup: filled
when the rerun ends.)*

## Work-precision plots

*(figures in `plots/`, from `tools/plot_wp.py`.)*

## Where the time goes

The five points at one thread (S14, one run per scenario,
`scans/S14-finals-1thread.csv`), the time summed over the six scenarios,
the speedup at 10 threads from S10's generated runs in the same session,
and the kernel cost model's share: the model multiplies the operation
counts by the per-cell kernel costs measured with `bin/kernel_bench` at
32-cell groups (forcing 225 ns, Jacobian 438 ns, factorization 676 ns,
triangular solve 279 ns) and compares with the measured nanoseconds per
cell and interval.

| Point | 1 thread, s | 10 threads, s | Speedup | Model / measured, range over the scenarios |
|---|---:|---:|---:|---|
| L2 | 5.96 | 1.17 | 5.1 | 0.70 to 0.79 |
| balanced | 10.96 | 2.24 | 4.9 | 0.75 to 0.78 |
| L3 | 12.78 | 2.76 | 4.6 | 0.76 to 0.78 |
| recommended | 15.02 | 3.12 | 4.8 | 0.74 to 0.78 |
| accurate | 17.11 | 3.53 | 4.8 | 0.75 to 0.79 |

- The four kernels account for about three quarters of the time at every
  point; the rest is the loop's own work (the error norm and the step
  control, the stage combinations, the extrapolation table, the rate
  constant updates, and copying the states in and out of the vector
  layout). The share is the same for the Rosenbrock and the extrapolation
  loops, so neither has an overhead worth attacking before the kernels.
- Ten threads buy 4.6 to 5.1 times on the M5's four performance and six
  efficiency cores, the same ratio Tier 1's benchmark showed.
- Per cell and interval at one thread, L2 costs 5.0 to 8.4 microseconds on
  the 10-minute and 5-minute scenarios and 18 on the 30-minute one; the
  recommended point 14 to 21 and 39.

## Generated against generic kernels

S10 (`scans/S10-generic.csv`): each point once on the six scenarios with
the generated kernels and with the generic (index-driven) kernels of the
same library, 10 threads, in one session.

| Point | Generated, s | Generic, s | Ratio |
|---|---:|---:|---:|
| L2 | 1.17 | 6.46 | 5.5 |
| balanced | 2.24 | 12.55 | 5.6 |
| L3 | 2.76 | 14.70 | 5.3 |
| recommended | 3.12 | 18.28 | 5.9 |
| accurate | 3.53 | 19.53 | 5.5 |

The two paths take exactly the same steps (identical internal step counts
on every scenario) and reach the same digits to the third decimal; the
outputs differ at roundoff (the generated kernels sum in a different
order), so their SHA-256 differ. The generic path is the fallback for a
mechanism without generated kernels, as in Tier 1; the Tier 3 loops do not
change the ratio between them.

## Build and run times

*(make env, make bench, make test-cpu with clang and gcc-15.)*

# Tier 1

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

**The final build** (minimum-fill ordering, 32-cell groups, the second
stage's forcing formed block by block), 2026-10-17 08:02 to 08:37; all 50
runs within 1.04e-13 of MICM:

| Program | 10 threads, median ms (range) | vs MICM | vs MICM's fastest | 1 thread, median ms (range) | vs MICM | vs MICM's fastest |
|---|---:|---:|---:|---:|---:|---:|
| MICM, default | 26,824 (25,890–27,505) | 1.00× | | 145,280 (144,403–147,173) | 1.00× | |
| MICM's fastest: vector8, in-place Mozart, 32-cell groups | 9,892 (9,285–9,949) | 2.71× | 1.00× | 50,900 (50,547–51,059) | 2.85× | 1.00× |
| this port, default (32-cell groups) | 4,594 (4,541–4,798) | 5.84× | 2.15× | 25,558 (25,452–25,589) | 5.68× | 1.99× |
| the same configuration, timed again as `--group-size 32` | 4,558 (4,551–4,701) | 5.89× | 2.17× | 25,536 (25,412–25,569) | 5.69× | 1.99× |
| this port, 16-cell groups | 4,776 (4,665–4,898) | 5.62× | 2.07× | 36,877 (34,961–37,289) | 3.94× | 1.38× |

The two earlier sessions, each five rounds against the same MICM programs:
the minimum-fill build without the block-by-block forcing (07:02 to 07:37)
gave 4,757 ms at 10 threads against 9,731 (2.05×) and 25,786 ms at one
thread against 50,735 (1.97×); the build before the minimum-fill ordering,
with MICM's order (05:01 to 05:37), 4,850 against 9,782 (2.02×) and 26,364
against 50,792 (1.93×), its 16-cell groups 33,324 and its 8-cell groups
5,030 at 10 threads (1.94×).

Against MICM's fastest configuration the port is 2.15× at 10 threads, the
ranking configuration, and 1.99× at 1 thread. The 10-thread speedups of the
two codes over their own 1-thread times are close (5.1× and 5.6×), so both
are bound by the same thing there: the mixed cores' throughput.

### The one-thread target

The one-thread ratio is 1.99× in the final session (1.97× and 1.93× in the
two before it, as the code and the machine's state changed); 2× would need
another 0.4 %, less than the sessions' own spread. Everything identified
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
  the one small fusion kept;
- larger kernel chunks (`KERNEL_CHUNK` 64 and 128, fewer calls per block):
  1.3 % fewer instructions, the same time at one thread, 2 to 3 % slower at
  ten.

Not yet swept, for a follow-up session: the fast mode's split-sum threshold
(`SPLIT_MIN`, 12 terms) and group sizes between 32 and 64.

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
