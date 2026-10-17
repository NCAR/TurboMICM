# Results

Measured on the experiment's machine: Apple M5 (4 performance and 6
efficiency cores), 16 GB, macOS 27.0.1, Apple clang 21.0.0 and GCC 15.2.0,
`OPTFLAGS=-O3 -DNDEBUG` for the port and for MICM. Times are the programs'
own `elapsed_ms`; the experimenter's timing session is the ranking.

## Test cases

`make test-cpu` prints one line per id of `reference/micm_cases.txt` (680
ids, each once), then the port's two extra checks on stderr.

| Compiler | PASS | FAIL | N/A |
|---|---:|---:|---:|
| Apple clang 21 | (final) | | |
| GCC 15 | (final) | | |

The N/A cases are MICM tests of C++ types alone (static assertions on view
categories and template policies, the C++ overload a builder picks), with
nothing that runs. The extra checks: the generated TS1 kernels against the
generic kernels on the 100-cell reference input (largest relative difference
3.5e-15), and the experiment's gate on that input (3.8e-14 in fast mode, 0
in exact mode).

## TS1 gate

`tools/check_benchmark.py` on 10,000 cells (seed 2026), 100 steps of 30 s,
against MICM's final state; the gate is |y - y_MICM| <= 1e-30 + 1e-10 |y_MICM|.

| Setting | Arguments | Gate | Largest relative difference | Internal steps |
|---|---|---|---:|---:|
| default (all threads, groups of 8) | | PASS | 1.1e-13 | 1,125,000 |
| 1 thread | `--threads 1` | PASS | 1.1e-13 | 1,125,000 |
| groups of 4 | `--group-size 4` | PASS | 1.1e-13 | 2,250,000 |
| groups of 16 | `--group-size 16` | PASS | 1.1e-13 | 562,500 |
| groups of 32 | `--group-size 32` | PASS | 1.1e-13 | 281,700 |
| generic kernels | `--generic` | PASS | 0 | 1,125,000 |
| static schedule | `--schedule static` | PASS | 1.1e-13 | 1,125,000 |
| guided schedule | `--schedule guided` | PASS | 1.1e-13 | 1,125,000 |
| 2 threads, groups of 3 | `--threads 2 --group-size 3` | PASS | 1.1e-13 | 3,000,600 |

Every grouping takes 9 internal steps per group per 30 s step on this input,
as MICM does (MICM's 32-cell groups: 281,700). The generic kernels reproduce
MICM exactly; the generated kernels in fast mode differ by the reassociated
sums and the reciprocal multiply, 1.1e-13 at most. Exact mode
(`make KERNEL_MODE=exact`) differs by 0 on the 100-cell input and on this one.

## TS1 timing

(final session on a quiet machine: 5 rotating rounds of
`admin/time_entries.py` at 10 and 1 threads, MICM's default and fastest
configurations beside this entry)

## Where the time goes

Instructions per cell and internal step (`/usr/bin/time -l`, one thread,
differenced between 20- and 40-step runs so setup cancels): 51,800 in the
default build. Static counts of the generated kernels per block of two
cells: LU 31,900 (12,600 loads, 5,800 stores), Jacobian 15,600, forcing
8,400 (run twice per step), solve 9,400 (run three times), 92,600 per
internal step, 46,300 per cell; the rest is the library's stage
combinations, error norm, and clamp. The code runs at 3.4 instructions per
cycle on a performance core with groups of 8 (4 blocks in lock step); the
526 KB of kernel code exceed the instruction cache, and groups of 4 (2 blocks
per loop body) fall to 1.2-2.4 instructions per cycle.

Generator variants, instructions per cell and internal step (the gate holds
for each):

| Variant | Instructions | Change | Kernel code (KB) |
|---|---:|---:|---:|
| panel 4, interleave 4, row budget 64 | 51,532 | | 513 |
| `--panel 8` | 51,253 | -0.5 % | 508 |
| `--panel 16` | 50,682 | -1.6 % | 503 |
| `--panel 32` | 50,414 | -2.2 % | 501 |
| `--panel 64` | 50,154 | -2.7 % | 500 |
| `--panel 211` (one panel) | 52,403 | +1.7 % | 515 |
| `--interleave 2` | 51,768 | +0.5 % | 513 |
| `--interleave 8` | 52,628 | +2.1 % | 516 |
| `--row-budget 32` | 52,017 | +0.9 % | 515 |
| `--row-budget 128` | 51,818 | +0.6 % | 513 |

The panel size changes no arithmetic (every element takes its updates in
pivot order), so the results are identical; larger panels load and store
each trailing element fewer times until the panel's own rows and columns,
updated pivot by pivot, dominate.
