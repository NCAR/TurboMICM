# The Fable entry on Derecho

The experiment's reference C port of MICM (`../fable`, the Fable entry, tuned on
an Apple M5 with clang) built, verified, and tuned for NCAR's Derecho: CPU nodes
of two AMD EPYC 7763 (Milan, Zen 3) processors, 128 cores with two hardware
threads each, 8 NUMA domains, AVX2 and FMA, 512 KB of L2 per core and 32 MB of
L3 per 8-core complex, with the GCC that MUSICA builds with there (gcc 14.3.0
from ncarenv/25.10). The question: does the port stay 2x faster than MICM's
fastest configuration on this machine, at a whole node and at one thread?

Everything here lives on branch `derecho/fable`, under `derecho/`, so the other
branches' work on `fable/` and the run root is untouched.

## Layout

| Path | Contents |
|---|---|
| `fable/` | a copy of `../fable` at d7e9722 with Derecho defaults and tuning (below); its own DESIGN.md, DECISIONS.md, and RESULTS.md are the entry's |
| `env.sh` | the modules (ncarenv/25.10, gcc/14.3.0), the PBS account, and OpenMP binding; `source derecho/env.sh` |
| `tools/submit.sh` | `NAME NCPUS WALLTIME COMMAND...`: runs a command as a develop-queue (cpudev) job from `derecho/`, log in `logs/NAME.log` |
| `tools/micm_reference.sh` | builds MICM's TS1 reference driver with g++ and writes the gate's reference final state for an input |
| `tools/micm_configs.sh` | the bash and g++ port of `admin/micm_configs.sh`: MICM's 24 CPU configurations built, checked, and timed |
| `tools/build_variant.sh` | `NAME [MAKE_VARIABLE=value ...]`: builds a variant of `fable/` in `work/variants/NAME` |
| `tools/verify_port.sh` | build, the 680 ported test cases, and the gate at several run-time settings |
| `tools/environment.sh` | records the node, CPU, kernel, compiler, and modules of a timing session |
| `pbs/` | the timing sessions' scripts, one per session |
| `results/<date>-<topic>/` | each session's `timings.csv`, `summary.md`, environment, and notes |
| `micm/`, `work/`, `logs/` | not committed: MICM at e15898ce, inputs, builds, variants, job logs |

## Building and running

```
source derecho/env.sh
cd derecho/fable && make env && make cpu          # gcc 14, width 4, Derecho flags
make test-cpu                                     # the 680 cases, from MICM's source tree in ../micm
cd .. && tools/submit.sh name 128 01:00:00 COMMAND  # a whole dev node for COMMAND
```

`fable/Makefile`'s Derecho defaults: `CC=gcc`, `CXX=g++`,
`OPTFLAGS=-O3 -DNDEBUG -march=znver3 -ffp-contract=fast`, `WIDTH=4`,
`MICM=../micm`, `REFERENCE=../../reference`. The entry's rules allow only `CC`
and `OPTFLAGS` as compiler choices, and these are them.

## What Derecho needed

- **FMA contraction.** The kernels are written for fused multiply-adds, which
  clang contracts by default. GCC in `-std=c11` mode does not
  (`-ffp-contract=off` is its standards-mode default), so `OPTFLAGS` adds
  `-ffp-contract=fast`; `-march=znver3` enables AVX2 and FMA at all.
- **Vector width 4.** The entry's default width of 2 matches NEON's two
  doubles; AVX2 holds four, so `WIDTH=4`.
- **POSIX in the benchmark driver.** `bench/ts1_bench.c` uses `getline` and
  `ssize_t`, which glibc hides under `-std=c11`; the driver now defines
  `_POSIX_C_SOURCE 200809L` first. (The entry built on macOS with GCC 15
  because Apple's libc exposes them regardless.)
- **Linking MICM's C++ reference.** ncarenv/25.10 leaves `LIBRARY_PATH` at
  gcc 12.5.0's lib64 even with gcc/14.3.0 loaded, so g++ 14 links against the
  older libstdc++ and fails (`_M_replace_cold`, `__cxa_call_terminate`);
  `env.sh` unsets it.
- **PBS.** The dev nodes advertise their 128 cores as `ncpus` (hardware
  threads are not scheduled separately), so `ncpus=128` takes a whole node,
  and the job may use all 256 hardware threads; a user may have 256 `ncpus`
  running at once, two whole nodes. Jobs started within a minute throughout.
  Output streams only with `-k oed`. The default memory limit is 10 GB, which
  killed 24 parallel g++ compiles of MICM's reference (`MEM=200GB` in
  `submit.sh`).
- **MICM's threads and libgomp.** The reference driver threads with
  `std::thread`. Built with `-fopenmp` (which the experiment's build never
  used), it links libgomp, and libgomp pins the main thread to one core at
  startup when `OMP_PROC_BIND` is set, as `env.sh` sets it for the port; every
  `std::thread` inherits that mask. MICM then ran its 128 threads on one core
  (247 s against 363 s on one thread; 2, 4, and 8 threads all took the same
  time). MICM is built without `-fopenmp`.
- **Conda.** `make env` uses conda from PATH or Derecho's system miniforge;
  the conda module is not loaded.

## The generated kernels on Zen 3

The generator writes each lane of a block as a loop `for (c = 0; c < W; ++c)`
over W doubles. GCC unrolls these before vectorizing and then rebuilds vectors
with its SLP vectorizer, with poor results: at width 4 the TS1 kernels were
151,000 instructions, with 4-double loads split into two halves and an insert,
3,600 scalar operations, and 41,000 stack references in the 64-pivot LU
panels. Emitting each lane loop as `_Pragma("omp simd")` (OpenMP is part of
the entry's language) hands them to the loop vectorizer whole. Static counts of
`ts1.o` (`tools/asm_stats.sh`), gcc 14.3.0, `-O3 -march=znver3
-ffp-contract=fast`, every build within 4e-14 of MICM on the 100-cell input:

| Build | Instructions | ymm FMA | xmm FMA | ymm mul/add/sub | Scalar | Stack refs | Shuffles | Text (KB) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| width 4, plain lane loops | 151,067 | 1,564 | 0 | 6,773 | 4,005 | 41,448 | 3,223 | 1,097 |
| width 4, omp simd lanes (the default) | 75,348 | 9,597 | 0 | 7,632 | 1 | 6,338 | 0 | 496 |
| width 2, omp simd | 76,818 | 1 | 9,597 | 0 (7,632 xmm) | 1 | 7,383 | 0 | 506 |
| width 8, omp simd | 142,466 | 19,193 | 0 | 15,306 | 1 | 32,695 | 0 | 980 |
| width 4, omp simd, no `-ffp-contract=fast` | 84,294 | 0 | 0 | 26,641 | 1 | 7,045 | 0 | 529 |
| width 4, omp simd, plain `-O3` (SSE2) | 189,606 | 0 | 0 | 0 (53,322 xmm) | 1 | 40,682 | 0 | 1,180 |

Width 4 does twice width 2's work per instruction at the same count; width 8
doubles the code and spills five times as much (16 vector registers). GCC's
Zen tuning (`avoid_fma_chains`) still writes 7,632 of the width-4 build's
multiply-adds as separate multiplies and adds; session 2 measures
`-mtune-ctrl=^avoid_fma_chains,^avoid_fma256_chains`.

## Results

All timings are from whole dev nodes (`ncpus=128`, nothing else on the node),
`admin/time_entries.py` with rotating rounds, 10,000 cells of seed 2026, 100
steps of 30 s, every run gated against MICM's own Derecho result for the input
(`tools/micm_reference.sh`; MICM's result does not depend on its thread count
or grouping, and differs from its Mac result by 5e-14). Threads: 128 is one
OpenMP thread per core (`OMP_PLACES=cores`, `OMP_PROC_BIND=close`); MICM's
`std::thread`s are unpinned. MICM's programs here are its default (standard
ordering, in-place Mozart LU, one cell group per thread) and the two leading
configurations of the Mac sweep: vector ordering of 8 or 4 with 32-cell groups.
Run-to-run spread: the port's runs repeat within 1 % at both thread counts;
MICM's 128-thread runs vary 3–10 % and its one-thread runs 1 %; an identical
port binary run under two names differed by 3 % at one thread.

### Session 1: width and flags (`results/2026-10-18-session-1-width-flags`)

Node dec2443, 3 rounds. The Derecho defaults are width 4, `-O3 -DNDEBUG
-march=znver3 -ffp-contract=fast`, omp simd lane loops, 64-byte aligned arrays,
32-cell groups.

| Program | 128 threads, median ms (range) | vs MICM vector8 | 1 thread, median ms (range) | vs MICM vector8 |
|---|---:|---:|---:|---:|
| MICM default | 3,399 (3,148–3,638) | 0.58x | 366,157 (363,028–367,767) | 0.39x |
| MICM vector8, 32-cell groups | 1,968 (1,872–1,986) | 1.00x | 141,656 (141,496–142,548) | 1.00x |
| MICM vector4, 32-cell groups | 2,237 (2,180–2,267) | 0.88x | 162,960 (162,561–162,993) | 0.87x |
| **port, defaults (width 4)** | **589 (587–592)** | **3.34x** | **45,176 (45,130–45,365)** | **3.14x** |
| port, width 2 | 808 (807–811) | 2.43x | 65,629 (65,597–65,943) | 2.16x |
| port, width 8 | 730 (729–731) | 2.70x | 54,882 (53,897–55,124) | 2.58x |
| port, width 4, no `-ffp-contract=fast` | 624 (620–625) | 3.15x | 47,753 (47,687–47,842) | 2.97x |
| port, width 4, plain `-O3` (generic x86-64) | 1,118 (1,114–1,121) | 1.76x | 93,567 (92,743–96,188) | 1.51x |

Width 4 is 1.37x width 2 at 128 threads and 1.45x at one; width 8 spills.
FMA contraction is worth 6 %, and `-march=znver3` over generic x86-64 1.9x
(the generic build has no FMA and 128-bit vectors).

### Session 2: code generation (`results/2026-10-18-session-2-codegen`)

Node dec0335, 3 rounds, the same MICM programs (MICM default 3,582 and
366,280 ms; vector8 1,870 and 141,722 ms; vector4 2,208 and 163,413 ms).

| Port build | 128 threads, median ms (range) | vs MICM vector8 | 1 thread, median ms (range) | vs MICM vector8 |
|---|---:|---:|---:|---:|
| defaults | 596 (594–598) | 3.14x | 45,262 (45,208–45,475) | 3.13x |
| plain lane loops (`--no-lane-pragma`) | 1,078 (1,073–1,079) | 1.74x | 93,083 (92,545–93,532) | 1.52x |
| `-mtune-ctrl=^avoid_fma_chains,^avoid_fma256_chains` (byte-identical code) | 593 (589–598) | 3.15x | 46,770 (45,509–47,385) | 3.03x |
| LU panels of 16 pivots | 588 (584–590) | 3.18x | 44,494 (44,255–44,501) | 3.19x |
| LU panels of 32 pivots | 595 (592–598) | 3.14x | 44,968 (44,321–45,217) | 3.15x |
| malloc instead of aligned_alloc (`TMC_ALIGN=0`) | 616 (614–621) | 3.03x | 47,920 (47,842–48,562) | 2.96x |

The omp simd lane loops are the Derecho result: without them the port is
1.74x and 1.52x MICM's fastest, below the target; with them 3.1x. Aligned
arrays are worth 4–6 %. The panel size changes 1–2 %, inside the spread of
the identical-binary repeat. GCC's FMA-chain tuning does not touch this code.

### Session 3: grouping, schedule, and hardware threads (`results/2026-10-18-session-3-grouping`)

Node dec1360, 3 rounds, the session-1 build (`w4`), run-time options only. At a
whole node the 10,000 cells make 313 groups of 32 for 128 threads, 2.4 per
thread, and the slowest thread does three; groups of 16 (625, 4.9 per thread)
balance better. One group per thread (79 cells, MICM's default grouping) does
not pay: a 79-cell group's working set is 2.4 MB. The static schedule, which
keeps each group on the thread that first touched its memory, gains nothing.
At 256 threads (both hardware threads of every core) MICM gains 25–30 %, the
compute-bound port loses 20 %.

| Program | 128 threads, median ms (range) | vs MICM vector8 g32 | 256 threads, median ms (range) | vs MICM vector8 g32 |
|---|---:|---:|---:|---:|
| MICM default | 3,642 (3,605–3,644) | 0.54x | 2,535 (2,534–2,574) | 0.62x |
| MICM vector8, 32-cell groups | 1,976 (1,974–1,988) | 1.00x | 1,568 (1,567–1,569) | 1.00x |
| MICM vector4, 32-cell groups | 2,299 (2,258–2,355) | 0.86x | 1,728 (1,722–1,762) | 0.91x |
| port, groups of 16 | **515 (514–518)** | **3.83x** | 618 (616–638) | 2.54x |
| port, groups of 32 (the entry's default) | 592 (591–596) | 3.33x | 775 (774–786) | 2.02x |
| port, groups of 64 | 754 (754–765) | 2.62x | 1,490 (1,274–1,534) | 1.05x |
| port, groups of 79 (one per thread at 128) | 618 (614–628) | 3.19x | 2,033 (1,991–2,091) | 0.77x |
| port, groups of 128 | 1,439 (1,403–1,489) | 1.37x | 7,158 (6,632–7,798) | 0.22x |
| port, groups of 16, static schedule | 523 (521–548) | 3.78x | 665 (664–671) | 2.36x |
| port, groups of 32, static schedule | 640 (640–652) | 3.09x | 814 (813–815) | 1.93x |
| port, groups of 79, static schedule | 593 (589–598) | 3.33x | 2,806 (2,724–2,818) | 0.56x |

The benchmark's default grouping on Derecho therefore adapts to the thread
count: groups of 32 cells, halved while there would be fewer than four groups
per thread (`bench/ts1_bench.c`), which is 16 at 128 threads and 8 at 256 for
10,000 cells, and 32 at one thread, where (part B, below) 32 stays best.

At one thread (part B, node dec1360, 3 rounds): groups of 32 take 45,252 ms
(44,395–45,279), 64 take 46,141 (45,829–46,208), and 128 take 46,447
(46,253–46,943); MICM's default 366,801, vector8 g32 140,675, vector4 g32
162,498. The port with groups of 32 is 3.11x MICM's vector8 there.

### MICM's 24 configurations at 128 threads (`results/2026-10-18-micm-configurations`)

`tools/micm_configs.sh` on node dec1596: the 24 configurations of
`micm_bench` (matrix layout standard or vector 1, 2, 4, 8, 128; LU in place or
as separate factors; Mozart or Doolittle), each with one cell group per thread
and with 32-cell groups, 3 rounds, g++ 14.3.0 with the port's flags. Every one
of the 144 runs passed the gate (`configs-t128.md` has the full tables; Mozart
configurations with standard or vector1 ordering reproduce MICM's default
result exactly, the others differ by 1e-13 at most). The fastest, median ms:

| MICM configuration | Grouping | Median ms (range) | vs MICM default |
|---|---|---:|---:|
| vector8, in-place Mozart | one group per thread | **1,803 (1,659–1,929)** | 1.96x |
| vector8, in-place Mozart | 32 cells | 1,864 (1,841–1,985) | 1.89x |
| vector8, in-place Doolittle | 32 cells | 1,865 (1,723–1,961) | 1.89x |
| vector8, in-place Doolittle | one per thread | 1,954 (1,895–2,051) | 1.81x |
| vector8, separate Mozart | 32 cells | 2,080 (2,047–2,102) | 1.70x |
| vector4, in-place Mozart | one per thread | 2,198 (2,180–2,208) | 1.61x |
| standard, in-place Mozart (the default) | one per thread | 3,528 (3,188–3,608) | 1.00x |

As on the Mac, vector8 with in-place Mozart LU is MICM's fastest; on Derecho
its one-group-per-thread grouping (78 cells per thread, balanced) edges the
32-cell grouping. Against it the port's groups of 16 (515 ms on dec1360, in
another session) are 3.5x, and its groups of 32 (592 ms) 3.0x.

### Verification of the Derecho build

`tools/verify_port.sh` on the committed defaults (width 4, omp simd lanes,
aligned arrays, fast norm, adaptive grouping), gcc 14.3.0, no warnings:

- **The 680 cases** (`results/2026-10-18-verify-default/tests.txt`): 673
  PASS, 0 FAIL, 7 N/A, every id of `reference/micm_cases.txt` exactly once;
  the same seven N/A as on the Mac (C++-only cases). The entry's two extra
  checks pass: generated against generic kernels on the 100-cell input,
  3.2e-15; the gate on that input, 4.0e-14.
- **The gate on the ranking input** (seed 2026, against MICM's Derecho result):
  PASS at every setting, largest relative difference 1.27e-13 (generic
  kernels 1.4e-13): the default (groups of 16 at 128 threads, 511 ms), one
  thread (44,685 ms), groups of 16, 64, and 128, the generic kernels, and the
  static schedule.
- **A held-out input** (seed 7, never used in tuning;
  `results/2026-10-18-verify-heldout-seed7`): the suite again 673/0/7, and the
  gate PASS at the default, one thread, groups of 16 and 64, and the generic
  kernels, largest relative difference 9.3e-14; 43,630 ms at one thread and
  502 ms with groups of 16 at 128 threads.

### MICM's leading configurations at 256 threads (`results/2026-10-18-micm-configurations/timing-t256`)

`pbs/micm_t256.sh` on node dec1360, 3 rounds: both hardware threads of every
core. MICM gains 20–30 % from them (the port loses, session 3), so MICM's
fastest on a whole node is at 256 threads:

| MICM configuration | Grouping | 256 threads, median ms (range) | 128 threads (sweep) |
|---|---|---:|---:|
| vector8, in-place Mozart | one group per thread | **1,403 (1,401–1,407)** | 1,803 |
| vector8, in-place Doolittle | one per thread | 1,448 (1,443–1,457) | 1,954 |
| vector4, in-place Mozart | one per thread | 1,526 (1,525–1,535) | 2,198 |
| vector8, in-place Mozart | 32 cells | 1,574 (1,572–1,593) | 1,864 |
| standard, in-place Mozart (the default) | one per thread | 2,532 (2,527–2,598) | 3,528 |

### Session 4: the one-thread knobs (`results/2026-10-18-session-4-one-thread-knobs`)

Node dec2446, 3 rounds. A `perf` profile of one thread (node dec2433, 1,000
cells; `results/2026-10-18-profile-w4`) put the time at LU 19 %, Jacobian
assembly 17 %, forward and backward solves 26 %, forcing gather 10 %, the
library's error norm 6 %, with 1.17 instructions per cycle and 35 % of L1
data accesses missing: a block's LU matrix is 73 KB against a 32 KB L1, so the
kernels stream from L2. The knobs tried against that profile:

| Port build | 128 threads, median ms | 1 thread, median ms (range) | vs MICM vector8 g32 (141,139 ms) |
|---|---:|---:|---:|
| before the fast norm (session-1 build) | 595 | 45,086 (44,769–45,214) | 3.13x |
| **fast norm** (per-lane partial sums; now the default) | 588 | 44,139 (43,926–44,183) | 3.20x |
| fast norm, solve rows split from 8 terms (`--split-min 8`) | 590 | 43,825 (43,806–43,889) | 3.22x |
| fast norm, split from 6 terms | 593 | 43,924 (43,837–44,184) | 3.21x |
| fast norm, 8 gathers side by side (`--interleave 8`) | 588 | 44,477 (44,316–44,546) | 3.17x |
| fast norm, 2 gathers side by side | 592 | 44,476 (44,151–45,050) | 3.17x |
| fast norm, two blocks per solve iteration (`--solve-pairs`) | 602 | 44,139 (44,032–44,244) | 3.20x |
| fast norm, groups of 16 cells | | 43,291 (43,281–43,423) | 3.26x |

The error norm's serial chain of 6,688 dependent adds per internal step was
worth 2 %; groups of 16 are worth another 2 % at one thread, unlike on the M5,
where 16 cost 30 %; the split threshold, interleave, and paired solves are
within the noise and stay at the entry's values. The paired solves, which
give the core two independent chains, change nothing, consistent with the
profile: the solves wait on L2, not on the chains.

(The final session follows as it completes.)
