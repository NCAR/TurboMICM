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

(Sessions 3 and 4, the MICM sweep, the verification, and the final session
follow below as they complete.)
