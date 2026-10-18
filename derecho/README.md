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

(in progress; see `results/`)
