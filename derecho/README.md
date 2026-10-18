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
  threads are not scheduled separately), so `ncpus=128` takes a whole node;
  jobs started within a minute throughout. Output streams only with
  `-k oed`.
- **Conda.** `make env` uses conda from PATH or Derecho's system miniforge;
  the conda module is not loaded.

## Results

(in progress; see `results/`)
