# Verification of tier2/sol/max

Branch `tier2/sol/max`, commit ca0119e Record compiler gates, full-size timings, build costs, and day-night scores; compilers clang, gcc-15; OPTFLAGS `-O3 -DNDEBUG`; 18 threads; 2026-10-22 17:32.

## Build (fresh worktree)

| Compiler | make env | make cpu | make bench | Clean build | Warnings | Generated C |
|---|---|---|---|---:|---:|---|
| clang | exit 0 (12 s) | exit 0 (9.4 s) | exit 0 (21.9 s) | 31.3 s | 0 | 6 files, 88,473 lines, 2.86 MB |
| gcc-15 | exit 0 (9 s) | exit 0 (15.1 s) | exit 0 (41.0 s) | 56.1 s | 0 | 6 files, 88,473 lines, 2.86 MB |

clang: generated C: `generated/example_no_user.c`, `generated/example_robertson.c`, `generated/example_user.c`, `generated/ts1.c`, `generated/ts1_specialized.c`, `generated/validation.c`

gcc-15: generated C: `generated/example_no_user.c`, `generated/example_robertson.c`, `generated/example_user.c`, `generated/ts1.c`, `generated/ts1_specialized.c`, `generated/validation.c`

## MICM's test cases (the Tier 1 gate)

| Compiler | make test-cpu | PASS | FAIL | N/A | Ids missing | Ids repeated | Unknown ids |
|---|---|---:|---:|---:|---:|---:|---:|
| clang | exit 0 (5 s) | 675 | 0 | 5 | 0 | 0 | 0 |
| gcc-15 | exit 0 (5 s) | 675 | 0 | 5 | 0 | 0 | 0 |

## Tier 1 TS1 gate: bin/ts1_bench_cpu (10,000 cells, seed 2026, 100 steps of 30 s)

| Compiler | Exit | elapsed_ms | compare_final.py (1e-10) |
|---|---|---:|---|
| clang | 0 | 2871.420 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.26e-13 (BCARY, cell 4179) (rtol 1e-10, atol 1e-30) |
| gcc-15 | 0 | 3087.044 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.16e-13 (BCARYO2VBS, cell 6455) (rtol 1e-10, atol 1e-30) |

## Specialized TS1 gate: bin/ts1_bench_specialized (10,000 cells, 100 steps of 30 s)

| Compiler | Input | Exit | elapsed_ms | Internal steps | compare_final.py (1e-10) |
|---|---|---|---:|---:|---|
| clang | seed 2026 | 0 | 2347.715 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.26e-13 (BCARY, cell 4179) (rtol 1e-10, atol 1e-30) |
| clang | input 1 | 0 | 2373.939 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 9.98e-14 (BCARYO2VBS, cell 5192) (rtol 1e-10, atol 1e-30) |
| gcc-15 | seed 2026 | 0 | 2247.223 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.16e-13 (BCARYO2VBS, cell 6455) (rtol 1e-10, atol 1e-30) |
| gcc-15 | input 1 | 0 | 2304.990 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.48e-13 (BCARYO2VBS, cell 9174) (rtol 1e-10, atol 1e-30) |

## The generator's tests: make test-specialized

| Compiler | Exit | Time | PASS lines | FAIL lines | Last lines |
|---|---|---:|---:|---:|---|
| clang | 0 | 32 s | 42 | 0 | ts1_variant cells=7 PASS rates, chemistry, LU, solves, adaptive integration / ts1_variant cells=32 PASS rates, chemistry, LU, solves, adaptive integration / ts1_variant cells=33 PASS rates, chemistry, |
| gcc-15 | 0 | 56 s | 42 | 0 | ts1_variant cells=7 PASS rates, chemistry, LU, solves, adaptive integration / ts1_variant cells=32 PASS rates, chemistry, LU, solves, adaptive integration / ts1_variant cells=33 PASS rates, chemistry, |

## Day-night scenarios (--threads 18, --rtol 1e-6 --atol 1e-15, 32-cell groups)

Agreement is score.py against MICM's result at the same settings, with `--level 10`; accuracy is score.py against the reference solution, for scale.

| Compiler | Scenario | Exit | elapsed_ms | Internal steps (MICM's) | Agreement with MICM | Worst species | Gate | Level 10 | Accuracy (MICM's) |
|---|---|---|---:|---|---:|---|---|---|---|
| clang | arctic-summer-10min | 0 | 944.710 | 103211 (103211) | 14.297 | ISOP 13.834 | PASS | PASS | 5.045 (5.04) |
| clang | night-start-5min | 0 | 1212.399 | 132855 (132855) | 14.468 | ISOP 14.019 | PASS | PASS | 4.729 (4.73) |
| clang | south-48h-30min | 0 | 874.029 | 95163 (95163) | 14.338 | ISOP 13.848 | PASS | PASS | 4.966 (4.97) |
| clang | spring45-20min | 0 | 546.634 | 60435 (60435) | 14.35 | ISOP 13.752 | PASS | PASS | 4.842 (4.84) |
| clang | tropics-polluted-15min | 0 | 691.822 | 75293 (75293) | 14.339 | ISOP 13.619 | PASS | PASS | 5.1 (5.10) |
| clang | clean-10min | 0 | 765.426 | 84273 (84273) | 14.343 | ISOP 13.809 | PASS | PASS | 4.907 (4.91) |
| clang | diurnal-10min | 0 | 829.755 | 90980 (90980) | 14.285 | ISOP 13.71 | PASS | PASS | 5.052 (5.05) |
| clang | diurnal-30min | 0 | 443.623 | 48619 (48619) | 14.379 | ISOPNO3 13.775 | PASS | PASS | 4.748 (4.75) |
| clang | diurnal-5min | 0 | 1365.141 | 145891 (145891) | 14.242 | ISOPNO3 13.574 | PASS | PASS | 5.214 (5.21) |
| clang | polluted-10min | 0 | 968.468 | 106393 (106393) | 14.065 | ISOP 13.435 | PASS | PASS | 5.398 (5.40) |
| gcc-15 | arctic-summer-10min | 0 | 894.722 | 103211 (103211) | 14.297 | ISOP 13.827 | PASS | PASS | 5.045 (5.04) |
| gcc-15 | night-start-5min | 0 | 1169.134 | 132855 (132855) | 14.456 | ISOP 13.993 | PASS | PASS | 4.729 (4.73) |
| gcc-15 | south-48h-30min | 0 | 835.677 | 95163 (95163) | 14.338 | ISOP 13.789 | PASS | PASS | 4.966 (4.97) |
| gcc-15 | spring45-20min | 0 | 521.861 | 60435 (60435) | 14.352 | ISOP 13.816 | PASS | PASS | 4.842 (4.84) |
| gcc-15 | tropics-polluted-15min | 0 | 668.420 | 75293 (75293) | 14.332 | MTERP 13.665 | PASS | PASS | 5.1 (5.10) |
| gcc-15 | clean-10min | 0 | 735.838 | 84273 (84273) | 14.359 | ISOP 13.851 | PASS | PASS | 4.907 (4.91) |
| gcc-15 | diurnal-10min | 0 | 787.613 | 90980 (90980) | 14.286 | ISOP 13.649 | PASS | PASS | 5.052 (5.05) |
| gcc-15 | diurnal-30min | 0 | 422.807 | 48619 (48619) | 14.368 | ISOP 13.796 | PASS | PASS | 4.748 (4.75) |
| gcc-15 | diurnal-5min | 0 | 1288.146 | 145891 (145891) | 14.236 | ISOP 13.689 | PASS | PASS | 5.214 (5.21) |
| gcc-15 | polluted-10min | 0 | 933.691 | 106393 (106393) | 14.059 | ISOP 13.46 | PASS | PASS | 5.398 (5.40) |

- clang: level 10 and the gate on every scenario: yes (10 scenarios)
- gcc-15: level 10 and the gate on every scenario: yes (10 scenarios)

## Build audit

`make -n -B cpu bench CC=clang`: exit 0; flags used: `--cell-stride --file --fuse-kernels --fuse-stages --initialize-kernels --kernel-chunk --kernel-chunk 32 --kernel-chunk 32 --linear-chunk 1 --fuse-kernels 1 --fuse-stages 1 --specialize-rates 1 --linear unrolled --linear-loop rows --vector-error 1 --cell-stride 4 --initialize-kernels 1 --solver-tile 4 --linear --linear-chunk --linear-loop --name --prefix --solver-tile --specialize-rates --vector-error --yes -DCHEM_GROUP_SIZE=32 -DCHEM_SPECIALIZED -DNDEBUG -Igenerated -Iinclude -L~/homebrew/opt/libomp/lib -O3 -Wall -Wextra -Wpedantic -Xpreprocessor -c -dynamiclib -f -fPIC -fopenmp -isystem -lm -lomp -p -s -std=c11`

Flags outside the rules (review): none

- `bin/libchem.dylib` links: build/clang/libchem.dylib, ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/libspecialized.dylib` links: build/clang/libspecialized.dylib, ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_bench_cpu` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_bench_specialized` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_scenario_specialized` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/tutorials` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib

C sources scanned: 37 (tracked and generated, outside conda); lines matching a banned construct or starting a process: 0 (all in `logs/banned.txt`)


