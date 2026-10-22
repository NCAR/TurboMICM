# Verification of tier2/astra/max

Branch `tier2/astra/max`, commit cf16213 Tune generated kernels and document complete compiler and scenario validation; compilers clang, gcc-15; OPTFLAGS `-O3 -DNDEBUG`; 18 threads; 2026-10-22 11:18.

## Build (fresh worktree)

| Compiler | make env | make cpu | make bench | Clean build | Warnings | Generated C |
|---|---|---|---|---:|---:|---|
| clang | exit 0 (11 s) | exit 0 (4.0 s) | exit 0 (13.8 s) | 17.8 s | 0 | 4 files, 115,752 lines, 3.05 MB |
| gcc-15 | exit 0 (11 s) | exit 0 (6.0 s) | exit 0 (46.7 s) | 52.6 s | 0 | 4 files, 115,752 lines, 3.05 MB |

clang: generated C: `build/generated/ts1.c`, `build/generated/ts1_specialized.c`, `build/generated/ts1_specialized.h`, `build/generated/ts1_specialized_kernels.c`

gcc-15: generated C: `build/generated/ts1.c`, `build/generated/ts1_specialized.c`, `build/generated/ts1_specialized.h`, `build/generated/ts1_specialized_kernels.c`

## MICM's test cases (the Tier 1 gate)

| Compiler | make test-cpu | PASS | FAIL | N/A | Ids missing | Ids repeated | Unknown ids |
|---|---|---:|---:|---:|---:|---:|---:|
| clang | exit 0 (9 s) | 680 | 0 | 0 | 0 | 0 | 0 |
| gcc-15 | exit 0 (8 s) | 680 | 0 | 0 | 0 | 0 | 0 |

## Tier 1 TS1 gate: bin/ts1_bench_cpu (10,000 cells, seed 2026, 100 steps of 30 s)

| Compiler | Exit | elapsed_ms | compare_final.py (1e-10) |
|---|---|---:|---|
| clang | 0 | 3516.993 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.96e-13 (BCARY, cell 4179) (rtol 1e-10, atol 1e-30) |
| gcc-15 | 0 | 3916.596 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.75e-13 (BCARY, cell 4179) (rtol 1e-10, atol 1e-30) |

## Specialized TS1 gate: bin/ts1_bench_specialized (10,000 cells, 100 steps of 30 s)

| Compiler | Input | Exit | elapsed_ms | Internal steps | compare_final.py (1e-10) |
|---|---|---|---:|---:|---|
| clang | seed 2026 | 0 | 2947.367 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.96e-13 (BCARY, cell 4179) (rtol 1e-10, atol 1e-30) |
| clang | input 1 | 0 | 3045.176 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.4e-13 (BCARYO2VBS, cell 9174) (rtol 1e-10, atol 1e-30) |
| gcc-15 | seed 2026 | 0 | 3240.309 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.75e-13 (BCARY, cell 4179) (rtol 1e-10, atol 1e-30) |
| gcc-15 | input 1 | 0 | 3314.159 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.5e-13 (BCARYO2VBS, cell 2723) (rtol 1e-10, atol 1e-30) |

## The generator's tests: make test-specialized

| Compiler | Exit | Time | PASS lines | FAIL lines | Last lines |
|---|---|---:|---:|---:|---|
| clang | 0 | 4 s | 6 | 0 | PASS robertson_specialized.h: rates, forcing, Jacobian, LU, solve and adaptive trajectories (1e-10) / PASS all_rates_specialized.h: rates, forcing, Jacobian, LU, solve and adaptive trajectories (1e-10 |
| gcc-15 | 0 | 5 s | 6 | 0 | PASS robertson_specialized.h: rates, forcing, Jacobian, LU, solve and adaptive trajectories (1e-10) / PASS all_rates_specialized.h: rates, forcing, Jacobian, LU, solve and adaptive trajectories (1e-10 |

## Day-night scenarios (--threads 18, --rtol 1e-6 --atol 1e-15, 32-cell groups)

Agreement is score.py against MICM's result at the same settings, with `--level 10`; accuracy is score.py against the reference solution, for scale.

| Compiler | Scenario | Exit | elapsed_ms | Internal steps (MICM's) | Agreement with MICM | Worst species | Gate | Level 10 | Accuracy (MICM's) |
|---|---|---|---:|---|---:|---|---|---|---|
| clang | arctic-summer-10min | 0 | 1176.420 | 103211 (103211) | 14.298 | MVK 13.85 | PASS | PASS | 5.045 (5.04) |
| clang | night-start-5min | 0 | 1547.699 | 132855 (132855) | 14.448 | ISOP 13.994 | PASS | PASS | 4.729 (4.73) |
| clang | south-48h-30min | 0 | 1091.620 | 95163 (95163) | 14.343 | ISOPAO2 13.629 | PASS | PASS | 4.966 (4.97) |
| clang | spring45-20min | 0 | 699.919 | 60435 (60435) | 14.344 | ISOP 13.769 | PASS | PASS | 4.842 (4.84) |
| clang | tropics-polluted-15min | 0 | 875.116 | 75293 (75293) | 14.341 | ISOP 13.618 | PASS | PASS | 5.1 (5.10) |
| clang | clean-10min | 0 | 972.155 | 84273 (84273) | 14.334 | ISOP 13.813 | PASS | PASS | 4.907 (4.91) |
| clang | diurnal-10min | 0 | 1053.892 | 90980 (90980) | 14.288 | ISOP 13.717 | PASS | PASS | 5.052 (5.05) |
| clang | diurnal-30min | 0 | 564.240 | 48619 (48619) | 14.365 | ISOP 13.78 | PASS | PASS | 4.748 (4.75) |
| clang | diurnal-5min | 0 | 1756.387 | 145891 (145891) | 14.229 | ISOP 13.643 | PASS | PASS | 5.214 (5.21) |
| clang | polluted-10min | 0 | 1261.093 | 106393 (106393) | 14.071 | ISOPNO3 13.428 | PASS | PASS | 5.398 (5.40) |
| gcc-15 | arctic-summer-10min | 0 | 1294.205 | 103211 (103211) | 14.289 | ISOP 13.817 | PASS | PASS | 5.045 (5.04) |
| gcc-15 | night-start-5min | 0 | 1708.588 | 132855 (132855) | 14.442 | ISOP 13.989 | PASS | PASS | 4.729 (4.73) |
| gcc-15 | south-48h-30min | 0 | 1198.653 | 95163 (95163) | 14.335 | ISOP 13.772 | PASS | PASS | 4.966 (4.97) |
| gcc-15 | spring45-20min | 0 | 754.884 | 60435 (60435) | 14.345 | ISOP 13.763 | PASS | PASS | 4.842 (4.84) |
| gcc-15 | tropics-polluted-15min | 0 | 941.099 | 75293 (75293) | 14.333 | ISOP 13.51 | PASS | PASS | 5.1 (5.10) |
| gcc-15 | clean-10min | 0 | 1064.195 | 84273 (84273) | 14.335 | ISOP 13.811 | PASS | PASS | 4.907 (4.91) |
| gcc-15 | diurnal-10min | 0 | 1158.783 | 90980 (90980) | 14.292 | ISOP 13.687 | PASS | PASS | 5.052 (5.05) |
| gcc-15 | diurnal-30min | 0 | 616.489 | 48619 (48619) | 14.363 | ISOP 13.791 | PASS | PASS | 4.748 (4.75) |
| gcc-15 | diurnal-5min | 0 | 1875.245 | 145891 (145891) | 14.232 | ISOPNO3 13.537 | PASS | PASS | 5.214 (5.21) |
| gcc-15 | polluted-10min | 0 | 1378.158 | 106393 (106393) | 14.05 | ISOP 13.432 | PASS | PASS | 5.398 (5.40) |

- clang: level 10 and the gate on every scenario: yes (10 scenarios)
- gcc-15: level 10 and the gate on every scenario: yes (10 scenarios)

## Build audit

`make -n -B cpu bench CC=clang`: exit 0; flags used: `--block --block 32 --fixed 1 --lu mozart --fusion 4 --file --fixed --fusion --lu --name --prefix --yes -DCH_DEFAULT_GROUP_SIZE=32 -DCH_KERNEL_WIDTH=32 -DCH_SPECIALIZED -DNDEBUG -Ibuild/generated -Iinclude -Isrc -Itools -L~/homebrew/opt/libomp/lib -O3 -Wall -Wextra -Wpedantic -Xpreprocessor -c -dynamiclib -fPIC -fopenmp -isystem -lm -lomp -p -s -std=c11`

Flags outside the rules (review): none

- `bin/test_native_cpu` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_specialized_all_rates` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_specialized_empty` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_specialized_no_user` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_specialized_robertson` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_specialized_ts1` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_specialized_user` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_bench_cpu` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_bench_specialized` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_scenario_specialized` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/tutorial_cpu` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib

C sources scanned: 44 (tracked and generated, outside conda); lines matching a banned construct or starting a process: 0 (all in `logs/banned.txt`)


