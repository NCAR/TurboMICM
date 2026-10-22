# Verification of tier2/astra/high

Branch `tier2/astra/high`, commit f478080 Document specialization experiments and verify both compiler gates; compilers clang, gcc-15; OPTFLAGS `-O3 -DNDEBUG`; 18 threads; 2026-10-22 15:14.

## Build (fresh worktree)

| Compiler | make env | make cpu | make bench | Clean build | Warnings | Generated C |
|---|---|---|---|---:|---:|---|
| clang | exit 0 (12 s) | exit 0 (4.1 s) | exit 0 (20.9 s) | 25.0 s | 0 | 2 files, 51,433 lines, 2.43 MB |
| gcc-15 | exit 0 (10 s) | exit 0 (9.8 s) | exit 0 (27.2 s) | 37.0 s | 0 | 2 files, 51,433 lines, 2.43 MB |

clang: generated C: `build/specialized/ts1.c`, `build/ts1.c`

gcc-15: generated C: `build/specialized/ts1.c`, `build/ts1.c`

## MICM's test cases (the Tier 1 gate)

| Compiler | make test-cpu | PASS | FAIL | N/A | Ids missing | Ids repeated | Unknown ids |
|---|---|---:|---:|---:|---:|---:|---:|
| clang | exit 0 (11 s) | 680 | 0 | 0 | 0 | 0 | 0 |
| gcc-15 | exit 0 (16 s) | 680 | 0 | 0 | 0 | 0 | 0 |

## Tier 1 TS1 gate: bin/ts1_bench_cpu (10,000 cells, seed 2026, 100 steps of 30 s)

| Compiler | Exit | elapsed_ms | compare_final.py (1e-10) |
|---|---|---:|---|
| clang | 0 | 3909.883022 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 0 (rtol 1e-10, atol 1e-30) |
| gcc-15 | 0 | 4514.078000 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.43e-13 (BCARYO2VBS, cell 4054) (rtol 1e-10, atol 1e-30) |

## Specialized TS1 gate: bin/ts1_bench_specialized (10,000 cells, 100 steps of 30 s)

| Compiler | Input | Exit | elapsed_ms | Internal steps | compare_final.py (1e-10) |
|---|---|---|---:|---:|---|
| clang | seed 2026 | 0 | 3160.110950 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.31e-14 (BCARYO2VBS, cell 763) (rtol 1e-10, atol 1e-30) |
| clang | input 1 | 0 | 3230.654001 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.95e-14 (BCARY, cell 4191) (rtol 1e-10, atol 1e-30) |
| gcc-15 | seed 2026 | 0 | 3107.996000 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.43e-13 (BCARYO2VBS, cell 4054) (rtol 1e-10, atol 1e-30) |
| gcc-15 | input 1 | 0 | 3158.420000 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.52e-13 (BCARYO2VBS, cell 2723) (rtol 1e-10, atol 1e-30) |

## The generator's tests: make test-specialized

| Compiler | Exit | Time | PASS lines | FAIL lines | Last lines |
|---|---|---:|---:|---:|---|
| clang | 0 | 3 s | 6 | 0 | robertson specialized equivalence: PASS / extra specialized equivalence: PASS / mutated specialized equivalence: PASS |
| gcc-15 | 0 | 3 s | 6 | 0 | robertson specialized equivalence: PASS / extra specialized equivalence: PASS / mutated specialized equivalence: PASS |

## Day-night scenarios (--threads 18, --rtol 1e-6 --atol 1e-15, 32-cell groups)

Agreement is score.py against MICM's result at the same settings, with `--level 10`; accuracy is score.py against the reference solution, for scale.

| Compiler | Scenario | Exit | elapsed_ms | Internal steps (MICM's) | Agreement with MICM | Worst species | Gate | Level 10 | Accuracy (MICM's) |
|---|---|---|---:|---|---:|---|---|---|---|
| clang | arctic-summer-10min | 0 | 1396.778345 | 103211 (103211) | 14.295 | ISOP 13.828 | PASS | PASS | 5.045 (5.04) |
| clang | night-start-5min | 0 | 1807.821512 | 132855 (132855) | 14.455 | ISOP 14.016 | PASS | PASS | 4.729 (4.73) |
| clang | south-48h-30min | 0 | 1247.620344 | 95163 (95163) | 14.335 | ISOP 13.755 | PASS | PASS | 4.966 (4.97) |
| clang | spring45-20min | 0 | 816.117287 | 60435 (60435) | 14.343 | ISOP 13.809 | PASS | PASS | 4.842 (4.84) |
| clang | tropics-polluted-15min | 0 | 998.260021 | 75293 (75293) | 14.338 | ISOP 13.633 | PASS | PASS | 5.1 (5.10) |
| clang | clean-10min | 0 | 1132.906199 | 84273 (84273) | 14.35 | ISOP 13.842 | PASS | PASS | 4.907 (4.91) |
| clang | diurnal-10min | 0 | 1211.935997 | 90980 (90980) | 14.287 | ISOPNO3 13.649 | PASS | PASS | 5.052 (5.05) |
| clang | diurnal-30min | 0 | 643.823385 | 48619 (48619) | 14.369 | ISOP 13.81 | PASS | PASS | 4.748 (4.75) |
| clang | diurnal-5min | 0 | 1998.919487 | 145891 (145891) | 14.235 | ISOP 13.631 | PASS | PASS | 5.214 (5.21) |
| clang | polluted-10min | 0 | 1430.610895 | 106393 (106393) | 14.076 | ISOP 13.43 | PASS | PASS | 5.398 (5.40) |
| gcc-15 | arctic-summer-10min | 0 | 1317.308000 | 103211 (103211) | 14.304 | ISOP 13.831 | PASS | PASS | 5.045 (5.04) |
| gcc-15 | night-start-5min | 0 | 1732.008000 | 132855 (132855) | 14.456 | ISOP 13.998 | PASS | PASS | 4.729 (4.73) |
| gcc-15 | south-48h-30min | 0 | 1184.790000 | 95163 (95163) | 14.342 | ISOPO2VBS 13.726 | PASS | PASS | 4.966 (4.97) |
| gcc-15 | spring45-20min | 0 | 755.738000 | 60435 (60435) | 14.347 | ISOP 13.781 | PASS | PASS | 4.842 (4.84) |
| gcc-15 | tropics-polluted-15min | 0 | 950.970000 | 75293 (75293) | 14.337 | MTERP 13.626 | PASS | PASS | 5.1 (5.10) |
| gcc-15 | clean-10min | 0 | 1064.589000 | 84273 (84273) | 14.344 | ISOP 13.844 | PASS | PASS | 4.907 (4.91) |
| gcc-15 | diurnal-10min | 0 | 1155.179000 | 90980 (90980) | 14.296 | ISOP 13.688 | PASS | PASS | 5.052 (5.05) |
| gcc-15 | diurnal-30min | 0 | 602.426000 | 48619 (48619) | 14.366 | ISOPNO3 13.591 | PASS | PASS | 4.748 (4.75) |
| gcc-15 | diurnal-5min | 0 | 1892.329000 | 145891 (145891) | 14.243 | ISOP 13.701 | PASS | PASS | 5.214 (5.21) |
| gcc-15 | polluted-10min | 0 | 1356.426000 | 106393 (106393) | 14.083 | ISOPNO3 13.424 | PASS | PASS | 5.398 (5.40) |

- clang: level 10 and the gate on every scenario: yes (10 scenarios)
- gcc-15: level 10 and the gate on every scenario: yes (10 scenarios)

## Build audit

`make -n -B cpu bench CC=clang`: exit 0; flags used: `--chunk --depfile --file --fusion --prefix --symbol --width --yes -DCHEM_GROUP_SIZE=32 -DCHEM_SPECIALIZED -DNDEBUG -Iinclude -Isrc -L~/homebrew/opt/libomp/lib -L~/homebrew/opt/libomp/lib -lomp -lm -O3 -Wall -Wextra -Wpedantic -Xpreprocessor -c -fopenmp -isystem -lm -lomp -p -s -std=c11 -std=c11 -Wall -Wextra -Wpedantic -O3 -DNDEBUG -Xpreprocessor -fopenmp -isystem ~/homebrew/opt/libomp/include -Iinclude -Isrc -x`

Flags outside the rules (review): none

- `bin/test_config` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_config_extra` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_cpu` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_specialized` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_bench_cpu` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_bench_specialized` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_scenario_specialized` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/tutorial_cpu` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib

C sources scanned: 75 (tracked and generated, outside conda); lines matching a banned construct or starting a process: 0 (all in `logs/banned.txt`)


