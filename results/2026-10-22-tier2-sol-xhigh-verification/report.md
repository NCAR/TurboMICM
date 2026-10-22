# Verification of tier2/sol/xhigh

Branch `tier2/sol/xhigh`, commit 7af0738 Record cross-compiler gates, benchmark results and completed specialization audit; compilers clang, gcc-15; OPTFLAGS `-O3 -DNDEBUG`; 18 threads; 2026-10-22 18:43.

## Build (fresh worktree)

| Compiler | make env | make cpu | make bench | Clean build | Warnings | Generated C |
|---|---|---|---|---:|---:|---|
| clang | exit 0 (14 s) | exit 0 (1.9 s) | exit 0 (16.4 s) | 18.3 s | 0 | 2 files, 114,227 lines, 3.77 MB |
| gcc-15 | exit 0 (11 s) | exit 0 (4.8 s) | exit 0 (41.3 s) | 46.1 s | 0 | 2 files, 114,227 lines, 3.77 MB |

clang: generated C: `generated/sp_ts1.c`, `generated/ts1.c`

gcc-15: generated C: `generated/sp_ts1.c`, `generated/ts1.c`

## MICM's test cases (the Tier 1 gate)

| Compiler | make test-cpu | PASS | FAIL | N/A | Ids missing | Ids repeated | Unknown ids |
|---|---|---:|---:|---:|---:|---:|---:|
| clang | exit 0 (5 s) | 675 | 0 | 5 | 0 | 0 | 0 |
| gcc-15 | exit 0 (8 s) | 675 | 0 | 5 | 0 | 0 | 0 |

## Tier 1 TS1 gate: bin/ts1_bench_cpu (10,000 cells, seed 2026, 100 steps of 30 s)

| Compiler | Exit | elapsed_ms | compare_final.py (1e-10) |
|---|---|---:|---|
| clang | 0 | 3405.583 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 0 (rtol 1e-10, atol 1e-30) |
| gcc-15 | 0 | 3850.492 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.43e-13 (BCARYO2VBS, cell 4054) (rtol 1e-10, atol 1e-30) |

## Specialized TS1 gate: bin/ts1_bench_specialized (10,000 cells, 100 steps of 30 s)

| Compiler | Input | Exit | elapsed_ms | Internal steps | compare_final.py (1e-10) |
|---|---|---|---:|---:|---|
| clang | seed 2026 | 0 | 2667.143 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.31e-14 (BCARYO2VBS, cell 763) (rtol 1e-10, atol 1e-30) |
| clang | input 1 | 0 | 2728.576 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.95e-14 (BCARY, cell 4191) (rtol 1e-10, atol 1e-30) |
| gcc-15 | seed 2026 | 0 | 2916.941 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.43e-13 (BCARYO2VBS, cell 4054) (rtol 1e-10, atol 1e-30) |
| gcc-15 | input 1 | 0 | 2969.333 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.52e-13 (BCARYO2VBS, cell 2723) (rtol 1e-10, atol 1e-30) |

## The generator's tests: make test-specialized

| Compiler | Exit | Time | PASS lines | FAIL lines | Last lines |
|---|---|---:|---:|---:|---|
| clang | 0 | 2 s | 70 | 0 | Ran 4 tests in 0.490s /  / OK |
| gcc-15 | 0 | 4 s | 70 | 0 | Ran 4 tests in 0.490s /  / OK |

## Day-night scenarios (--threads 18, --rtol 1e-6 --atol 1e-15, 32-cell groups)

Agreement is score.py against MICM's result at the same settings, with `--level 10`; accuracy is score.py against the reference solution, for scale.

| Compiler | Scenario | Exit | elapsed_ms | Internal steps (MICM's) | Agreement with MICM | Worst species | Gate | Level 10 | Accuracy (MICM's) |
|---|---|---|---:|---|---:|---|---|---|---|
| clang | arctic-summer-10min | 0 | 1095.288 | 103211 (103211) | 14.295 | ISOP 13.828 | PASS | PASS | 5.045 (5.04) |
| clang | night-start-5min | 0 | 1421.875 | 132855 (132855) | 14.455 | ISOP 14.016 | PASS | PASS | 4.729 (4.73) |
| clang | south-48h-30min | 0 | 1004.194 | 95163 (95163) | 14.335 | ISOP 13.755 | PASS | PASS | 4.966 (4.97) |
| clang | spring45-20min | 0 | 641.116 | 60435 (60435) | 14.343 | ISOP 13.809 | PASS | PASS | 4.842 (4.84) |
| clang | tropics-polluted-15min | 0 | 820.690 | 75293 (75293) | 14.338 | ISOP 13.633 | PASS | PASS | 5.1 (5.10) |
| clang | clean-10min | 0 | 900.977 | 84273 (84273) | 14.35 | ISOP 13.842 | PASS | PASS | 4.907 (4.91) |
| clang | diurnal-10min | 0 | 977.144 | 90980 (90980) | 14.287 | ISOPNO3 13.649 | PASS | PASS | 5.052 (5.05) |
| clang | diurnal-30min | 0 | 513.293 | 48619 (48619) | 14.369 | ISOP 13.81 | PASS | PASS | 4.748 (4.75) |
| clang | diurnal-5min | 0 | 1598.145 | 145891 (145891) | 14.235 | ISOP 13.631 | PASS | PASS | 5.214 (5.21) |
| clang | polluted-10min | 0 | 1157.600 | 106393 (106393) | 14.076 | ISOP 13.43 | PASS | PASS | 5.398 (5.40) |
| gcc-15 | arctic-summer-10min | 0 | 1191.586 | 103211 (103211) | 14.304 | ISOP 13.831 | PASS | PASS | 5.045 (5.04) |
| gcc-15 | night-start-5min | 0 | 1564.793 | 132855 (132855) | 14.456 | ISOP 13.998 | PASS | PASS | 4.729 (4.73) |
| gcc-15 | south-48h-30min | 0 | 1109.165 | 95163 (95163) | 14.342 | ISOPO2VBS 13.726 | PASS | PASS | 4.966 (4.97) |
| gcc-15 | spring45-20min | 0 | 692.937 | 60435 (60435) | 14.347 | ISOP 13.781 | PASS | PASS | 4.842 (4.84) |
| gcc-15 | tropics-polluted-15min | 0 | 863.907 | 75293 (75293) | 14.337 | MTERP 13.626 | PASS | PASS | 5.1 (5.10) |
| gcc-15 | clean-10min | 0 | 972.542 | 84273 (84273) | 14.345 | ISOP 13.845 | PASS | PASS | 4.907 (4.91) |
| gcc-15 | diurnal-10min | 0 | 1060.123 | 90980 (90980) | 14.296 | ISOP 13.688 | PASS | PASS | 5.052 (5.05) |
| gcc-15 | diurnal-30min | 0 | 560.672 | 48619 (48619) | 14.366 | ISOPNO3 13.591 | PASS | PASS | 4.748 (4.75) |
| gcc-15 | diurnal-5min | 0 | 1700.356 | 145891 (145891) | 14.243 | ISOP 13.701 | PASS | PASS | 5.214 (5.21) |
| gcc-15 | polluted-10min | 0 | 1254.570 | 106393 (106393) | 14.083 | ISOPNO3 13.424 | PASS | PASS | 5.398 (5.40) |

- clang: level 10 and the gate on every scenario: yes (10 scenarios)
- gcc-15: level 10 and the gate on every scenario: yes (10 scenarios)

## Build audit

`make -n -B cpu bench CC=clang`: exit 0; flags used: `--algebra-batch --batch --file --fuse --joint-batch --lu-order --prefix --split --symbol --tile --width -DCK_DEFAULT_GROUP_SIZE=32 -DCK_SPECIALIZED -DNDEBUG -Igenerated -Iinclude -Iinclude -Igenerated -Isrc -Isrc -L~/homebrew/opt/libomp/lib -L~/homebrew/opt/libomp/lib -lomp -O3 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Xpreprocessor -Xpreprocessor -fopenmp -isystem ~/homebrew/opt/libomp/include -c -fopenmp -isystem -lm -lomp -p -s -std=c11 -x`

Flags outside the rules (review): none

- `bin/test_cpu` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_specialized` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_bench_cpu` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_bench_specialized` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_scenario_specialized` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib

C sources scanned: 65 (tracked and generated, outside conda); lines matching a banned construct or starting a process: 0 (all in `logs/banned.txt`)


