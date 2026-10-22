# Verification of tier2/astra/xhigh

Branch `tier2/astra/xhigh`, commit 1c9112c Verify both compilers, scenario gates and benchmark results; compilers clang, gcc-15; OPTFLAGS `-O3 -DNDEBUG`; 18 threads; 2026-10-22 13:32.

## Build (fresh worktree)

| Compiler | make env | make cpu | make bench | Clean build | Warnings | Generated C |
|---|---|---|---|---:|---:|---|
| clang | exit 0 (12 s) | exit 0 (50.2 s) | exit 0 (21.1 s) | 71.2 s | 0 | 2 files, 91,513 lines, 3.17 MB |
| gcc-15 | exit 0 (10 s) | exit 0 (37.8 s) | exit 0 (35.7 s) | 73.6 s | 0 | 2 files, 91,513 lines, 3.17 MB |

clang: generated C: `generated/ts1.c`, `generated/ts1_specialized.c`

gcc-15: generated C: `generated/ts1.c`, `generated/ts1_specialized.c`

## MICM's test cases (the Tier 1 gate)

| Compiler | make test-cpu | PASS | FAIL | N/A | Ids missing | Ids repeated | Unknown ids |
|---|---|---:|---:|---:|---:|---:|---:|
| clang | exit 0 (12 s) | 680 | 0 | 0 | 0 | 0 | 0 |
| gcc-15 | exit 0 (18 s) | 680 | 0 | 0 | 0 | 0 | 0 |

## Tier 1 TS1 gate: bin/ts1_bench_cpu (10,000 cells, seed 2026, 100 steps of 30 s)

| Compiler | Exit | elapsed_ms | compare_final.py (1e-10) |
|---|---|---:|---|
| clang | 0 | 2892.489910 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.03e-13 (BCARY, cell 4179) (rtol 1e-10, atol 1e-30) |
| gcc-15 | 0 | 3210.751000 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.14e-13 (BCARY, cell 6848) (rtol 1e-10, atol 1e-30) |

## Specialized TS1 gate: bin/ts1_bench_specialized (10,000 cells, 100 steps of 30 s)

| Compiler | Input | Exit | elapsed_ms | Internal steps | compare_final.py (1e-10) |
|---|---|---|---:|---:|---|
| clang | seed 2026 | 0 | 1964.349031 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.96e-13 (BCARY, cell 4179) (rtol 1e-10, atol 1e-30) |
| clang | input 1 | 0 | 1997.138977 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.4e-13 (BCARYO2VBS, cell 9174) (rtol 1e-10, atol 1e-30) |
| gcc-15 | seed 2026 | 0 | 2098.971000 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.75e-13 (BCARY, cell 4179) (rtol 1e-10, atol 1e-30) |
| gcc-15 | input 1 | 0 | 2128.223000 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.5e-13 (BCARYO2VBS, cell 2723) (rtol 1e-10, atol 1e-30) |

## The generator's tests: make test-specialized

| Compiler | Exit | Time | PASS lines | FAIL lines | Last lines |
|---|---|---:|---:|---:|---|
| clang | 0 | 4 s | 70 | 0 | empty cells=32 strict=1 PASS / empty cells=33 strict=0 PASS / empty cells=33 strict=1 PASS |
| gcc-15 | 0 | 4 s | 70 | 0 | empty cells=32 strict=1 PASS / empty cells=33 strict=0 PASS / empty cells=33 strict=1 PASS |

## Day-night scenarios (--threads 18, --rtol 1e-6 --atol 1e-15, 32-cell groups)

Agreement is score.py against MICM's result at the same settings, with `--level 10`; accuracy is score.py against the reference solution, for scale.

| Compiler | Scenario | Exit | elapsed_ms | Internal steps (MICM's) | Agreement with MICM | Worst species | Gate | Level 10 | Accuracy (MICM's) |
|---|---|---|---:|---|---:|---|---|---|---|
| clang | arctic-summer-10min | 0 | 795.277119 | 103211 (103211) | 14.299 | ISOP 13.843 | PASS | PASS | 5.045 (5.04) |
| clang | night-start-5min | 0 | 1034.361839 | 132855 (132855) | 14.45 | ISOP 14.0 | PASS | PASS | 4.729 (4.73) |
| clang | south-48h-30min | 0 | 731.298208 | 95163 (95163) | 14.337 | ISOP 13.769 | PASS | PASS | 4.966 (4.97) |
| clang | spring45-20min | 0 | 466.338873 | 60435 (60435) | 14.342 | ISOP 13.761 | PASS | PASS | 4.842 (4.84) |
| clang | tropics-polluted-15min | 0 | 578.296185 | 75293 (75293) | 14.331 | ISOPO2VBS 13.434 | PASS | PASS | 5.1 (5.10) |
| clang | clean-10min | 0 | 652.070045 | 84273 (84273) | 14.342 | ISOP 13.793 | PASS | PASS | 4.907 (4.91) |
| clang | diurnal-10min | 0 | 710.495949 | 90980 (90980) | 14.285 | ISOP 13.656 | PASS | PASS | 5.052 (5.05) |
| clang | diurnal-30min | 0 | 375.899076 | 48619 (48619) | 14.377 | ISOP 13.828 | PASS | PASS | 4.748 (4.75) |
| clang | diurnal-5min | 0 | 1138.568640 | 145891 (145891) | 14.23 | ISOP 13.688 | PASS | PASS | 5.214 (5.21) |
| clang | polluted-10min | 0 | 836.148739 | 106393 (106393) | 14.047 | ISOPNO3 13.305 | PASS | PASS | 5.398 (5.40) |
| gcc-15 | arctic-summer-10min | 0 | 838.794000 | 103211 (103211) | 14.289 | ISOP 13.817 | PASS | PASS | 5.045 (5.04) |
| gcc-15 | night-start-5min | 0 | 1097.667000 | 132855 (132855) | 14.442 | ISOP 13.989 | PASS | PASS | 4.729 (4.73) |
| gcc-15 | south-48h-30min | 0 | 753.295000 | 95163 (95163) | 14.335 | ISOP 13.772 | PASS | PASS | 4.966 (4.97) |
| gcc-15 | spring45-20min | 0 | 483.351000 | 60435 (60435) | 14.345 | ISOP 13.763 | PASS | PASS | 4.842 (4.84) |
| gcc-15 | tropics-polluted-15min | 0 | 607.535000 | 75293 (75293) | 14.333 | ISOP 13.51 | PASS | PASS | 5.1 (5.10) |
| gcc-15 | clean-10min | 0 | 696.517000 | 84273 (84273) | 14.335 | ISOP 13.811 | PASS | PASS | 4.907 (4.91) |
| gcc-15 | diurnal-10min | 0 | 741.427000 | 90980 (90980) | 14.292 | ISOP 13.687 | PASS | PASS | 5.052 (5.05) |
| gcc-15 | diurnal-30min | 0 | 385.932000 | 48619 (48619) | 14.363 | ISOP 13.791 | PASS | PASS | 4.748 (4.75) |
| gcc-15 | diurnal-5min | 0 | 1197.383000 | 145891 (145891) | 14.232 | ISOPNO3 13.537 | PASS | PASS | 5.214 (5.21) |
| gcc-15 | polluted-10min | 0 | 877.179000 | 106393 (106393) | 14.05 | ISOP 13.432 | PASS | PASS | 5.398 (5.40) |

- clang: level 10 and the gate on every scenario: yes (10 scenarios)
- gcc-15: level 10 and the gate on every scenario: yes (10 scenarios)

## Build audit

`make -n -B cpu bench CC=clang`: exit 0; flags used: `--algebra-batch --chemistry-batch --chemistry-batch 8 --algebra-batch 8 --function-lines 256 --tile-size 4  --file --function-lines --prefix --prune; --symbol --tile-size --yes; -DCHEM_GROUP_SIZE=32 -DCHEM_SPECIALIZED -DNDEBUG -Iexamples -Iinclude -L~/homebrew/opt/libomp/lib -L~/homebrew/opt/libomp/lib -lomp -lm -O3 -O3 -DNDEBUG -std=c11 -Wall -Wextra -Wpedantic -Xpreprocessor -fopenmp -isystem ~/homebrew/opt/libomp/include -Iinclude -Wall -Wextra -Wpedantic -Xpreprocessor -c -f -fopenmp -isystem -lm -lomp -p -s -std=c11`

Flags outside the rules (review): none

- `bin/test_cpu` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_specialized` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_bench_cpu` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_bench_specialized` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_scenario_specialized` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/tutorial` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib

C sources scanned: 82 (tracked and generated, outside conda); lines matching a banned construct or starting a process: 0 (all in `logs/banned.txt`)


