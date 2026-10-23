# Verification of tier2/luna/max

Branch `tier2/luna/max`, commit 3562d5a Generate specialized TS1 solver kernels; compilers clang, gcc-15; OPTFLAGS `-O3 -DNDEBUG`; 18 threads; 2026-10-22 21:51.

## Build (fresh worktree)

| Compiler | make env | make cpu | make bench | Clean build | Warnings | Generated C |
|---|---|---|---|---:|---:|---|
| clang | exit 0 (12 s) | exit 0 (3.3 s) | exit 0 (2.2 s) | 5.5 s | 0 | 2 files, 29,385 lines, 1.67 MB |
| gcc-15 | exit 0 (11 s) | exit 0 (4.2 s) | exit 0 (9.0 s) | 13.2 s | 0 | 2 files, 29,385 lines, 1.67 MB |

clang: generated C: `build/ts1_model.c`, `build/ts1_model_specialized.c`

gcc-15: generated C: `build/ts1_model.c`, `build/ts1_model_specialized.c`

## MICM's test cases (the Tier 1 gate)

| Compiler | make test-cpu | PASS | FAIL | N/A | Ids missing | Ids repeated | Unknown ids |
|---|---|---:|---:|---:|---:|---:|---:|
| clang | exit 0 (5 s) | 671 | 0 | 9 | 0 | 0 | 0 |
| gcc-15 | exit 0 (7 s) | 671 | 0 | 9 | 0 | 0 | 0 |

## Tier 1 TS1 gate: bin/ts1_bench_cpu (10,000 cells, seed 2026, 100 steps of 30 s)

| Compiler | Exit | elapsed_ms | compare_final.py (1e-10) |
|---|---|---:|---|
| clang | 0 | 14274.321 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.78e-13 (BCARY, cell 4179) (rtol 1e-10, atol 1e-30) |
| gcc-15 | 0 | 15509.089 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.2e-13 (BCARY, cell 4179) (rtol 1e-10, atol 1e-30) |

## Specialized TS1 gate: bin/ts1_bench_specialized (10,000 cells, 100 steps of 30 s)

| Compiler | Input | Exit | elapsed_ms | Internal steps | compare_final.py (1e-10) |
|---|---|---|---:|---:|---|
| clang | seed 2026 | 0 | 5393.324 | 16200 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.78e-13 (BCARY, cell 4179) (rtol 1e-10, atol 1e-30) |
| clang | input 1 | 0 | 5447.704 | 16200 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.74e-13 (BCARYO2VBS, cell 2723) (rtol 1e-10, atol 1e-30) |
| gcc-15 | seed 2026 | 0 | 6078.039 | 16200 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.2e-13 (BCARY, cell 4179) (rtol 1e-10, atol 1e-30) |
| gcc-15 | input 1 | 0 | 6104.566 | 16200 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.78e-13 (BCARYO2VBS, cell 2723) (rtol 1e-10, atol 1e-30) |

## The generator's tests: make test-specialized

| Compiler | Exit | Time | PASS lines | FAIL lines | Last lines |
|---|---|---:|---:|---:|---|
| clang | 0 | 3 s | 4 | 0 | PASS specialized solver agrees with generic solver: reaction rates no user defined / bin/test_specialized_tutorial_no_user / PASS specialized solver agrees with generic solver: reaction rates no user  |
| gcc-15 | 0 | 10 s | 4 | 0 | PASS specialized solver agrees with generic solver: reaction rates no user defined / bin/test_specialized_tutorial_no_user / PASS specialized solver agrees with generic solver: reaction rates no user  |

## Day-night scenarios (--threads 18, --rtol 1e-6 --atol 1e-15, 32-cell groups)

Agreement is score.py against MICM's result at the same settings, with `--level 10`; accuracy is score.py against the reference solution, for scale.

| Compiler | Scenario | Exit | elapsed_ms | Internal steps (MICM's) | Agreement with MICM | Worst species | Gate | Level 10 | Accuracy (MICM's) |
|---|---|---|---:|---|---:|---|---|---|---|
| clang | arctic-summer-10min | 0 | 2200.309 | 103211 (103211) | 14.288 | ISOP 13.807 | PASS | PASS | 5.045 (5.04) |
| clang | night-start-5min | 0 | 2878.731 | 132855 (132855) | 14.454 | ISOP 14.002 | PASS | PASS | 4.729 (4.73) |
| clang | south-48h-30min | 0 | 1965.450 | 95163 (95163) | 14.334 | ISOP 13.841 | PASS | PASS | 4.966 (4.97) |
| clang | spring45-20min | 0 | 1265.228 | 60435 (60435) | 14.344 | ISOP 13.792 | PASS | PASS | 4.842 (4.84) |
| clang | tropics-polluted-15min | 0 | 1586.414 | 75293 (75293) | 14.323 | ISOP 13.681 | PASS | PASS | 5.1 (5.10) |
| clang | clean-10min | 0 | 1808.889 | 84273 (84273) | 14.347 | ISOP 13.812 | PASS | PASS | 4.907 (4.91) |
| clang | diurnal-10min | 0 | 1956.949 | 90980 (90980) | 14.291 | ISOP 13.688 | PASS | PASS | 5.052 (5.05) |
| clang | diurnal-30min | 0 | 1016.829 | 48619 (48619) | 14.369 | ISOP 13.799 | PASS | PASS | 4.748 (4.75) |
| clang | diurnal-5min | 0 | 3186.686 | 145891 (145891) | 14.234 | ISOP 13.677 | PASS | PASS | 5.214 (5.21) |
| clang | polluted-10min | 0 | 2282.821 | 106393 (106393) | 14.044 | ISOP 13.417 | PASS | PASS | 5.398 (5.40) |
| gcc-15 | arctic-summer-10min | 0 | 2583.766 | 103211 (103211) | 14.29 | ISOP 13.799 | PASS | PASS | 5.045 (5.04) |
| gcc-15 | night-start-5min | 0 | 3410.644 | 132855 (132855) | 14.452 | ISOP 13.998 | PASS | PASS | 4.729 (4.73) |
| gcc-15 | south-48h-30min | 0 | 2326.867 | 95163 (95163) | 14.333 | ISOPAO2 13.726 | PASS | PASS | 4.966 (4.97) |
| gcc-15 | spring45-20min | 0 | 1492.616 | 60435 (60435) | 14.337 | ISOP 13.767 | PASS | PASS | 4.842 (4.84) |
| gcc-15 | tropics-polluted-15min | 0 | 1855.784 | 75293 (75293) | 14.336 | ISOP 13.641 | PASS | PASS | 5.1 (5.10) |
| gcc-15 | clean-10min | 0 | 2127.656 | 84273 (84273) | 14.354 | ISOP 13.835 | PASS | PASS | 4.907 (4.91) |
| gcc-15 | diurnal-10min | 0 | 2285.437 | 90980 (90980) | 14.296 | ISOP 13.704 | PASS | PASS | 5.052 (5.05) |
| gcc-15 | diurnal-30min | 0 | 1181.397 | 48619 (48619) | 14.372 | ISOP 13.803 | PASS | PASS | 4.748 (4.75) |
| gcc-15 | diurnal-5min | 0 | 3715.548 | 145891 (145891) | 14.228 | ISOPNO3 13.544 | PASS | PASS | 5.214 (5.21) |
| gcc-15 | polluted-10min | 0 | 2642.185 | 106393 (106393) | 14.038 | ISOP 13.414 | PASS | PASS | 5.398 (5.40) |

- clang: level 10 and the gate on every scenario: yes (10 scenarios)
- gcc-15: level 10 and the gate on every scenario: yes (10 scenarios)

## Build audit

`make -n -B cpu bench CC=clang`: exit 0; flags used: `--specialized -DNDEBUG -Ibuild -Iinclude -L~/homebrew/opt/libomp/lib -O3 -Wall -Wextra -Wpedantic -Xpreprocessor -c -fopenmp -isysroot -isystem -lm -lomp -p -std=c11`

Flags outside the rules (review): none

- `bin/test_cpu` links: /usr/lib/libSystem.B.dylib, ~/homebrew/opt/libomp/lib/libomp.dylib
- `bin/test_example_no_user` links: /usr/lib/libSystem.B.dylib, ~/homebrew/opt/libomp/lib/libomp.dylib
- `bin/test_example_user` links: /usr/lib/libSystem.B.dylib, ~/homebrew/opt/libomp/lib/libomp.dylib
- `bin/test_programs` links: /usr/lib/libSystem.B.dylib, ~/homebrew/opt/libomp/lib/libomp.dylib
- `bin/test_specialized_robertson` links: /usr/lib/libSystem.B.dylib, ~/homebrew/opt/libomp/lib/libomp.dylib
- `bin/test_specialized_ts1` links: /usr/lib/libSystem.B.dylib, ~/homebrew/opt/libomp/lib/libomp.dylib
- `bin/test_specialized_tutorial_no_user` links: /usr/lib/libSystem.B.dylib, ~/homebrew/opt/libomp/lib/libomp.dylib
- `bin/test_specialized_tutorial_user` links: /usr/lib/libSystem.B.dylib, ~/homebrew/opt/libomp/lib/libomp.dylib
- `bin/ts1_bench_cpu` links: /usr/lib/libSystem.B.dylib, ~/homebrew/opt/libomp/lib/libomp.dylib
- `bin/ts1_bench_specialized` links: /usr/lib/libSystem.B.dylib, ~/homebrew/opt/libomp/lib/libomp.dylib
- `bin/ts1_scenario_specialized` links: /usr/lib/libSystem.B.dylib, ~/homebrew/opt/libomp/lib/libomp.dylib

C sources scanned: 29 (tracked and generated, outside conda); lines matching a banned construct or starting a process: 0 (all in `logs/banned.txt`)


