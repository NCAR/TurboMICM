# Verification of tier2/sol/high

Branch `tier2/sol/high`, commit 07848a8 Cover negative rate coefficients and record complete compiler gates and timings; compilers clang, gcc-15; OPTFLAGS `-O3 -DNDEBUG`; 18 threads; 2026-10-22 20:45.

## Build (fresh worktree)

| Compiler | make env | make cpu | make bench | Clean build | Warnings | Generated C |
|---|---|---|---|---:|---:|---|
| clang | exit 0 (12 s) | exit 0 (3.3 s) | exit 0 (13.8 s) | 17.2 s | 0 | 2 files, 27,951 lines, 1.14 MB |
| gcc-15 | exit 0 (10 s) | exit 0 (6.8 s) | exit 0 (17.0 s) | 23.8 s | 0 | 2 files, 27,951 lines, 1.14 MB |

clang: generated C: `build/ts1.c`, `build/ts1_specialized.c`

gcc-15: generated C: `build/ts1.c`, `build/ts1_specialized.c`

## MICM's test cases (the Tier 1 gate)

| Compiler | make test-cpu | PASS | FAIL | N/A | Ids missing | Ids repeated | Unknown ids |
|---|---|---:|---:|---:|---:|---:|---:|
| clang | exit 0 (12 s) | 675 | 0 | 5 | 0 | 0 | 0 |
| gcc-15 | exit 0 (17 s) | 675 | 0 | 5 | 0 | 0 | 0 |

## Tier 1 TS1 gate: bin/ts1_bench_cpu (10,000 cells, seed 2026, 100 steps of 30 s)

| Compiler | Exit | elapsed_ms | compare_final.py (1e-10) |
|---|---|---:|---|
| clang | 0 | 3637.786 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.03e-13 (BCARY, cell 4179) (rtol 1e-10, atol 1e-30) |
| gcc-15 | 0 | 4191.404 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.14e-13 (BCARY, cell 6848) (rtol 1e-10, atol 1e-30) |

## Specialized TS1 gate: bin/ts1_bench_specialized (10,000 cells, 100 steps of 30 s)

| Compiler | Input | Exit | elapsed_ms | Internal steps | compare_final.py (1e-10) |
|---|---|---|---:|---:|---|
| clang | seed 2026 | 0 | 2717.628 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.96e-13 (BCARY, cell 4179) (rtol 1e-10, atol 1e-30) |
| clang | input 1 | 0 | 2790.813 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.4e-13 (BCARYO2VBS, cell 9174) (rtol 1e-10, atol 1e-30) |
| gcc-15 | seed 2026 | 0 | 2986.186 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.75e-13 (BCARY, cell 4179) (rtol 1e-10, atol 1e-30) |
| gcc-15 | input 1 | 0 | 3084.121 | 281700 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.5e-13 (BCARYO2VBS, cell 2723) (rtol 1e-10, atol 1e-30) |

## The generator's tests: make test-specialized

| Compiler | Exit | Time | PASS lines | FAIL lines | Last lines |
|---|---|---:|---:|---:|---|
| clang | 0 | 5 s | 36 | 0 | PASS shapes: 32 cells, sparse plan, dynamic rates, default/tight trajectories / PASS shapes: 33 cells, sparse plan, dynamic rates, default/tight trajectories / All specialized mechanisms agree with th |
| gcc-15 | 0 | 7 s | 36 | 0 | PASS shapes: 32 cells, sparse plan, dynamic rates, default/tight trajectories / PASS shapes: 33 cells, sparse plan, dynamic rates, default/tight trajectories / All specialized mechanisms agree with th |

## Day-night scenarios (--threads 18, --rtol 1e-6 --atol 1e-15, 32-cell groups)

Agreement is score.py against MICM's result at the same settings, with `--level 10`; accuracy is score.py against the reference solution, for scale.

| Compiler | Scenario | Exit | elapsed_ms | Internal steps (MICM's) | Agreement with MICM | Worst species | Gate | Level 10 | Accuracy (MICM's) |
|---|---|---|---:|---|---:|---|---|---|---|
| clang | arctic-summer-10min | 0 | 1109.326 | 103211 (103211) | 14.298 | ISOP 13.843 | PASS | PASS | 5.045 (5.04) |
| clang | night-start-5min | 0 | 1455.194 | 132855 (132855) | 14.45 | ISOP 13.999 | PASS | PASS | 4.729 (4.73) |
| clang | south-48h-30min | 0 | 1020.142 | 95163 (95163) | 14.337 | ISOP 13.769 | PASS | PASS | 4.966 (4.97) |
| clang | spring45-20min | 0 | 646.555 | 60435 (60435) | 14.34 | ISOP 13.761 | PASS | PASS | 4.842 (4.84) |
| clang | tropics-polluted-15min | 0 | 812.008 | 75293 (75293) | 14.331 | ISOPO2VBS 13.434 | PASS | PASS | 5.1 (5.10) |
| clang | clean-10min | 0 | 920.965 | 84273 (84273) | 14.342 | ISOP 13.792 | PASS | PASS | 4.907 (4.91) |
| clang | diurnal-10min | 0 | 990.186 | 90980 (90980) | 14.285 | ISOP 13.655 | PASS | PASS | 5.052 (5.05) |
| clang | diurnal-30min | 0 | 528.748 | 48619 (48619) | 14.376 | ISOP 13.828 | PASS | PASS | 4.748 (4.75) |
| clang | diurnal-5min | 0 | 1662.389 | 145891 (145891) | 14.23 | ISOP 13.688 | PASS | PASS | 5.214 (5.21) |
| clang | polluted-10min | 0 | 1182.701 | 106393 (106393) | 14.047 | ISOPNO3 13.305 | PASS | PASS | 5.398 (5.40) |
| gcc-15 | arctic-summer-10min | 0 | 1223.657 | 103211 (103211) | 14.289 | ISOP 13.816 | PASS | PASS | 5.045 (5.04) |
| gcc-15 | night-start-5min | 0 | 1616.771 | 132855 (132855) | 14.442 | ISOP 13.99 | PASS | PASS | 4.729 (4.73) |
| gcc-15 | south-48h-30min | 0 | 1148.074 | 95163 (95163) | 14.335 | ISOP 13.772 | PASS | PASS | 4.966 (4.97) |
| gcc-15 | spring45-20min | 0 | 716.288 | 60435 (60435) | 14.345 | ISOP 13.763 | PASS | PASS | 4.842 (4.84) |
| gcc-15 | tropics-polluted-15min | 0 | 903.194 | 75293 (75293) | 14.332 | ISOP 13.51 | PASS | PASS | 5.1 (5.10) |
| gcc-15 | clean-10min | 0 | 1015.324 | 84273 (84273) | 14.335 | ISOP 13.811 | PASS | PASS | 4.907 (4.91) |
| gcc-15 | diurnal-10min | 0 | 1091.458 | 90980 (90980) | 14.292 | ISOP 13.686 | PASS | PASS | 5.052 (5.05) |
| gcc-15 | diurnal-30min | 0 | 580.456 | 48619 (48619) | 14.363 | ISOP 13.791 | PASS | PASS | 4.748 (4.75) |
| gcc-15 | diurnal-5min | 0 | 1786.361 | 145891 (145891) | 14.233 | ISOPNO3 13.537 | PASS | PASS | 5.214 (5.21) |
| gcc-15 | polluted-10min | 0 | 1295.930 | 106393 (106393) | 14.051 | ISOP 13.44 | PASS | PASS | 5.398 (5.40) |

- clang: level 10 and the gate on every scenario: yes (10 scenarios)
- gcc-15: level 10 and the gate on every scenario: yes (10 scenarios)

## Build audit

`make -n -B cpu bench CC=clang`: exit 0; flags used: `--batch --file --fixed-width --prefix --symbol --yes; -DCHEM_GROUP_SIZE=32 -DCHEM_SPECIALIZED -DNDEBUG -Iinclude -Isrc -L~/homebrew/opt/libomp/lib -L~/homebrew/opt/libomp/lib -lomp -lm -O3 -Wall -Wextra -Wpedantic -Xpreprocessor -c -d -dynamiclib -fPIC -fopenmp -isystem -lm -lomp -p -s -std=c11 -std=c11 -Wall -Wextra -Wpedantic -Iinclude -Isrc -DCHEM_GROUP_SIZE=32 -Xpreprocessor -fopenmp -isystem ~/homebrew/opt/libomp/include -O3 -DNDEBUG`

Flags outside the rules (review): none

- `bin/cpu_tutorials` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/rate_tutorial` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/test_rates` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_bench_cpu` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_bench_specialized` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_scenario_specialized` links: ~/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib

C sources scanned: 43 (tracked and generated, outside conda); lines matching a banned construct or starting a process: 0 (all in `logs/banned.txt`)


