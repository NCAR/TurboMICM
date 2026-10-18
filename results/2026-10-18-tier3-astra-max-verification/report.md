# Verification of tier3/astra-max

Branch `tier3/astra/max`, commit 3cab602 Publish final work-precision results and reproducibility evidence; compilers clang, gcc-15; OPTFLAGS `-O3 -DNDEBUG`; 2026-10-18 11:27.

## Build

| Compiler | make env | make cpu bench |
|---|---|---|
| clang | exit 0 (12 s) | cpu exit 0 (2 s) bench exit 0 (0 s); 1 warnings |
| gcc-15 | exit 0 (10 s) | cpu exit 0 (4 s) bench exit 0 (0 s); 0 warnings |

## MICM's test cases

| Compiler | make test-cpu | PASS | FAIL | N/A | Ids missing | Ids repeated | Unknown ids |
|---|---|---:|---:|---:|---:|---:|---:|
| clang | exit 0 (6 s) | 680 | 0 | 0 | 0 | 0 | 0 |
| gcc-15 | exit 0 (5 s) | 680 | 0 | 0 | 0 | 0 | 0 |

## Tier 1 TS1 gate (10,000 cells, seed 2026, 100 steps of 30 s)

| Compiler | ts1_bench_cpu | elapsed_ms | compare_final.py (1e-10) |
|---|---|---:|---|
| clang | exit 0 | 8057.091 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.96e-13 (BCARY, cell 4179) (rtol 1e-10, atol 1e-30) |
| gcc-15 | exit 0 | 9473.133 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.75e-13 (BCARY, cell 4179) (rtol 1e-10, atol 1e-30) |

## Operating points on the scenarios (--threads 10, 1,000 cells)

- `economy`: `bin/ts1_scenario --method rodas3 --option rtol=0.003 --option atol=1e-16 --option h_max=0 --option group_size=32 SCENARIO`
- `fast3`: `bin/ts1_scenario --method rodas3 --option rtol=0.003 --option atol=1e-18 --option h_max=300 --option group_size=32 SCENARIO`
- `balanced`: `bin/ts1_scenario --method rodas3 --option rtol=0.003 --option atol=1e-18 --option h_max=200 --option group_size=32 SCENARIO`
- `conservative`: `bin/ts1_scenario --method rodas3 --option rtol=0.003 --option atol=1e-18 --option h_max=100 --option group_size=32 SCENARIO`

| Compiler | Point | Scenario | Exit | elapsed_ms | Digits | Worst species | Gate | Level 2 | Level 3 |
|---|---|---|---|---:|---:|---|---|---|---|
| clang | economy | clean-10min | 0 | 507.559 | 2.622 | ISOP 1.562 | PASS | PASS | FAIL |
| clang | economy | diurnal-10min | 0 | 521.676 | 2.59 | ISOP 1.446 | PASS | PASS | FAIL |
| clang | economy | diurnal-30min | 0 | 384.910 | 2.266 | ISOP 1.271 | PASS | PASS | FAIL |
| clang | economy | diurnal-5min | 0 | 717.493 | 2.794 | BRO 1.949 | PASS | PASS | FAIL |
| clang | economy | polluted-10min | 0 | 619.880 | 2.456 | ISOP 1.362 | PASS | PASS | FAIL |
| clang | economy | tier1-30s | 0 | 225.048 | 5.357 | C6H5O2 4.667 | PASS | PASS | PASS |
| clang | fast3 | clean-10min | 0 | 711.047 | 3.224 | ISOP 2.259 | PASS | PASS | PASS |
| clang | fast3 | diurnal-10min | 0 | 760.065 | 3.17 | ISOP 2.157 | PASS | PASS | PASS |
| clang | fast3 | diurnal-30min | 0 | 673.965 | 3.343 | ISOP 2.299 | PASS | PASS | PASS |
| clang | fast3 | diurnal-5min | 0 | 823.774 | 3.139 | ISOP 2.099 | PASS | PASS | PASS |
| clang | fast3 | polluted-10min | 0 | 894.609 | 3.086 | ISOP 2.06 | PASS | PASS | PASS |
| clang | fast3 | tier1-30s | 0 | 243.687 | 5.463 | C6H5O2 4.813 | PASS | PASS | PASS |
| clang | balanced | clean-10min | 0 | 845.005 | 3.553 | ISOP 2.667 | PASS | PASS | PASS |
| clang | balanced | diurnal-10min | 0 | 887.174 | 3.466 | ISOP 2.607 | PASS | PASS | PASS |
| clang | balanced | diurnal-30min | 0 | 829.950 | 3.612 | ISOP 2.721 | PASS | PASS | PASS |
| clang | balanced | diurnal-5min | 0 | 1203.338 | 3.449 | ISOP 2.58 | PASS | PASS | PASS |
| clang | balanced | polluted-10min | 0 | 1021.064 | 3.433 | ISOP 2.509 | PASS | PASS | PASS |
| clang | balanced | tier1-30s | 0 | 246.469 | 5.463 | C6H5O2 4.813 | PASS | PASS | PASS |
| clang | conservative | clean-10min | 0 | 1367.481 | 4.125 | ISOP 3.349 | PASS | PASS | PASS |
| clang | conservative | diurnal-10min | 0 | 1374.311 | 3.975 | HYDRALD 3.395 | PASS | PASS | PASS |
| clang | conservative | diurnal-30min | 0 | 1382.509 | 4.008 | HYDRALD 3.433 | PASS | PASS | PASS |
| clang | conservative | diurnal-5min | 0 | 1457.044 | 3.867 | BRO 3.065 | PASS | PASS | PASS |
| clang | conservative | polluted-10min | 0 | 1515.411 | 3.947 | HOBR 3.199 | PASS | PASS | PASS |
| clang | conservative | tier1-30s | 0 | 245.954 | 5.463 | C6H5O2 4.813 | PASS | PASS | PASS |

- `economy`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: NO
- `fast3`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
- `balanced`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
- `conservative`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes

## Build audit

`make -n -B cpu bench CC=clang`: exit 0; flags used: `--file --name --prefix --yes -DCH_DEFAULT_GROUP_SIZE=32 -DNDEBUG -Ibuild/generated -Iinclude -L/opt/homebrew/opt/libomp/lib -O3 -Wall -Wextra -Wpedantic -Xpreprocessor -c -dynamiclib -fPIC -fopenmp -isystem -lm -lomp -p -s -std=c11`

Flags outside the rules (review): none

- `bin/test_native_cpu` links: /opt/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_bench_cpu` links: /opt/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_scenario` links: /opt/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/tutorial_cpu` links: /opt/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib

C sources scanned: 25 (tracked and generated, outside conda); lines matching a banned construct: 0 (all in `logs/banned.txt`)


