# Verification of tier3/sol-max

Branch `tier3/sol/max`, commit c122d12 Finalize compiler gates, operating-point benchmarks, and research reports; compilers clang, gcc-15; OPTFLAGS `-O3 -DNDEBUG`; 2026-10-18 06:15.

## Build

| Compiler | make env | make cpu bench |
|---|---|---|
| clang | exit 0 (11 s) | cpu exit 0 (7 s) bench exit 0 (0 s); 0 warnings |
| gcc-15 | exit 0 (10 s) | cpu exit 0 (12 s) bench exit 0 (0 s); 0 warnings |

## MICM's test cases

| Compiler | make test-cpu | PASS | FAIL | N/A | Ids missing | Ids repeated | Unknown ids |
|---|---|---:|---:|---:|---:|---:|---:|
| clang | exit 0 (5 s) | 675 | 0 | 5 | 0 | 0 | 0 |
| gcc-15 | exit 0 (4 s) | 675 | 0 | 5 | 0 | 0 | 0 |

## Tier 1 TS1 gate (10,000 cells, seed 2026, 100 steps of 30 s)

| Compiler | ts1_bench_cpu | elapsed_ms | compare_final.py (1e-10) |
|---|---|---:|---|
| clang | exit 0 | 6830.027 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.26e-13 (BCARY, cell 4179) (rtol 1e-10, atol 1e-30) |
| gcc-15 | exit 0 | 7394.601 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.16e-13 (BCARYO2VBS, cell 6455) (rtol 1e-10, atol 1e-30) |

## Operating points on the scenarios (--threads 10, 1,000 cells)

- `fast`: `bin/ts1_scenario --method rodas3 --option rtol=0.01 --option group=16 --option schedule=dynamic SCENARIO`
- `balanced`: `bin/ts1_scenario --method rodas3 --option rtol=0.005 --option group=16 --option schedule=dynamic SCENARIO`
- `accurate`: `bin/ts1_scenario --method rodas3 --option rtol=0.0005 --option group=16 --option schedule=dynamic SCENARIO`
- `recommended`: `bin/ts1_scenario --method rodas3 --option rtol=0.0003 --option group=16 --option schedule=dynamic SCENARIO`
- `precise`: `bin/ts1_scenario --method rodas3 --option rtol=0.0001 --option group=16 --option schedule=dynamic SCENARIO`

| Compiler | Point | Scenario | Exit | elapsed_ms | Digits | Worst species | Gate | Level 2 | Level 3 |
|---|---|---|---|---:|---:|---|---|---|---|
| clang | fast | clean-10min | 0 | 310.861 | 2.413 | ISOP 1.528 | PASS | PASS | FAIL |
| clang | fast | diurnal-10min | 0 | 330.629 | 2.509 | ISOP 1.401 | PASS | PASS | FAIL |
| clang | fast | diurnal-30min | 0 | 254.652 | 2.122 | ISOP 1.064 | PASS | PASS | FAIL |
| clang | fast | diurnal-5min | 0 | 481.060 | 2.041 | BRO 1.153 | PASS | PASS | FAIL |
| clang | fast | polluted-10min | 0 | 376.999 | 2.401 | ISOP 1.305 | PASS | PASS | FAIL |
| clang | fast | tier1-30s | 0 | 159.751 | 5.062 | C6H5O2 4.374 | PASS | PASS | PASS |
| clang | balanced | clean-10min | 0 | 355.862 | 2.644 | ISOP 1.603 | PASS | PASS | FAIL |
| clang | balanced | diurnal-10min | 0 | 383.162 | 2.628 | ISOP 1.523 | PASS | PASS | FAIL |
| clang | balanced | diurnal-30min | 0 | 295.998 | 2.328 | ISOP 1.386 | PASS | PASS | FAIL |
| clang | balanced | diurnal-5min | 0 | 501.037 | 2.29 | BRO 1.413 | PASS | PASS | FAIL |
| clang | balanced | polluted-10min | 0 | 462.279 | 2.566 | ISOP 1.578 | PASS | PASS | FAIL |
| clang | balanced | tier1-30s | 0 | 164.599 | 5.322 | C6H5O2 4.682 | PASS | PASS | PASS |
| clang | accurate | clean-10min | 0 | 878.561 | 3.299 | ISOP 2.398 | PASS | PASS | PASS |
| clang | accurate | diurnal-10min | 0 | 938.624 | 3.273 | ISOP 2.341 | PASS | PASS | PASS |
| clang | accurate | diurnal-30min | 0 | 708.979 | 3.701 | ISOP 2.404 | PASS | PASS | PASS |
| clang | accurate | diurnal-5min | 0 | 1074.620 | 3.414 | ISOP 2.388 | PASS | PASS | PASS |
| clang | accurate | polluted-10min | 0 | 1252.526 | 3.344 | ISOP 2.425 | PASS | PASS | PASS |
| clang | accurate | tier1-30s | 0 | 205.052 | 5.822 | ISOP 4.995 | PASS | PASS | PASS |
| clang | recommended | clean-10min | 0 | 1150.955 | 3.665 | ISOP 2.696 | PASS | PASS | PASS |
| clang | recommended | diurnal-10min | 0 | 1209.909 | 3.626 | ISOP 2.593 | PASS | PASS | PASS |
| clang | recommended | diurnal-30min | 0 | 865.257 | 3.884 | ISOP 2.632 | PASS | PASS | PASS |
| clang | recommended | diurnal-5min | 0 | 1384.917 | 3.571 | ISOP 2.554 | PASS | PASS | PASS |
| clang | recommended | polluted-10min | 0 | 1597.177 | 3.759 | ISOP 2.663 | PASS | PASS | PASS |
| clang | recommended | tier1-30s | 0 | 227.166 | 5.842 | ISOP 5.003 | PASS | PASS | PASS |
| clang | precise | clean-10min | 0 | 1872.346 | 4.422 | ISOP 3.261 | PASS | PASS | PASS |
| clang | precise | diurnal-10min | 0 | 1891.649 | 4.458 | ISOP 3.161 | PASS | PASS | PASS |
| clang | precise | diurnal-30min | 0 | 1243.733 | 4.296 | ISOP 3.08 | PASS | PASS | PASS |
| clang | precise | diurnal-5min | 0 | 2397.262 | 4.25 | ISOP 3.069 | PASS | PASS | PASS |
| clang | precise | polluted-10min | 0 | 2297.091 | 4.428 | ISOP 3.217 | PASS | PASS | PASS |
| clang | precise | tier1-30s | 0 | 267.327 | 5.862 | ISOP 5.017 | PASS | PASS | PASS |

- `fast`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: NO
- `balanced`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: NO
- `accurate`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
- `recommended`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
- `precise`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes

## Build audit

`make -n -B cpu bench CC=clang`: exit 0; flags used: `--file --kernel-chunk --kernel-chunk 32 --name --prefix --yes -DCHEM_GROUP_SIZE=32 -DNDEBUG -Igenerated -Iinclude -L/opt/homebrew/opt/libomp/lib -O3 -Wall -Wextra -Wpedantic -Xpreprocessor -c -dynamiclib -f -fPIC -fopenmp -isystem -lm -lomp -p -s -std=c11`

Flags outside the rules (review): none

- `bin/libchem.dylib` links: build/clang/libchem.dylib, /opt/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_bench_cpu` links: /opt/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_scenario` links: /opt/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/tutorials` links: /opt/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib

C sources scanned: 31 (tracked and generated, outside conda); lines matching a banned construct: 0 (all in `logs/banned.txt`)


