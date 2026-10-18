# Verification of tier3/sol-xhigh

Branch `tier3/sol/xhigh`, commit 63d246f Complete integration research, work-precision reports and final validation audit; compilers clang, gcc-15; OPTFLAGS `-O3 -DNDEBUG`; 2026-10-18 04:17.

## Build

| Compiler | make env | make cpu bench |
|---|---|---|
| clang | exit 0 (11 s) | cpu exit 0 (1 s) bench exit 0 (0 s); 0 warnings |
| gcc-15 | exit 0 (12 s) | cpu exit 0 (2 s) bench exit 0 (0 s); 0 warnings |

## MICM's test cases

| Compiler | make test-cpu | PASS | FAIL | N/A | Ids missing | Ids repeated | Unknown ids |
|---|---|---:|---:|---:|---:|---:|---:|
| clang | exit 0 (4 s) | 675 | 0 | 5 | 0 | 0 | 0 |
| gcc-15 | exit 0 (6 s) | 675 | 0 | 5 | 0 | 0 | 0 |

## Tier 1 TS1 gate (10,000 cells, seed 2026, 100 steps of 30 s)

| Compiler | ts1_bench_cpu | elapsed_ms | compare_final.py (1e-10) |
|---|---|---:|---|
| clang | exit 0 | 8140.174 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 0 (rtol 1e-10, atol 1e-30) |
| gcc-15 | exit 0 | 9379.147 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 1.43e-13 (BCARYO2VBS, cell 4054) (rtol 1e-10, atol 1e-30) |

## Operating points on the scenarios (--threads 10, 1,000 cells)

- `L2`: `bin/ts1_scenario --method rodas31 --option rtol=0.002 --option atol=1e-18 --option group_size=32 --option schedule=2 --option schedule_chunk=1 SCENARIO`
- `L2_margin`: `bin/ts1_scenario --method rodas31 --option rtol=0.001 --option atol=1e-18 --option group_size=32 --option schedule=2 --option schedule_chunk=1 SCENARIO`
- `L3`: `bin/ts1_scenario --method rodas31 --option rtol=0.0003 --option atol=1e-16 --option group_size=32 --option schedule=2 --option schedule_chunk=1 SCENARIO`
- `recommended`: `bin/ts1_scenario --method rodas31 --option rtol=0.0003 --option atol=1e-18 --option group_size=32 --option schedule=2 --option schedule_chunk=1 SCENARIO`
- `L4`: `bin/ts1_scenario --method rodas4 --option rtol=0.0001 --option atol=1e-18 --option group_size=32 --option schedule=2 --option schedule_chunk=1 SCENARIO`

| Compiler | Point | Scenario | Exit | elapsed_ms | Digits | Worst species | Gate | Level 2 | Level 3 |
|---|---|---|---|---:|---:|---|---|---|---|
| clang | L2 | clean-10min | 0 | 440.174 | 2.805 | ISOP 1.696 | PASS | PASS | FAIL |
| clang | L2 | diurnal-10min | 0 | 489.904 | 2.842 | ISOP 1.708 | PASS | PASS | FAIL |
| clang | L2 | diurnal-30min | 0 | 353.812 | 2.541 | ISOP 1.445 | PASS | PASS | FAIL |
| clang | L2 | diurnal-5min | 0 | 685.028 | 3.3 | ISOP 2.297 | PASS | PASS | PASS |
| clang | L2 | polluted-10min | 0 | 572.838 | 2.843 | ISOP 1.933 | PASS | PASS | FAIL |
| clang | L2 | tier1-30s | 0 | 240.751 | 5.741 | ISOP 5.127 | PASS | PASS | PASS |
| clang | L2_margin | clean-10min | 0 | 584.844 | 3.146 | ISOP 2.208 | PASS | PASS | PASS |
| clang | L2_margin | diurnal-10min | 0 | 632.872 | 3.099 | ISOP 2.112 | PASS | PASS | PASS |
| clang | L2_margin | diurnal-30min | 0 | 466.127 | 2.799 | ISOP 1.689 | PASS | PASS | FAIL |
| clang | L2_margin | diurnal-5min | 0 | 806.304 | 3.444 | ISOP 2.335 | PASS | PASS | PASS |
| clang | L2_margin | polluted-10min | 0 | 752.736 | 3.04 | ISOP 2.144 | PASS | PASS | PASS |
| clang | L2_margin | tier1-30s | 0 | 228.088 | 5.852 | ISOP 5.127 | PASS | PASS | PASS |
| clang | L3 | clean-10min | 0 | 792.738 | 3.527 | ISOP 2.533 | PASS | PASS | PASS |
| clang | L3 | diurnal-10min | 0 | 879.022 | 3.501 | ISOP 2.491 | PASS | PASS | PASS |
| clang | L3 | diurnal-30min | 0 | 665.710 | 3.614 | ISOP 2.322 | PASS | PASS | PASS |
| clang | L3 | diurnal-5min | 0 | 1100.711 | 3.697 | ISOP 2.836 | PASS | PASS | PASS |
| clang | L3 | polluted-10min | 0 | 1037.131 | 3.341 | ISOP 2.472 | PASS | PASS | PASS |
| clang | L3 | tier1-30s | 0 | 254.279 | 5.935 | ISOP 5.128 | PASS | PASS | PASS |
| clang | recommended | clean-10min | 0 | 1094.896 | 3.84 | ISOP 2.929 | PASS | PASS | PASS |
| clang | recommended | diurnal-10min | 0 | 1155.193 | 3.794 | ISOP 2.899 | PASS | PASS | PASS |
| clang | recommended | diurnal-30min | 0 | 867.540 | 4.012 | ISOP 2.734 | PASS | PASS | PASS |
| clang | recommended | diurnal-5min | 0 | 1378.804 | 3.872 | ISOP 2.969 | PASS | PASS | PASS |
| clang | recommended | polluted-10min | 0 | 1513.035 | 3.705 | ISOP 2.86 | PASS | PASS | PASS |
| clang | recommended | tier1-30s | 0 | 310.820 | 5.959 | ISOP 5.14 | PASS | PASS | PASS |
| clang | L4 | clean-10min | 0 | 2991.398 | 4.882 | ISOP 3.797 | PASS | PASS | PASS |
| clang | L4 | diurnal-10min | 0 | 3079.571 | 4.9 | ISOP 3.764 | PASS | PASS | PASS |
| clang | L4 | diurnal-30min | 0 | 1621.564 | 5.224 | ISOP 3.877 | PASS | PASS | PASS |
| clang | L4 | diurnal-5min | 0 | 4500.320 | 4.657 | ISOP 3.614 | PASS | PASS | PASS |
| clang | L4 | polluted-10min | 0 | 3573.186 | 4.994 | ISOP 3.903 | PASS | PASS | PASS |
| clang | L4 | tier1-30s | 0 | 391.493 | 7.992 | C6H5O2 7.297 | PASS | PASS | PASS |

- `L2`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: NO
- `L2_margin`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: NO
- `L3`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
- `recommended`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
- `L4`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes

## Build audit

`make -n -B cpu bench CC=clang`: exit 0; flags used: `--file --prefix --symbol -DCK_DEFAULT_GROUP_SIZE=32 -DNDEBUG -Igenerated -Iinclude -Iinclude -Igenerated -L/opt/homebrew/opt/libomp/lib -L/opt/homebrew/opt/libomp/lib -lomp -O3 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Xpreprocessor -Xpreprocessor -fopenmp -isystem /opt/homebrew/opt/libomp/include -c -fopenmp -isystem -lm -lomp -p -s -std=c11 -x`

Flags outside the rules (review): none

- `bin/test_cpu` links: /opt/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_bench_cpu` links: /opt/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_scenario` links: /opt/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib

C sources scanned: 57 (tracked and generated, outside conda); lines matching a banned construct: 0 (all in `logs/banned.txt`)


