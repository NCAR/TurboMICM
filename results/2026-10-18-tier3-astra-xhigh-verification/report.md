# Verification of tier3/astra-xhigh

Branch `tier3/astra/xhigh`, commit db44838 Complete operating-point evaluation, work-precision figures and final research reports; compilers clang, gcc-15; OPTFLAGS `-O3 -DNDEBUG`; 2026-10-18 08:32.

## Build

| Compiler | make env | make cpu bench |
|---|---|---|
| clang | exit 0 (12 s) | cpu exit 0 (46 s) bench exit 0 (0 s); 0 warnings |
| gcc-15 | exit 0 (12 s) | cpu exit 0 (35 s) bench exit 0 (0 s); 0 warnings |

## MICM's test cases

| Compiler | make test-cpu | PASS | FAIL | N/A | Ids missing | Ids repeated | Unknown ids |
|---|---|---:|---:|---:|---:|---:|---:|
| clang | exit 0 (11 s) | 680 | 0 | 0 | 0 | 0 | 0 |
| gcc-15 | exit 0 (17 s) | 680 | 0 | 0 | 0 | 0 | 0 |

## Tier 1 TS1 gate (10,000 cells, seed 2026, 100 steps of 30 s)

| Compiler | ts1_bench_cpu | elapsed_ms | compare_final.py (1e-10) |
|---|---|---:|---|
| clang | exit 0 | 6871.623993 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.03e-13 (BCARY, cell 4179) (rtol 1e-10, atol 1e-30) |
| gcc-15 | exit 0 | 7713.189000 | PASS: PASS: 2110000 values in 10000 cells x 211 columns; largest relative difference 2.14e-13 (BCARY, cell 6848) (rtol 1e-10, atol 1e-30) |

## Operating points on the scenarios (--threads 10, 1,000 cells)

- `fast`: `bin/ts1_scenario --method rodas3 --option rtol=0.01 --option atol=1e-18 --option group-size=32 --option rate-cache=1 --option step-block=4 --option schedule=dynamic --option chunk=0 --option jacobian-age=1 --option controller=0 --option norm=1 --option warm-start=1 --option reuse-window=0 --option clip-stages=0 SCENARIO`
- `balanced`: `bin/ts1_scenario --method rodas3 --option rtol=0.003 --option atol=1e-18 --option group-size=32 --option rate-cache=1 --option step-block=4 --option schedule=dynamic --option chunk=0 --option jacobian-age=1 --option controller=0 --option norm=1 --option warm-start=1 --option reuse-window=0 --option clip-stages=0 SCENARIO`
- `accurate`: `bin/ts1_scenario --method extrap --option rtol=0.004 --option atol=1e-18 --option group-size=32 --option order=4 --option rate-cache=1 --option step-block=4 --option schedule=dynamic --option chunk=0 --option jacobian-age=1 --option controller=0 --option norm=1 --option warm-start=1 --option reuse-window=0 --option clip-stages=0 SCENARIO`
- `recommended`: `bin/ts1_scenario --method extrap --option rtol=0.002 --option atol=1e-18 --option group-size=32 --option order=4 --option rate-cache=1 --option step-block=4 --option schedule=dynamic --option chunk=0 --option jacobian-age=1 --option controller=0 --option norm=1 --option warm-start=1 --option reuse-window=0 --option clip-stages=0 SCENARIO`
- `strict`: `bin/ts1_scenario --method extrap --option rtol=0.001 --option atol=1e-18 --option group-size=32 --option order=4 --option rate-cache=1 --option step-block=4 --option schedule=dynamic --option chunk=0 --option jacobian-age=1 --option controller=0 --option norm=1 --option warm-start=1 --option reuse-window=0 --option clip-stages=0 SCENARIO`

| Compiler | Point | Scenario | Exit | elapsed_ms | Digits | Worst species | Gate | Level 2 | Level 3 |
|---|---|---|---|---:|---:|---|---|---|---|
| clang | fast | clean-10min | 0 | 269.258976 | 2.537 | ISOP 1.525 | PASS | PASS | FAIL |
| clang | fast | diurnal-10min | 0 | 282.574177 | 2.527 | ISOP 1.393 | PASS | PASS | FAIL |
| clang | fast | diurnal-30min | 0 | 210.593939 | 2.102 | ISOP 1.035 | PASS | PASS | FAIL |
| clang | fast | diurnal-5min | 0 | 427.569628 | 2.041 | BRO 1.153 | PASS | PASS | FAIL |
| clang | fast | polluted-10min | 0 | 336.705208 | 2.401 | ISOP 1.3 | PASS | PASS | FAIL |
| clang | fast | tier1-30s | 0 | 153.798103 | 5.029 | C6H5O2 4.31 | PASS | PASS | PASS |
| clang | balanced | clean-10min | 0 | 390.110016 | 2.765 | ISOP 1.792 | PASS | PASS | FAIL |
| clang | balanced | diurnal-10min | 0 | 410.676003 | 2.728 | ISOP 1.681 | PASS | PASS | FAIL |
| clang | balanced | diurnal-30min | 0 | 326.029301 | 2.739 | ISOP 1.646 | PASS | PASS | FAIL |
| clang | balanced | diurnal-5min | 0 | 510.091782 | 3.162 | ISOP 2.096 | PASS | PASS | PASS |
| clang | balanced | polluted-10min | 0 | 513.029337 | 2.645 | ISOP 1.729 | PASS | PASS | FAIL |
| clang | balanced | tier1-30s | 0 | 160.169125 | 5.486 | C6H5O2 4.836 | PASS | PASS | PASS |
| clang | accurate | clean-10min | 0 | 526.969910 | 3.201 | BRO 2.501 | PASS | PASS | PASS |
| clang | accurate | diurnal-10min | 0 | 573.259830 | 3.213 | OCLO 2.39 | PASS | PASS | PASS |
| clang | accurate | diurnal-30min | 0 | 427.609921 | 3.67 | ISOP 2.529 | PASS | PASS | PASS |
| clang | accurate | diurnal-5min | 0 | 1005.594015 | 3.448 | BRO 2.55 | PASS | PASS | PASS |
| clang | accurate | polluted-10min | 0 | 720.103979 | 3.927 | ISOP 2.794 | PASS | PASS | PASS |
| clang | accurate | tier1-30s | 0 | 354.341030 | 6.085 | MEKO2 5.497 | PASS | PASS | PASS |
| clang | recommended | clean-10min | 0 | 602.373838 | 4.243 | ISOP 3.084 | PASS | PASS | PASS |
| clang | recommended | diurnal-10min | 0 | 694.551468 | 4.196 | ISOP 3.015 | PASS | PASS | PASS |
| clang | recommended | diurnal-30min | 0 | 512.052298 | 3.824 | ISOP 2.689 | PASS | PASS | PASS |
| clang | recommended | diurnal-5min | 0 | 980.376244 | 3.448 | BRO 2.55 | PASS | PASS | PASS |
| clang | recommended | polluted-10min | 0 | 885.501146 | 4.137 | ISOP 2.959 | PASS | PASS | PASS |
| clang | recommended | tier1-30s | 0 | 360.485077 | 6.334 | MEKO2 5.814 | PASS | PASS | PASS |
| clang | strict | clean-10min | 0 | 785.354853 | 4.553 | ISOP 3.403 | PASS | PASS | PASS |
| clang | strict | diurnal-10min | 0 | 845.487118 | 4.562 | ISOP 3.334 | PASS | PASS | PASS |
| clang | strict | diurnal-30min | 0 | 702.108860 | 4.064 | ISOP 2.919 | PASS | PASS | PASS |
| clang | strict | diurnal-5min | 0 | 1083.094835 | 4.824 | ISOP 3.674 | PASS | PASS | PASS |
| clang | strict | polluted-10min | 0 | 1158.963203 | 4.578 | ISOP 3.412 | PASS | PASS | PASS |
| clang | strict | tier1-30s | 0 | 408.417225 | 6.663 | MEKO2 6.066 | PASS | PASS | PASS |

- `fast`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: NO
- `balanced`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: NO
- `accurate`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
- `recommended`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
- `strict`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes

## Build audit

`make -n -B cpu bench CC=clang`: exit 0; flags used: `--file --prefix --symbol -DNDEBUG -Iexamples -Iinclude -L/opt/homebrew/opt/libomp/lib -L/opt/homebrew/opt/libomp/lib -lomp -lm -O3 -O3 -DNDEBUG -std=c11 -Wall -Wextra -Wpedantic -Xpreprocessor -fopenmp -isystem /opt/homebrew/opt/libomp/include -Iinclude -Wall -Wextra -Wpedantic -Xpreprocessor -c -fopenmp -isystem -lm -lomp -p -s -std=c11`

Flags outside the rules (review): none

- `bin/test_cpu` links: /opt/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_bench_cpu` links: /opt/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/ts1_scenario` links: /opt/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib
- `bin/tutorial` links: /opt/homebrew/opt/libomp/lib/libomp.dylib, /usr/lib/libSystem.B.dylib

C sources scanned: 77 (tracked and generated, outside conda); lines matching a banned construct: 0 (all in `logs/banned.txt`)


