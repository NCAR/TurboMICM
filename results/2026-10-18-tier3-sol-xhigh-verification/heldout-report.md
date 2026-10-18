# Verification of tier3/sol-xhigh

Branch `tier3/sol/xhigh`, commit 63d246f Complete integration research, work-precision reports and final validation audit; compilers clang; OPTFLAGS `-O3 -DNDEBUG`; 2026-10-18 04:20.

## Operating points on the scenarios (--threads 10, 1,000 cells)

- `L2`: `bin/ts1_scenario --method rodas31 --option rtol=0.002 --option atol=1e-18 --option group_size=32 --option schedule=2 --option schedule_chunk=1 SCENARIO`
- `L2_margin`: `bin/ts1_scenario --method rodas31 --option rtol=0.001 --option atol=1e-18 --option group_size=32 --option schedule=2 --option schedule_chunk=1 SCENARIO`
- `L3`: `bin/ts1_scenario --method rodas31 --option rtol=0.0003 --option atol=1e-16 --option group_size=32 --option schedule=2 --option schedule_chunk=1 SCENARIO`
- `recommended`: `bin/ts1_scenario --method rodas31 --option rtol=0.0003 --option atol=1e-18 --option group_size=32 --option schedule=2 --option schedule_chunk=1 SCENARIO`
- `L4`: `bin/ts1_scenario --method rodas4 --option rtol=0.0001 --option atol=1e-18 --option group_size=32 --option schedule=2 --option schedule_chunk=1 SCENARIO`

| Compiler | Point | Scenario | Exit | elapsed_ms | Digits | Worst species | Gate | Level 2 | Level 3 |
|---|---|---|---|---:|---:|---|---|---|---|
| clang | L2 | clean-48h-5min | 0 | 1213.189 | 2.474 | PHENO 1.583 | PASS | PASS | FAIL |
| clang | L2 | polluted-48h-30min | 0 | 576.061 | 2.455 | ISOP 1.71 | PASS | PASS | FAIL |
| clang | L2 | summer-15min | 0 | 341.654 | 2.607 | ISOP 1.646 | PASS | PASS | FAIL |
| clang | L2 | tropics-20min | 0 | 374.145 | 2.53 | ISOP 1.466 | PASS | PASS | FAIL |
| clang | L2 | winter60-10min | 0 | 497.264 | 2.672 | ISOP 1.613 | PASS | PASS | FAIL |
| clang | L2_margin | clean-48h-5min | 0 | 1411.457 | 3.586 | ISOP 2.449 | PASS | PASS | PASS |
| clang | L2_margin | polluted-48h-30min | 0 | 801.916 | 2.82 | ISOP 1.976 | PASS | PASS | FAIL |
| clang | L2_margin | summer-15min | 0 | 497.506 | 2.854 | ISOP 1.938 | PASS | PASS | FAIL |
| clang | L2_margin | tropics-20min | 0 | 501.234 | 2.876 | ISOP 1.703 | PASS | PASS | FAIL |
| clang | L2_margin | winter60-10min | 0 | 650.042 | 2.978 | ISOP 1.859 | PASS | PASS | FAIL |
| clang | L3 | clean-48h-5min | 0 | 1846.184 | 3.824 | ISOP 2.89 | PASS | PASS | PASS |
| clang | L3 | polluted-48h-30min | 0 | 1186.120 | 3.429 | ISOP 2.308 | PASS | PASS | PASS |
| clang | L3 | summer-15min | 0 | 782.862 | 3.387 | ISOP 2.502 | PASS | PASS | PASS |
| clang | L3 | tropics-20min | 0 | 732.803 | 3.618 | ISOP 2.361 | PASS | PASS | PASS |
| clang | L3 | winter60-10min | 0 | 857.100 | 3.708 | ISOP 2.575 | PASS | PASS | PASS |
| clang | recommended | clean-48h-5min | 0 | 2232.491 | 3.995 | ISOP 3.03 | PASS | PASS | PASS |
| clang | recommended | polluted-48h-30min | 0 | 1606.608 | 3.938 | ISOP 2.746 | PASS | PASS | PASS |
| clang | recommended | summer-15min | 0 | 1031.138 | 3.927 | ISOP 2.845 | PASS | PASS | PASS |
| clang | recommended | tropics-20min | 0 | 975.317 | 4.061 | ISOP 2.794 | PASS | PASS | PASS |
| clang | recommended | winter60-10min | 0 | 1078.639 | 4.131 | ISOP 2.98 | PASS | PASS | PASS |
| clang | L4 | clean-48h-5min | 0 | 8488.023 | 4.722 | ISOP 3.614 | PASS | PASS | PASS |
| clang | L4 | polluted-48h-30min | 0 | 3342.155 | 5.263 | ISOP 3.816 | PASS | PASS | PASS |
| clang | L4 | summer-15min | 0 | 2705.530 | 5.122 | ISOP 4.087 | PASS | PASS | PASS |
| clang | L4 | tropics-20min | 0 | 2022.558 | 5.231 | ISOP 4.09 | PASS | PASS | PASS |
| clang | L4 | winter60-10min | 0 | 2120.566 | 5.477 | ISOP 4.466 | PASS | PASS | PASS |

- `L2`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: NO
- `L2_margin`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: NO
- `L3`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
- `recommended`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
- `L4`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes

