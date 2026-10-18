# Verification of tier3/sol-max

Branch `tier3/sol/max`, commit c122d12 Finalize compiler gates, operating-point benchmarks, and research reports; compilers clang; OPTFLAGS `-O3 -DNDEBUG`; 2026-10-18 06:17.

## Operating points on the scenarios (--threads 10, 1,000 cells)

- `fast`: `bin/ts1_scenario --method rodas3 --option rtol=0.01 --option group=16 --option schedule=dynamic SCENARIO`
- `balanced`: `bin/ts1_scenario --method rodas3 --option rtol=0.005 --option group=16 --option schedule=dynamic SCENARIO`
- `accurate`: `bin/ts1_scenario --method rodas3 --option rtol=0.0005 --option group=16 --option schedule=dynamic SCENARIO`
- `recommended`: `bin/ts1_scenario --method rodas3 --option rtol=0.0003 --option group=16 --option schedule=dynamic SCENARIO`
- `precise`: `bin/ts1_scenario --method rodas3 --option rtol=0.0001 --option group=16 --option schedule=dynamic SCENARIO`

| Compiler | Point | Scenario | Exit | elapsed_ms | Digits | Worst species | Gate | Level 2 | Level 3 |
|---|---|---|---|---:|---:|---|---|---|---|
| clang | fast | clean-48h-5min | 0 | 859.688 | 2.303 | NO 1.409 | PASS | PASS | FAIL |
| clang | fast | polluted-48h-30min | 0 | 425.872 | 2.216 | ISOP 1.091 | PASS | PASS | FAIL |
| clang | fast | summer-15min | 0 | 257.721 | 2.341 | ISOP 1.22 | PASS | PASS | FAIL |
| clang | fast | tropics-20min | 0 | 282.832 | 2.183 | ISOP 1.157 | PASS | PASS | FAIL |
| clang | fast | winter60-10min | 0 | 352.838 | 2.277 | ISOP 1.195 | PASS | PASS | FAIL |
| clang | balanced | clean-48h-5min | 0 | 898.301 | 2.461 | PHENO 1.569 | PASS | PASS | FAIL |
| clang | balanced | polluted-48h-30min | 0 | 519.241 | 2.444 | ISOP 1.397 | PASS | PASS | FAIL |
| clang | balanced | summer-15min | 0 | 307.377 | 2.477 | ISOP 1.512 | PASS | PASS | FAIL |
| clang | balanced | tropics-20min | 0 | 329.262 | 2.337 | ISOP 1.388 | PASS | PASS | FAIL |
| clang | balanced | winter60-10min | 0 | 401.692 | 2.473 | ISOP 1.472 | PASS | PASS | FAIL |
| clang | accurate | clean-48h-5min | 0 | 1768.067 | 3.539 | ISOP 2.449 | PASS | PASS | PASS |
| clang | accurate | polluted-48h-30min | 0 | 1368.098 | 3.736 | ISOP 2.398 | PASS | PASS | PASS |
| clang | accurate | summer-15min | 0 | 868.055 | 3.25 | ISOP 2.431 | PASS | PASS | PASS |
| clang | accurate | tropics-20min | 0 | 813.352 | 3.609 | ISOP 2.388 | PASS | PASS | PASS |
| clang | accurate | winter60-10min | 0 | 880.955 | 3.664 | ISOP 2.554 | PASS | PASS | PASS |
| clang | recommended | clean-48h-5min | 0 | 2317.495 | 3.691 | ISOP 2.63 | PASS | PASS | PASS |
| clang | recommended | polluted-48h-30min | 0 | 1674.450 | 3.96 | ISOP 2.608 | PASS | PASS | PASS |
| clang | recommended | summer-15min | 0 | 1132.947 | 3.8 | ISOP 2.657 | PASS | PASS | PASS |
| clang | recommended | tropics-20min | 0 | 998.250 | 3.925 | ISOP 2.63 | PASS | PASS | PASS |
| clang | recommended | winter60-10min | 0 | 1058.776 | 3.948 | ISOP 2.8 | PASS | PASS | PASS |
| clang | precise | clean-48h-5min | 0 | 4300.419 | 4.314 | ISOP 3.109 | PASS | PASS | PASS |
| clang | precise | polluted-48h-30min | 0 | 2336.098 | 4.387 | ISOP 3.076 | PASS | PASS | PASS |
| clang | precise | summer-15min | 0 | 1751.050 | 4.423 | ISOP 3.181 | PASS | PASS | PASS |
| clang | precise | tropics-20min | 0 | 1466.256 | 4.368 | ISOP 3.118 | PASS | PASS | PASS |
| clang | precise | winter60-10min | 0 | 1502.846 | 4.439 | ISOP 3.306 | PASS | PASS | PASS |

- `fast`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: NO
- `balanced`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: NO
- `accurate`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
- `recommended`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
- `precise`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes

