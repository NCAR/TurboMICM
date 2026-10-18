# Verification of tier3/astra-max

Branch `tier3/astra/max`, commit 3cab602 Publish final work-precision results and reproducibility evidence; compilers clang, gcc-15; OPTFLAGS `-O3 -DNDEBUG`; 2026-10-18 11:30.

## Operating points on the scenarios (--threads 10, 1,000 cells)

- `economy`: `bin/ts1_scenario --method rodas3 --option rtol=0.003 --option atol=1e-16 --option h_max=0 --option group_size=32 SCENARIO`
- `fast3`: `bin/ts1_scenario --method rodas3 --option rtol=0.003 --option atol=1e-18 --option h_max=300 --option group_size=32 SCENARIO`
- `balanced`: `bin/ts1_scenario --method rodas3 --option rtol=0.003 --option atol=1e-18 --option h_max=200 --option group_size=32 SCENARIO`
- `conservative`: `bin/ts1_scenario --method rodas3 --option rtol=0.003 --option atol=1e-18 --option h_max=100 --option group_size=32 SCENARIO`

| Compiler | Point | Scenario | Exit | elapsed_ms | Digits | Worst species | Gate | Level 2 | Level 3 |
|---|---|---|---|---:|---:|---|---|---|---|
| clang | economy | clean-48h-5min | 0 | 1250.396 | 2.462 | PHENO 1.569 | PASS | PASS | FAIL |
| clang | economy | polluted-48h-30min | 0 | 693.537 | 2.399 | ISOP 1.302 | PASS | PASS | FAIL |
| clang | economy | summer-15min | 0 | 411.641 | 2.397 | ISOP 1.273 | PASS | PASS | FAIL |
| clang | economy | tropics-20min | 0 | 413.578 | 2.278 | ISOP 1.271 | PASS | PASS | FAIL |
| clang | economy | winter60-10min | 0 | 518.132 | 2.424 | ISOP 1.35 | PASS | PASS | FAIL |
| clang | fast3 | clean-48h-5min | 0 | 1357.605 | 2.463 | PHENO 1.569 | PASS | PASS | FAIL |
| clang | fast3 | polluted-48h-30min | 0 | 1240.098 | 3.364 | ISOP 2.219 | PASS | PASS | PASS |
| clang | fast3 | summer-15min | 0 | 653.072 | 3.305 | ISOP 2.324 | PASS | PASS | PASS |
| clang | fast3 | tropics-20min | 0 | 695.529 | 3.263 | ISOP 2.233 | PASS | PASS | PASS |
| clang | fast3 | winter60-10min | 0 | 763.253 | 2.94 | ISOP 1.866 | PASS | PASS | FAIL |
| clang | balanced | clean-48h-5min | 0 | 2022.466 | 3.294 | PHENO 2.618 | PASS | PASS | PASS |
| clang | balanced | polluted-48h-30min | 0 | 1564.116 | 3.653 | ISOP 2.699 | PASS | PASS | PASS |
| clang | balanced | summer-15min | 0 | 834.816 | 3.585 | ISOP 2.807 | PASS | PASS | PASS |
| clang | balanced | tropics-20min | 0 | 845.286 | 3.54 | ISOP 2.666 | PASS | PASS | PASS |
| clang | balanced | winter60-10min | 0 | 904.828 | 3.273 | ISOP 2.263 | PASS | PASS | PASS |
| clang | conservative | clean-48h-5min | 0 | 2701.883 | 3.942 | HO2NO2 3.26 | PASS | PASS | PASS |
| clang | conservative | polluted-48h-30min | 0 | 2673.102 | 4.074 | MPAN 3.483 | PASS | PASS | PASS |
| clang | conservative | summer-15min | 0 | 1312.495 | 4.05 | HYDRALD 3.479 | PASS | PASS | PASS |
| clang | conservative | tropics-20min | 0 | 1369.855 | 3.96 | HYDRALD 3.384 | PASS | PASS | PASS |
| clang | conservative | winter60-10min | 0 | 1423.560 | 3.775 | ISOP 2.951 | PASS | PASS | PASS |

- `economy`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: NO
- `fast3`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: NO
- `balanced`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
- `conservative`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes

