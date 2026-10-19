# Verification of tier3/fable-max

Branch `tier3/fable/max`, commit ec34c99 Final: build times with both compilers, two gcc warnings fixed, the flags rule clears objects on a compiler change; compilers clang, gcc-15; OPTFLAGS `-O3 -DNDEBUG`; 2026-10-19 07:43.

## Operating points on the scenarios (--threads 10, 1,000 cells)

- `L2`: `bin/ts1_scenario --method rodas3 --option rtol=5e-3 --option atol=1e-20 --option quantize=1 --option controller=micm --option safety=0.9 SCENARIO`
- `balanced`: `bin/ts1_scenario --method rodas3 --option rtol=2e-3 --option atol=1e-20 --option quantize=1 --option controller=h211b --option safety=0.9 SCENARIO`
- `L3`: `bin/ts1_scenario --method extrap --option order=3 --option sequence=shifted --option rtol=1e-3 --option atol=1e-20 --option quantize=1 --option controller=micm --option safety=0.9 SCENARIO`
- `recommended`: `bin/ts1_scenario --method extrap --option order=4 --option sequence=shifted --option rtol=5e-4 --option atol=1e-22 --option quantize=1 --option controller=micm --option safety=0.9 SCENARIO`
- `accurate`: `bin/ts1_scenario --method extrap --option order=4 --option sequence=shifted --option rtol=2e-4 --option atol=1e-20 --option quantize=1 --option controller=micm --option safety=0.9 SCENARIO`

| Compiler | Point | Scenario | Exit | elapsed_ms | Digits | Worst species | Gate | Level 2 | Level 3 |
|---|---|---|---|---:|---:|---|---|---|---|
| clang | L2 | clean-48h-5min | 0 | 536.221 | 2.462 | PHENO 1.569 | PASS | PASS | FAIL |
| clang | L2 | polluted-48h-30min | 0 | 315.705 | 2.608 | ISOP 1.514 | PASS | PASS | FAIL |
| clang | L2 | summer-15min | 0 | 185.001 | 2.482 | ISOP 1.545 | PASS | PASS | FAIL |
| clang | L2 | tropics-20min | 0 | 187.643 | 2.304 | ISOP 1.465 | PASS | PASS | FAIL |
| clang | L2 | winter60-10min | 0 | 214.389 | 2.273 | ISOP 1.16 | PASS | PASS | FAIL |
| clang | balanced | clean-48h-5min | 0 | 853.867 | 2.706 | PHENO 1.788 | PASS | PASS | FAIL |
| clang | balanced | polluted-48h-30min | 0 | 638.433 | 3.347 | ISOP 2.264 | PASS | PASS | PASS |
| clang | balanced | summer-15min | 0 | 380.727 | 3.128 | ISOP 2.328 | PASS | PASS | PASS |
| clang | balanced | tropics-20min | 0 | 364.667 | 3.048 | ISOP 2.281 | PASS | PASS | PASS |
| clang | balanced | winter60-10min | 0 | 408.297 | 3.186 | ISOP 2.267 | PASS | PASS | PASS |
| clang | L3 | clean-48h-5min | 0 | 1233.573 | 3.684 | MALO2 2.953 | PASS | PASS | PASS |
| clang | L3 | polluted-48h-30min | 0 | 662.321 | 3.555 | ISOP 2.646 | PASS | PASS | PASS |
| clang | L3 | summer-15min | 0 | 388.468 | 3.738 | ISOP 2.667 | PASS | PASS | PASS |
| clang | L3 | tropics-20min | 0 | 406.428 | 3.661 | ISOP 2.471 | PASS | PASS | PASS |
| clang | L3 | winter60-10min | 0 | 515.625 | 3.725 | ISOP 2.639 | PASS | PASS | PASS |
| clang | recommended | clean-48h-5min | 0 | 1580.369 | 4.059 | NO 3.201 | PASS | PASS | PASS |
| clang | recommended | polluted-48h-30min | 0 | 657.000 | 3.92 | BR 3.209 | PASS | PASS | PASS |
| clang | recommended | summer-15min | 0 | 375.021 | 4.336 | ISOP 3.231 | PASS | PASS | PASS |
| clang | recommended | tropics-20min | 0 | 426.144 | 4.426 | ISOP 3.46 | PASS | PASS | PASS |
| clang | recommended | winter60-10min | 0 | 543.636 | 4.336 | ISOP 3.248 | PASS | PASS | PASS |
| clang | accurate | clean-48h-5min | 0 | 1675.731 | 4.266 | MALO2 3.537 | PASS | PASS | PASS |
| clang | accurate | polluted-48h-30min | 0 | 838.242 | 4.482 | ISOP 3.502 | PASS | PASS | PASS |
| clang | accurate | summer-15min | 0 | 497.183 | 4.46 | OCLO 3.418 | PASS | PASS | PASS |
| clang | accurate | tropics-20min | 0 | 528.013 | 4.55 | ISOP 3.472 | PASS | PASS | PASS |
| clang | accurate | winter60-10min | 0 | 672.423 | 4.661 | ISOP 3.598 | PASS | PASS | PASS |

- `L2`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: NO
- `balanced`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: NO
- `L3`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
- `recommended`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
- `accurate`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes

