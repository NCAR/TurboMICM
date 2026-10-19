# Verification of tier3/fable-max

Branch `tier3/fable/max`, commit cb0f8a3 L3r: ROS3 under H211b at 8e-3 / 1e-20 replaces the balanced point as a level-3 point (S15, S15b, S16); the judgment in DECISIONS.md; plots; compilers clang, gcc-15; OPTFLAGS `-O3 -DNDEBUG`; 2026-10-19 08:12.

## Operating points on the scenarios (--threads 10, 1,000 cells)

- `L2`: `bin/ts1_scenario --method rodas3 --option rtol=5e-3 --option atol=1e-20 --option quantize=1 --option controller=micm --option safety=0.9 SCENARIO`
- `L3r`: `bin/ts1_scenario --method ros3 --option rtol=8e-3 --option atol=1e-20 --option quantize=1 --option controller=h211b --option safety=0.9 SCENARIO`
- `L3`: `bin/ts1_scenario --method extrap --option order=3 --option sequence=shifted --option rtol=1e-3 --option atol=1e-20 --option quantize=1 --option controller=micm --option safety=0.9 SCENARIO`
- `recommended`: `bin/ts1_scenario --method extrap --option order=4 --option sequence=shifted --option rtol=5e-4 --option atol=1e-22 --option quantize=1 --option controller=micm --option safety=0.9 SCENARIO`
- `accurate`: `bin/ts1_scenario --method extrap --option order=4 --option sequence=shifted --option rtol=2e-4 --option atol=1e-20 --option quantize=1 --option controller=micm --option safety=0.9 SCENARIO`

| Compiler | Point | Scenario | Exit | elapsed_ms | Digits | Worst species | Gate | Level 2 | Level 3 |
|---|---|---|---|---:|---:|---|---|---|---|
| clang | L2 | clean-48h-5min | 0 | 532.592 | 2.462 | PHENO 1.569 | PASS | PASS | FAIL |
| clang | L2 | polluted-48h-30min | 0 | 319.955 | 2.608 | ISOP 1.514 | PASS | PASS | FAIL |
| clang | L2 | summer-15min | 0 | 179.233 | 2.482 | ISOP 1.545 | PASS | PASS | FAIL |
| clang | L2 | tropics-20min | 0 | 188.494 | 2.304 | ISOP 1.465 | PASS | PASS | FAIL |
| clang | L2 | winter60-10min | 0 | 209.190 | 2.273 | ISOP 1.16 | PASS | PASS | FAIL |
| clang | L3r | clean-48h-5min | 0 | 967.273 | 3.593 | HYDRALD 2.827 | PASS | PASS | PASS |
| clang | L3r | polluted-48h-30min | 0 | 730.098 | 3.522 | ISOP 2.326 | PASS | PASS | PASS |
| clang | L3r | summer-15min | 0 | 420.494 | 3.308 | HYDRALD 2.597 | PASS | PASS | PASS |
| clang | L3r | tropics-20min | 0 | 399.871 | 3.183 | HYDRALD 2.421 | PASS | PASS | PASS |
| clang | L3r | winter60-10min | 0 | 424.454 | 3.418 | ISOP 2.442 | PASS | PASS | PASS |
| clang | L3 | clean-48h-5min | 0 | 1229.263 | 3.684 | MALO2 2.953 | PASS | PASS | PASS |
| clang | L3 | polluted-48h-30min | 0 | 661.572 | 3.555 | ISOP 2.646 | PASS | PASS | PASS |
| clang | L3 | summer-15min | 0 | 387.313 | 3.738 | ISOP 2.667 | PASS | PASS | PASS |
| clang | L3 | tropics-20min | 0 | 397.381 | 3.661 | ISOP 2.471 | PASS | PASS | PASS |
| clang | L3 | winter60-10min | 0 | 518.799 | 3.725 | ISOP 2.639 | PASS | PASS | PASS |
| clang | recommended | clean-48h-5min | 0 | 1593.306 | 4.059 | NO 3.201 | PASS | PASS | PASS |
| clang | recommended | polluted-48h-30min | 0 | 661.207 | 3.92 | BR 3.209 | PASS | PASS | PASS |
| clang | recommended | summer-15min | 0 | 368.178 | 4.336 | ISOP 3.231 | PASS | PASS | PASS |
| clang | recommended | tropics-20min | 0 | 423.322 | 4.426 | ISOP 3.46 | PASS | PASS | PASS |
| clang | recommended | winter60-10min | 0 | 536.939 | 4.336 | ISOP 3.248 | PASS | PASS | PASS |
| clang | accurate | clean-48h-5min | 0 | 1685.963 | 4.266 | MALO2 3.537 | PASS | PASS | PASS |
| clang | accurate | polluted-48h-30min | 0 | 823.587 | 4.482 | ISOP 3.502 | PASS | PASS | PASS |
| clang | accurate | summer-15min | 0 | 502.805 | 4.46 | OCLO 3.418 | PASS | PASS | PASS |
| clang | accurate | tropics-20min | 0 | 515.733 | 4.55 | ISOP 3.472 | PASS | PASS | PASS |
| clang | accurate | winter60-10min | 0 | 665.117 | 4.661 | ISOP 3.598 | PASS | PASS | PASS |

- `L2`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: NO
- `L3r`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
- `L3`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
- `recommended`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
- `accurate`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes

