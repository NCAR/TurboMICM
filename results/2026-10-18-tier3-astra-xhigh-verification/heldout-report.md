# Verification of tier3/astra-xhigh

Branch `tier3/astra/xhigh`, commit db44838 Complete operating-point evaluation, work-precision figures and final research reports; compilers clang; OPTFLAGS `-O3 -DNDEBUG`; 2026-10-18 08:36.

## Operating points on the scenarios (--threads 10, 1,000 cells)

- `fast`: `bin/ts1_scenario --method rodas3 --option rtol=0.01 --option atol=1e-18 --option group-size=32 --option rate-cache=1 --option step-block=4 --option schedule=dynamic --option chunk=0 --option jacobian-age=1 --option controller=0 --option norm=1 --option warm-start=1 --option reuse-window=0 --option clip-stages=0 SCENARIO`
- `balanced`: `bin/ts1_scenario --method rodas3 --option rtol=0.003 --option atol=1e-18 --option group-size=32 --option rate-cache=1 --option step-block=4 --option schedule=dynamic --option chunk=0 --option jacobian-age=1 --option controller=0 --option norm=1 --option warm-start=1 --option reuse-window=0 --option clip-stages=0 SCENARIO`
- `accurate`: `bin/ts1_scenario --method extrap --option rtol=0.004 --option atol=1e-18 --option group-size=32 --option order=4 --option rate-cache=1 --option step-block=4 --option schedule=dynamic --option chunk=0 --option jacobian-age=1 --option controller=0 --option norm=1 --option warm-start=1 --option reuse-window=0 --option clip-stages=0 SCENARIO`
- `recommended`: `bin/ts1_scenario --method extrap --option rtol=0.002 --option atol=1e-18 --option group-size=32 --option order=4 --option rate-cache=1 --option step-block=4 --option schedule=dynamic --option chunk=0 --option jacobian-age=1 --option controller=0 --option norm=1 --option warm-start=1 --option reuse-window=0 --option clip-stages=0 SCENARIO`
- `strict`: `bin/ts1_scenario --method extrap --option rtol=0.001 --option atol=1e-18 --option group-size=32 --option order=4 --option rate-cache=1 --option step-block=4 --option schedule=dynamic --option chunk=0 --option jacobian-age=1 --option controller=0 --option norm=1 --option warm-start=1 --option reuse-window=0 --option clip-stages=0 SCENARIO`

| Compiler | Point | Scenario | Exit | elapsed_ms | Digits | Worst species | Gate | Level 2 | Level 3 |
|---|---|---|---|---:|---:|---|---|---|---|
| clang | fast | clean-48h-5min | 0 | 797.592163 | 2.333 | NO 1.449 | PASS | PASS | FAIL |
| clang | fast | polluted-48h-30min | 0 | 392.999887 | 2.21 | ISOP 1.133 | PASS | PASS | FAIL |
| clang | fast | summer-15min | 0 | 226.917744 | 2.303 | ISOP 1.236 | PASS | PASS | FAIL |
| clang | fast | tropics-20min | 0 | 237.386703 | 2.19 | ISOP 1.13 | PASS | PASS | FAIL |
| clang | fast | winter60-10min | 0 | 299.978971 | 2.27 | ISOP 1.174 | PASS | PASS | FAIL |
| clang | balanced | clean-48h-5min | 0 | 893.284559 | 2.463 | PHENO 1.569 | PASS | PASS | FAIL |
| clang | balanced | polluted-48h-30min | 0 | 613.862514 | 2.776 | ISOP 1.63 | PASS | PASS | FAIL |
| clang | balanced | summer-15min | 0 | 350.373268 | 2.533 | ISOP 1.652 | PASS | PASS | FAIL |
| clang | balanced | tropics-20min | 0 | 358.461857 | 2.538 | ISOP 1.629 | PASS | PASS | FAIL |
| clang | balanced | winter60-10min | 0 | 432.734728 | 2.674 | ISOP 1.679 | PASS | PASS | FAIL |
| clang | accurate | clean-48h-5min | 0 | 1779.496908 | 3.661 | HO2NO2 2.827 | PASS | PASS | PASS |
| clang | accurate | polluted-48h-30min | 0 | 753.181219 | 3.559 | ISOP 2.628 | PASS | PASS | PASS |
| clang | accurate | summer-15min | 0 | 495.028257 | 3.935 | OCLO 2.754 | PASS | PASS | PASS |
| clang | accurate | tropics-20min | 0 | 464.788198 | 3.653 | ISOP 2.705 | PASS | PASS | PASS |
| clang | accurate | winter60-10min | 0 | 673.997164 | 3.936 | ISOP 2.833 | PASS | PASS | PASS |
| clang | recommended | clean-48h-5min | 0 | 1888.580322 | 3.814 | HO2NO2 3.066 | PASS | PASS | PASS |
| clang | recommended | polluted-48h-30min | 0 | 959.461451 | 3.673 | ISOP 2.791 | PASS | PASS | PASS |
| clang | recommended | summer-15min | 0 | 558.098793 | 4.036 | ISOP 3.071 | PASS | PASS | PASS |
| clang | recommended | tropics-20min | 0 | 576.612949 | 4.017 | ISOP 2.807 | PASS | PASS | PASS |
| clang | recommended | winter60-10min | 0 | 774.138927 | 4.124 | ISOP 3.014 | PASS | PASS | PASS |
| clang | strict | clean-48h-5min | 0 | 2003.452778 | 3.861 | MALO2 3.14 | PASS | PASS | PASS |
| clang | strict | polluted-48h-30min | 0 | 1287.728071 | 4.228 | ISOP 2.984 | PASS | PASS | PASS |
| clang | strict | summer-15min | 0 | 704.458714 | 4.092 | ISOP 3.202 | PASS | PASS | PASS |
| clang | strict | tropics-20min | 0 | 741.677999 | 4.219 | ISOP 3.03 | PASS | PASS | PASS |
| clang | strict | winter60-10min | 0 | 942.986488 | 4.269 | ISOP 3.145 | PASS | PASS | PASS |

- `fast`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: NO
- `balanced`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: NO
- `accurate`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
- `recommended`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes
- `strict`: gate on every scenario: yes; level 2 on every scenario: yes; level 3 on every scenario: yes

