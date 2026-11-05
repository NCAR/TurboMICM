# Agreement of the C library with the archived Python solver

`agree.py`, 4 November 2026, one thread, the 64 cells of each scenario in two
groups of 32 (the archive's `data/scenarios/*`), against the rows of the
archive's `results/ts1/runs.csv` (`alg` at rtol 1e-2 and 5e-3 with
atol = rtol x 1e-14; `restart` at 3e-2 and 3e-4 with atol 1e-15). Steps,
rejections, forcing evaluations, factorizations, and solves are identical in
all 44 runs; digits agree to the scorer's three decimals. The library's
Jacobian count exceeds the Python count by the number of rejections, since
its in-place factorization re-evaluates the Jacobian after a rejection (the
prototype keeps it); the paper notes this convention.

| Configuration | Scenario | Steps | Rejected | Forcings | Solves | Digits C | Digits Python |
|---|---|---:|---:|---:|---:|---:|---:|
| retaining 1e-2 | tier1-30s | 269 | 0 | 538 | 807 | 5.820 | 5.8197 |
| retaining 1e-2 | diurnal-30min | 699 | 96 | 1302 | 2097 | 2.699 | 2.6990 |
| retaining 1e-2 | diurnal-10min | 808 | 74 | 1542 | 2424 | 3.360 | 3.3602 |
| retaining 1e-2 | diurnal-5min | 969 | 44 | 1894 | 2907 | 3.440 | 3.4398 |
| retaining 1e-2 | polluted-10min | 990 | 112 | 1868 | 2970 | 3.255 | 3.2554 |
| retaining 1e-2 | clean-10min | 765 | 64 | 1466 | 2295 | 3.448 | 3.4483 |
| retaining 1e-2 | summer-15min | 724 | 88 | 1360 | 2172 | 3.102 | 3.1019 |
| retaining 1e-2 | tropics-20min | 725 | 98 | 1352 | 2175 | 2.900 | 2.9000 |
| retaining 1e-2 | winter60-10min | 823 | 80 | 1566 | 2469 | 3.146 | 3.1463 |
| retaining 1e-2 | polluted-48h-30min | 1295 | 204 | 2386 | 3885 | 2.875 | 2.8751 |
| retaining 1e-2 | clean-48h-5min | 1817 | 80 | 3554 | 5451 | 3.496 | 3.4955 |
| retaining 5e-3 | tier1-30s | 282 | 0 | 564 | 846 | 5.840 | 5.8395 |
| retaining 5e-3 | diurnal-30min | 1102 | 176 | 2028 | 3306 | 3.286 | 3.2862 |
| retaining 5e-3 | diurnal-10min | 1168 | 153 | 2183 | 3504 | 3.560 | 3.5595 |
| retaining 5e-3 | diurnal-5min | 1269 | 108 | 2430 | 3807 | 3.508 | 3.5081 |
| retaining 5e-3 | polluted-10min | 1623 | 250 | 2996 | 4869 | 3.422 | 3.4215 |
| retaining 5e-3 | clean-10min | 1092 | 141 | 2043 | 3276 | 3.567 | 3.5665 |
| retaining 5e-3 | summer-15min | 1130 | 178 | 2082 | 3390 | 3.313 | 3.3132 |
| retaining 5e-3 | tropics-20min | 1150 | 186 | 2114 | 3450 | 3.206 | 3.2059 |
| retaining 5e-3 | winter60-10min | 1226 | 165 | 2287 | 3678 | 3.401 | 3.4012 |
| retaining 5e-3 | polluted-48h-30min | 2306 | 411 | 4201 | 6918 | 3.421 | 3.4208 |
| retaining 5e-3 | clean-48h-5min | 2325 | 187 | 4463 | 6975 | 3.582 | 3.5822 |
| restarting 3e-2 | tier1-30s | 1826 | 0 | 3652 | 5478 | 5.382 | 5.3821 |
| restarting 3e-2 | diurnal-30min | 989 | 0 | 1978 | 2967 | 2.316 | 2.3165 |
| restarting 3e-2 | diurnal-10min | 2732 | 0 | 5464 | 8196 | 3.396 | 3.3962 |
| restarting 3e-2 | diurnal-5min | 5291 | 0 | 10582 | 15873 | 3.314 | 3.3139 |
| restarting 3e-2 | polluted-10min | 2783 | 0 | 5566 | 8349 | 3.185 | 3.1855 |
| restarting 3e-2 | clean-10min | 2718 | 0 | 5436 | 8154 | 3.369 | 3.3686 |
| restarting 3e-2 | summer-15min | 1868 | 0 | 3736 | 5604 | 3.150 | 3.1502 |
| restarting 3e-2 | tropics-20min | 1429 | 0 | 2858 | 4287 | 2.628 | 2.6280 |
| restarting 3e-2 | winter60-10min | 2711 | 0 | 5422 | 8133 | 3.030 | 3.0300 |
| restarting 3e-2 | polluted-48h-30min | 1953 | 0 | 3906 | 5859 | 2.379 | 2.3792 |
| restarting 3e-2 | clean-48h-5min | 10524 | 0 | 21048 | 31572 | 3.465 | 3.4645 |
| restarting 3e-4 | tier1-30s | 1857 | 0 | 3714 | 5571 | 6.479 | 6.4795 |
| restarting 3e-4 | diurnal-30min | 1484 | 3 | 2965 | 4452 | 3.336 | 3.3359 |
| restarting 3e-4 | diurnal-10min | 3422 | 0 | 6844 | 10266 | 4.029 | 4.0287 |
| restarting 3e-4 | diurnal-5min | 6173 | 0 | 12346 | 18519 | 4.363 | 4.3632 |
| restarting 3e-4 | polluted-10min | 3609 | 0 | 7218 | 10827 | 3.867 | 3.8670 |
| restarting 3e-4 | clean-10min | 3374 | 3 | 6745 | 10122 | 3.981 | 3.9815 |
| restarting 3e-4 | summer-15min | 2511 | 0 | 5022 | 7533 | 3.744 | 3.7445 |
| restarting 3e-4 | tropics-20min | 1983 | 0 | 3966 | 5949 | 3.426 | 3.4256 |
| restarting 3e-4 | winter60-10min | 3219 | 3 | 6435 | 9657 | 3.711 | 3.7114 |
| restarting 3e-4 | polluted-48h-30min | 2908 | 0 | 5816 | 8724 | 3.519 | 3.5186 |
| restarting 3e-4 | clean-48h-5min | 12198 | 0 | 24396 | 36594 | 4.278 | 4.2783 |
