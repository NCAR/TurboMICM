# Tier 3 ranking: 11 scenarios, 10 threads

Time: the sum over scenarios of each point's median elapsed_ms, in seconds.

## Level 2

| Rank | Entry | Fastest point reaching it everywhere | Time (s) | disclosed-10000 (s) | heldout-10000 (s) | Speedup over MICM |
|---:|---|---|---:|---:|---:|---:|
| 1 | fable-max | L2 | 24.55 | 11.07 | 13.48 | 18.93 |
| 2 | astra-xhigh | fast | 34.10 | 15.86 | 18.24 | 13.63 |
| 3 | sol-max | fast | 39.45 | 18.74 | 20.71 | 11.78 |
| 4 | sol-xhigh | L3 | 95.69 | 44.70 | 51.00 | 4.86 |
| 5 | micm | L2 | 464.74 | 218.37 | 246.37 | 1.00 |

## Level 3

| Rank | Entry | Fastest point reaching it everywhere | Time (s) | disclosed-10000 (s) | heldout-10000 (s) | Speedup over MICM |
|---:|---|---|---:|---:|---:|---:|
| 1 | fable-max | L3r | 49.08 | 21.96 | 27.12 | 11.19 |
| 2 | astra-xhigh | accurate | 73.46 | 33.99 | 39.47 | 7.47 |
| 3 | sol-xhigh | L3 | 95.69 | 44.70 | 51.00 | 5.74 |
| 4 | micm | L3 | 549.07 | 256.28 | 292.79 | 1.00 |
| | sol-max | none | |  |  |  |

## Every point

| Point | Level 2 | Level 3 | Lowest overall | Lowest worst species | Time (s) | disclosed-10000 (s) | heldout-10000 (s) | Runs | Scenarios with differing outputs |
|---|---|---|---:|---|---:|---:|---:|---:|---:|
| sol-xhigh/L3 | yes | yes | 3.343 | ISOP 2.311 | 95.69 | 44.70 | 51.00 | 33 | 0 |
| sol-max/fast | yes | no | 2.041 | ISOP 1.064 | 39.45 | 18.74 | 20.71 | 33 | 0 |
| astra-xhigh/fast | yes | no | 2.041 | ISOP 1.035 | 34.10 | 15.86 | 18.24 | 33 | 0 |
| astra-xhigh/accurate | yes | yes | 3.201 | OCLO 2.390 | 73.46 | 33.99 | 39.47 | 33 | 0 |
| fable-max/L2 | yes | no | 2.273 | ISOP 1.160 | 24.55 | 11.07 | 13.48 | 33 | 0 |
| fable-max/L3r | yes | yes | 3.183 | ISOP 2.326 | 49.08 | 21.96 | 27.12 | 33 | 0 |
| fable-max/L3 | yes | yes | 3.555 | ISOP 2.471 | 54.60 | 25.09 | 29.51 | 33 | 0 |
| fable-max/recommended | yes | yes | 3.619 | OCLO 2.692 | 62.14 | 28.88 | 33.27 | 33 | 0 |
| fable-max/accurate | yes | yes | 3.867 | BRO 2.997 | 71.91 | 32.89 | 39.02 | 33 | 0 |
| micm/L2 | yes | no | 2.487 | ISOP 1.242 | 464.74 | 218.37 | 246.37 | 33 | 0 |
| micm/L3 | yes | yes | 3.288 | ISOP 2.053 | 549.07 | 256.28 | 292.79 | 33 | 0 |

## Median elapsed time by scenario (s)

| Point | disclosed-10000/clean-10min | disclosed-10000/diurnal-10min | disclosed-10000/diurnal-30min | disclosed-10000/diurnal-5min | disclosed-10000/polluted-10min | disclosed-10000/tier1-30s | heldout-10000/clean-48h-5min | heldout-10000/polluted-48h-30min | heldout-10000/summer-15min | heldout-10000/tropics-20min | heldout-10000/winter60-10min |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| sol-xhigh/L3 | 7.533 | 8.270 | 6.444 | 10.152 | 9.868 | 2.430 | 17.434 | 11.368 | 7.335 | 6.962 | 7.898 |
| sol-max/fast | 2.995 | 3.192 | 2.587 | 4.681 | 3.727 | 1.557 | 8.078 | 4.115 | 2.457 | 2.727 | 3.331 |
| astra-xhigh/fast | 2.451 | 2.644 | 2.012 | 4.078 | 3.269 | 1.407 | 7.489 | 3.616 | 2.144 | 2.185 | 2.802 |
| astra-xhigh/accurate | 5.057 | 5.507 | 4.145 | 9.044 | 6.898 | 3.340 | 17.175 | 7.049 | 4.443 | 4.434 | 6.368 |
| fable-max/L2 | 1.685 | 1.946 | 1.588 | 2.685 | 2.306 | 0.861 | 5.032 | 2.964 | 1.713 | 1.750 | 2.020 |
| fable-max/L3r | 3.867 | 4.077 | 3.485 | 4.730 | 4.926 | 0.875 | 8.921 | 6.685 | 3.929 | 3.703 | 3.879 |
| fable-max/L3 | 3.685 | 4.332 | 3.668 | 6.129 | 5.190 | 2.089 | 11.379 | 6.055 | 3.591 | 3.716 | 4.769 |
| fable-max/recommended | 4.318 | 4.514 | 3.654 | 7.800 | 5.743 | 2.848 | 14.907 | 5.986 | 3.476 | 3.951 | 4.946 |
| fable-max/accurate | 4.931 | 5.329 | 4.727 | 8.281 | 6.688 | 2.936 | 15.742 | 7.685 | 4.644 | 4.843 | 6.110 |
| micm/L2 | 35.406 | 36.854 | 13.513 | 71.228 | 37.598 | 23.773 | 140.284 | 26.152 | 24.868 | 19.154 | 35.914 |
| micm/L3 | 42.699 | 44.276 | 19.467 | 78.608 | 46.995 | 24.233 | 155.633 | 36.927 | 32.661 | 25.546 | 42.025 |

