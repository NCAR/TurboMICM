# Tier 3 ranking: 11 scenarios, 10 threads

Time: the sum over scenarios of each point's median elapsed_ms, in seconds.

## Level 2

| Rank | Entry | Fastest point reaching it everywhere | Time (s) | disclosed-10000 (s) | heldout-10000 (s) | Speedup over MICM |
|---:|---|---|---:|---:|---:|---:|
| 1 | astra-xhigh | fast | 34.72 | 16.15 | 18.58 | 13.40 |
| 2 | sol-max | fast | 39.80 | 18.79 | 21.01 | 11.69 |
| 3 | sol-xhigh | L2 | 54.52 | 25.79 | 28.74 | 8.54 |
| 4 | astra-max | economy | 62.13 | 29.02 | 33.12 | 7.49 |
| 5 | micm | L2 | 465.40 | 221.06 | 244.34 | 1.00 |

## Level 3

| Rank | Entry | Fastest point reaching it everywhere | Time (s) | disclosed-10000 (s) | heldout-10000 (s) | Speedup over MICM |
|---:|---|---|---:|---:|---:|---:|
| 1 | astra-xhigh | accurate | 74.47 | 34.42 | 40.05 | 7.29 |
| 2 | sol-xhigh | L3 | 93.67 | 43.58 | 50.09 | 5.79 |
| 3 | sol-max | accurate | 106.57 | 50.65 | 55.92 | 5.09 |
| 4 | astra-max | balanced | 109.29 | 48.52 | 60.77 | 4.97 |
| 5 | micm | L3 | 542.66 | 254.96 | 287.70 | 1.00 |

## Every point

| Point | Level 2 | Level 3 | Lowest overall | Lowest worst species | Time (s) | disclosed-10000 (s) | heldout-10000 (s) | Runs | Scenarios with differing outputs |
|---|---|---|---:|---|---:|---:|---:|---:|---:|
| sol-xhigh/L2 | yes | no | 2.456 | ISOP 1.445 | 54.52 | 25.79 | 28.74 | 55 | 0 |
| sol-xhigh/L2_margin | yes | no | 2.800 | ISOP 1.689 | 68.34 | 32.21 | 36.13 | 55 | 0 |
| sol-xhigh/L3 | yes | yes | 3.343 | ISOP 2.311 | 93.67 | 43.58 | 50.09 | 55 | 0 |
| sol-xhigh/recommended | yes | yes | 3.707 | ISOP 2.734 | 123.26 | 58.27 | 64.99 | 55 | 0 |
| sol-xhigh/L4 | yes | yes | 4.659 | ISOP 3.616 | 322.50 | 149.81 | 172.69 | 55 | 0 |
| sol-max/fast | yes | no | 2.041 | ISOP 1.064 | 39.80 | 18.79 | 21.01 | 55 | 0 |
| sol-max/balanced | yes | no | 2.290 | ISOP 1.386 | 44.75 | 21.12 | 23.63 | 55 | 0 |
| sol-max/accurate | yes | yes | 3.251 | ISOP 2.341 | 106.57 | 50.65 | 55.92 | 55 | 0 |
| sol-max/recommended | yes | yes | 3.571 | ISOP 2.554 | 137.36 | 65.36 | 72.00 | 55 | 0 |
| sol-max/precise | yes | yes | 4.251 | ISOP 3.069 | 218.38 | 103.30 | 115.08 | 55 | 0 |
| astra-xhigh/fast | yes | no | 2.041 | ISOP 1.035 | 34.72 | 16.15 | 18.58 | 55 | 0 |
| astra-xhigh/balanced | yes | no | 2.463 | PHENO 1.569 | 47.60 | 22.38 | 25.22 | 55 | 0 |
| astra-xhigh/accurate | yes | yes | 3.201 | OCLO 2.390 | 74.47 | 34.42 | 40.05 | 55 | 0 |
| astra-xhigh/recommended | yes | yes | 3.448 | BRO 2.550 | 86.32 | 40.20 | 46.12 | 55 | 0 |
| astra-xhigh/strict | yes | yes | 3.861 | ISOP 2.919 | 104.15 | 48.52 | 55.63 | 55 | 0 |
| astra-max/economy | yes | no | 2.267 | ISOP 1.271 | 62.13 | 29.02 | 33.12 | 55 | 0 |
| astra-max/fast3 | yes | no | 2.463 | PHENO 1.569 | 86.41 | 39.89 | 46.52 | 55 | 0 |
| astra-max/balanced | yes | yes | 3.274 | ISOP 2.263 | 109.29 | 48.52 | 60.77 | 55 | 0 |
| astra-max/conservative | yes | yes | 3.776 | ISOP 2.951 | 163.69 | 71.45 | 92.24 | 55 | 0 |
| micm/L2 | yes | no | 2.487 | ISOP 1.242 | 465.40 | 221.06 | 244.34 | 55 | 0 |
| micm/L3 | yes | yes | 3.288 | ISOP 2.053 | 542.66 | 254.96 | 287.70 | 55 | 0 |
| micm/L3b | yes | yes | 3.425 | ISOP 2.206 | 551.06 | 258.62 | 292.44 | 55 | 0 |

## Median elapsed time by scenario (s)

| Point | disclosed-10000/clean-10min | disclosed-10000/diurnal-10min | disclosed-10000/diurnal-30min | disclosed-10000/diurnal-5min | disclosed-10000/polluted-10min | disclosed-10000/tier1-30s | heldout-10000/clean-48h-5min | heldout-10000/polluted-48h-30min | heldout-10000/summer-15min | heldout-10000/tropics-20min | heldout-10000/winter60-10min |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| sol-xhigh/L2 | 4.052 | 4.482 | 3.416 | 6.382 | 5.271 | 2.186 | 11.694 | 5.389 | 3.288 | 3.628 | 4.737 |
| sol-xhigh/L2_margin | 5.401 | 5.778 | 4.448 | 7.447 | 6.871 | 2.261 | 13.206 | 7.502 | 4.677 | 4.760 | 5.985 |
| sol-xhigh/L3 | 7.313 | 8.056 | 6.385 | 10.056 | 9.422 | 2.346 | 17.285 | 11.100 | 7.114 | 6.841 | 7.752 |
| sol-xhigh/recommended | 10.029 | 10.617 | 8.200 | 12.873 | 13.769 | 2.778 | 21.528 | 14.967 | 9.691 | 9.010 | 9.797 |
| sol-xhigh/L4 | 28.554 | 28.635 | 15.208 | 41.428 | 32.559 | 3.427 | 78.599 | 31.474 | 25.056 | 18.763 | 18.796 |
| sol-max/fast | 3.031 | 3.242 | 2.589 | 4.659 | 3.758 | 1.511 | 8.173 | 4.168 | 2.483 | 2.783 | 3.405 |
| sol-max/balanced | 3.470 | 3.757 | 2.953 | 4.854 | 4.528 | 1.559 | 8.514 | 5.011 | 3.008 | 3.217 | 3.881 |
| sol-max/accurate | 8.781 | 9.456 | 7.209 | 10.653 | 12.597 | 1.950 | 17.084 | 13.653 | 8.524 | 8.000 | 8.664 |
| sol-max/recommended | 11.671 | 12.320 | 8.828 | 14.068 | 16.282 | 2.188 | 22.960 | 17.180 | 11.414 | 9.950 | 10.499 |
| sol-max/precise | 19.503 | 19.872 | 12.692 | 24.962 | 23.757 | 2.514 | 43.280 | 24.225 | 17.704 | 14.776 | 15.100 |
| astra-xhigh/fast | 2.540 | 2.726 | 2.037 | 4.137 | 3.310 | 1.397 | 7.522 | 3.730 | 2.195 | 2.242 | 2.888 |
| astra-xhigh/balanced | 3.742 | 3.986 | 3.162 | 4.977 | 5.016 | 1.494 | 8.377 | 5.894 | 3.370 | 3.430 | 4.146 |
| astra-xhigh/accurate | 5.141 | 5.599 | 4.164 | 9.248 | 6.986 | 3.280 | 17.371 | 7.229 | 4.489 | 4.482 | 6.476 |
| astra-xhigh/recommended | 5.859 | 6.856 | 5.109 | 9.944 | 8.984 | 3.443 | 18.639 | 9.276 | 5.244 | 5.430 | 7.533 |
| astra-xhigh/strict | 7.560 | 8.442 | 6.684 | 10.955 | 11.266 | 3.611 | 19.988 | 12.699 | 6.703 | 7.156 | 9.087 |
| astra-max/economy | 4.778 | 5.116 | 3.682 | 7.207 | 6.074 | 2.159 | 12.996 | 6.970 | 4.080 | 3.972 | 5.102 |
| astra-max/fast3 | 6.773 | 7.443 | 6.445 | 8.211 | 8.696 | 2.322 | 13.884 | 12.290 | 6.295 | 6.680 | 7.366 |
| astra-max/balanced | 8.209 | 8.789 | 7.930 | 11.293 | 9.992 | 2.307 | 20.397 | 15.364 | 8.053 | 8.125 | 8.829 |
| astra-max/conservative | 13.164 | 13.734 | 13.020 | 14.487 | 14.734 | 2.313 | 27.169 | 25.841 | 12.427 | 13.083 | 13.720 |
| micm/L2 | 37.076 | 37.979 | 13.553 | 71.958 | 37.820 | 22.676 | 135.771 | 26.924 | 25.373 | 19.360 | 36.914 |
| micm/L3 | 43.450 | 44.478 | 19.354 | 78.112 | 45.894 | 23.671 | 150.345 | 37.384 | 32.266 | 26.055 | 41.650 |
| micm/L3b | 43.931 | 44.483 | 20.458 | 79.746 | 45.872 | 24.129 | 152.811 | 38.785 | 33.099 | 26.484 | 41.257 |

