# MICM's fastest layout on the eleven Tier 3 scenarios (5 November 2026)

Per-scenario medians (s) over the rounds; digits on the first 64 cells; gate on all cells.

| Scenario | micm8/C2 | micm8/C3 | retaining/C2 | retaining/C3 | micm8/ret C2 | micm8/ret C3 |
|---|---:|---:|---:|---:|---:|---:|
| tier1-30s | 6.67 (5.63) | 6.68 (6.48) | 0.54 (5.82) | 0.58 (5.84) | 12.34 | 11.59 |
| diurnal-30min | 3.83 (2.52) | 5.37 (3.34) | 1.39 (2.70) | 2.26 (3.29) | 2.75 | 2.38 |
| diurnal-10min | 10.71 (3.46) | 13.03 (4.03) | 1.74 (3.36) | 2.52 (3.56) | 6.15 | 5.18 |
| diurnal-5min | 21.35 (3.40) | 24.35 (4.36) | 2.19 (3.44) | 2.75 (3.51) | 9.74 | 8.86 |
| polluted-10min | 11.33 (3.30) | 13.75 (3.87) | 2.17 (3.25) | 3.47 (3.42) | 5.23 | 3.96 |
| clean-10min | 10.87 (3.44) | 13.20 (3.98) | 1.69 (3.45) | 2.33 (3.57) | 6.45 | 5.66 |
| summer-15min | 7.35 (3.25) | 9.58 (3.74) | 1.59 (3.10) | 2.39 (3.31) | 4.63 | 4.00 |
| tropics-20min | 5.54 (2.71) | 7.43 (3.43) | 1.52 (2.90) | 2.36 (3.21) | 3.64 | 3.14 |
| winter60-10min | 10.61 (3.14) | 12.49 (3.71) | 1.75 (3.15) | 2.52 (3.40) | 6.07 | 4.96 |
| polluted-48h-30min | 7.91 (2.62) | 10.97 (3.52) | 2.75 (2.88) | 4.84 (3.42) | 2.88 | 2.27 |
| clean-48h-5min | 43.01 (3.69) | 48.26 (4.28) | 4.18 (3.50) | 5.28 (3.58) | 10.30 | 9.14 |

Rounds per program and scenario: [3]; every run passing the gate: True

| Sum over the eleven scenarios (s) | 5 Nov | 4 Nov | 5 Nov / 4 Nov |
|---|---:|---:|---:|
| retaining/C2 (anchor) | 21.50 | 21.12 | 1.018 |
| retaining/C3 (anchor) | 31.30 | 30.60 | 1.023 |
| micm8/C2 | 139.18 | | |
| micm8/C3 | 165.12 | | |
| micm/C2 | | 333.14 | |
| micm/C3 | | 402.24 | |
| loop/C2 | | 74.47 | |
| loop/C3 | | 90.85 | |

| Ratio | All | Per-scenario range |
|---|---:|---:|
| micm8/C2 / retaining/C2, same session | 6.47 | 2.75 to 12.34 |
| micm8/C3 / retaining/C3, same session | 5.27 | 2.27 to 11.59 |
| micm/C2 (default layout) / micm8/C2, anchored | 2.44 | 2.38 to 2.59 |
| micm/C3 (default layout) / micm8/C3, anchored | 2.49 | 2.37 to 2.59 |
| micm8/C2 / loop/C2 (library kernels), anchored | 1.84 | 1.78 to 1.93 |
| micm8/C3 / loop/C3 (library kernels), anchored | 1.78 | 1.73 to 1.84 |
