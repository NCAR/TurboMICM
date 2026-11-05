# MICM's fastest layout on the eleven Tier 3 scenarios, 2026-11-05

The 4 November session (`../2026-11-04-tier3-paper-configuration/`) timed MICM's
loop in its default matrix layout only: standard ordering with in-place Mozart
LU, the `CpuSolverBuilderInPlace` defaults. On this machine MICM's fastest CPU
configuration, vector ordering of 8 cells with in-place Mozart LU in 32-cell
groups, ran 2.4 times faster than that layout on the Tier 2 benchmark
(`../2026-10-22-tier2-sol-max-verification/`: 4,475.8 against 10,625.4 ms at
18 threads). So that the short manuscript does not overstate the gain over
MICM, this session times the fastest layout on the eleven scenarios, at the
same tolerances, with the retaining solver again as the anchor to the first
session.

## Build

`micm_scenario_configs.cpp` is `tier3/scenarios/micm_scenario.cpp` (unchanged
since 2026-10-16) passed through `admin/micm_configs_source.py`, which replaces
only the solver type aliases with the `TS1_VECTOR`, `TS1_LU_SEPARATE`, and
`TS1_LU_DOOLITTLE` macros. Both binaries were built from it with Apple clang
21.0.0 (clang-2100.3.34.2) against MICM `e15898ce` (v3.7.0-259, the checkout in
`~/NSF-NCAR/ODE/benchmarks/micm`):

    SDKROOT=$(xcrun --sdk macosx --show-sdk-path) clang++ -std=c++20 -O3 -DNDEBUG \
        [-DTS1_VECTOR=8] -I $MICM/include -I $MICM/test/benchmark \
        micm_scenario_configs.cpp -o micm_scenario_{standard,vector8}

## Check (archived 64-cell Diurnal 10 min, rtol 3e-4, atol 1e-15, 32-cell groups)

- `micm_scenario_standard` reproduces the 4 November binary
  (`tier3-paper-work/bin/micm_scenario`) byte for byte: identical snapshots and
  statistics (3,422 steps, none rejected, as the Python restarting loop).
- `micm_scenario_vector8` takes the same 3,422 steps and scores the same
  (4.03 digits overall, worst species 2.80, gate and level 3 pass); its largest
  relative difference from the standard layout, over values of at least
  1e-14, is 3.1e-14.

## Session

`time_fastest.sh` is `time_paper.sh` with the programs file as a fourth
argument; `programs_fastest.txt` lists `micm8/C2` and `micm8/C3` (the vector8
binary at rtol 1e-2 and 3e-4, atol 1e-15, 32-cell groups, the 4 November MICM
settings) and `retaining/C2` and `retaining/C3` (the 4 November command lines;
the binary was rebuilt at 00:55 on 5 November from `7a23d26`, which changes
only the scenario driver's option parsing for backward Euler). Three rotating
rounds at 10 threads, 5 s pauses, the 10,000-cell scenarios of the first
session, every run gated and the first 64 cells scored against the archived
references. `compare.py` writes `compare.md` from this session's `out/` and the
first session's.

## Results (`compare.md`; session 02:41:41 to 03:19:38, `time.log`)

All 132 runs pass the gate. On every scenario at both criteria the vector8
layout takes the same internal steps and rejections and reaches the same
digits as the 4 November default layout. Sums over the eleven scenarios of the
per-scenario medians:

| Program | Criterion 2 (s) | Criterion 3 (s) |
|---|---:|---:|
| MICM, vector8, 32-cell groups | 139.18 | 165.12 |
| retaining solver (anchor; 4 November: 21.12, 30.60) | 21.50 | 31.30 |

| Ratio | Criterion 2 | Criterion 3 |
|---|---:|---:|
| MICM vector8 / retaining solver, same session | 6.47 (2.75 to 12.34) | 5.27 (2.27 to 11.59) |
| MICM default layout / MICM vector8, anchored | 2.44 | 2.49 |
| MICM vector8 / library restarting loop, anchored | 1.84 | 1.78 |

The short manuscript's Table 5 now gives MICM's loop in this layout, 6.5 and
5.3 times the retaining solver's time, in place of the default layout's 15.8
and 13.1; at criterion 3, 5.3 = 3.0 (step control) x 1.8 (kernels and data
layout).
