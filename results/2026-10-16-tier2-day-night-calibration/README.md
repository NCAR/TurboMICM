# Calibrating Tier 2's day-night gate, 2026-10-16

How closely can two correct implementations of MICM's three-stage Rosenbrock
agree over a day and a night? MICM's scenario tool, `tier3/scenarios/micm_scenario.cpp`,
was built twice with `-O3 -DNDEBUG`, once with Apple clang 21 and once with
gcc-15, whose different contraction of multiply-adds perturbs the arithmetic at
roundoff. Both ran the five disclosed day-night scenarios (1,000 cells each) at
the settings of Tier 2's day-night tests:

    micm_scenario --threads 10 --group-size 32 --atol 1e-15 SCENARIO

(relative tolerance MICM's default, 1e-6; absolute 1e-15 mol m-3 for every
species, instead of MICM's default 1e-3, with which MICM fails three of these
five scenarios). Agreement is `score.py` with the clang result as the
reference; accuracy is `score.py` against the converged reference solution.

| Scenario | Internal steps (both builds) | Agreement, overall digits | Worst species | Accuracy vs reference | MICM (clang), 10 threads |
|---|---:|---:|---|---:|---:|
| diurnal-30min | 48,619 | 14.37 | ISOP 13.81 | 4.75 | 4.1 s |
| diurnal-10min | 90,980 | 14.28 | ISOP 13.70 | 5.05 | 7.7 s |
| diurnal-5min | 145,891 | 14.24 | ISOP 13.69 | 5.21 | 12.8 s |
| polluted-10min | 106,393 | 14.07 | ISOP 13.47 | 5.40 | 9.7 s |
| clean-10min | 84,273 | 14.35 | ISOP 13.84 | 4.91 | 7.5 s |

- Both builds take exactly the same internal steps on every scenario, none
  rejected: roundoff does not change the step sequence over 24 hours.
- They agree to 14 digits overall and at least 13.47 for every species.
- So Tier 2's gate, agreement with MICM to 10 digits overall and 9 for every
  species with 32-cell groups, leaves more than three digits of margin for
  arithmetic reordered by specialization.

`calibration.log` has the raw output.
