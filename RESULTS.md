# Results

All runs used `-O3 -DNDEBUG` with the local Apple M5 Pro, clang 21.0.0 or GCC
15.3.0, and the repository's OpenMP setup. The generated TS1 source was
2,404,083 bytes. The clean `make bench` build took 62.03 seconds with clang and
62.32 seconds with GCC; this includes Python generation and compilation of the
generic library, generic benchmark, specialized benchmark, and scenario runner.

## Regression gates

| Compiler | Tier 1 cases | TS1 serial / parallel | Specialized checks |
| --- | --- | --- | --- |
| clang 21 | 662 PASS, 0 FAIL, 18 N/A of 680 | PASS / PASS | TS1, no-user-rate config, user-rate config PASS |
| GCC 15.3 | 662 PASS, 0 FAIL, 18 N/A of 680 | PASS / PASS | TS1, no-user-rate config, user-rate config PASS |

The 18 N/A cases are C++-specific API or type behavior with no C equivalent,
as listed by `tools/test_cases.py`. The generated comparisons use three cells
and check every state value against the generic library at `1e-10` relative
tolerance with a `1e-30` absolute floor. TS1 has 209 integrated species and 547
reactions; the two configuration mechanisms have 7 species and 7 or 10
reactions.

## TS1 benchmark

The ranking-shaped comparison used clang 21 for specialized C and Tier 1, and
clang++ for MICM, all with `-O3 -DNDEBUG`. It used 10,000 cells from
`make_cells_csv.py 10000 2026`, 100 timed intervals of 30 seconds, one untimed
warm-up interval, and the compiler's default 18-thread OpenMP team. Each solver
reported 16,200 internal steps, 16,200 accepted steps, and no rejected steps.
One-thread runs reported 900 internal and accepted steps. The single-run timings
were:

| Threads | Specialized C | Tier 1 generic C | MICM | Specialized speedup vs Tier 1 |
| ---: | ---: | ---: | ---: | ---: |
| 18 | 9,975.272 ms | 10,721.189 ms | 12,673.117 ms | 7.0% |
| 1 | 47,993.195 ms | 58,589.504 ms | 137,279.541 ms | 18.1% |

`compare_final.py` passed all 2,110,000 cell/species values against MICM at
both thread counts. The largest relative difference was `2.15e-13` for the
specialized result (BCARY, cell 4179) and `1.77e-13` for Tier 1 (BCARY, cell
6848). These are local single-run measurements; the experiment's hidden seed
and five-run median were not measured here.

## Day-night scenarios

Each scenario ran with 32-cell groups, `rtol=1e-6`, and `atol=1e-15`. Both
compilers passed `score.py --level 10` against MICM, with zero NaN or infinite
values and zero values below `-1e-20`. The reference column is the agreement
with the high-accuracy reference snapshots. Scenario elapsed times are local
single runs and include photolysis updates.

| Scenario | MICM steps | clang vs MICM: overall (worst species) | GCC vs MICM: overall (worst species) | Reference agreement | clang / GCC elapsed |
| --- | ---: | --- | --- | ---: | ---: |
| diurnal-30min | 48,619 | 14.26 (ISOP 13.69) | 14.26 (ISOP 13.67) | 4.75 digits | 631.525 / 737.229 ms |
| diurnal-10min | 90,980 | 14.18 (ISOP 13.64) | 14.18 (ISOP 13.62) | 5.05 digits | 1,200.017 / 1,427.346 ms |
| diurnal-5min | 145,891 | 14.15 (ISOP 13.59) | 14.15 (ISOP 13.60) | 5.21 digits | 1,954.433 / 2,360.184 ms |
| polluted-10min | 106,393 | 14.01 (ISOP 13.42) | 13.99 (ISOPNO3 13.31) | 5.40 digits | 1,388.764 / 1,647.069 ms |
| clean-10min | 84,273 | 14.19 (ISOP 13.60) | 14.20 (ISOP 13.62) | 4.91 digits | 1,114.975 / 1,346.685 ms |

The reference scores show the distance from MICM's gate-tolerance three-stage
solution to the supplied high-accuracy reference. The specialized runs retain
at least 13.31 digits for every species against MICM, above the gate requirement
of 9.
