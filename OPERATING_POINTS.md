# Operating points

Five configurations of `bin/ts1_scenario` for the experimenter to evaluate,
each with the accuracy it reached on every disclosed scenario and the
accuracy expected on a scenario it has not seen (ACCURACY.md gives the target
and the margin rule, RESULTS.md the measurements, SCANS.md how the points
were found). Build with `make bench` (`CC=clang` by default); every point
runs any scenario at any thread count, with 32-cell groups and dynamic
scheduling (the defaults). `N` is the thread count and `SCENARIO` the
scenario's `scenario.txt`; the output is `ts1_snapshots.csv` in the working
directory.

| Point | Command | Declared level | Lowest overall / worst species over the six disclosed scenarios | Expected on an undisclosed scenario |
|---|---|---|---|---|
| L2 | `bin/ts1_scenario --threads N --method rodas3 --option rtol=5e-3 --option atol=1e-20 --option quantize=1 --option controller=micm --option safety=0.9 SCENARIO` | 2 | 2.42 / isoprene 1.29 (polluted-10min); 2.42 to 2.73 overall and 1.29 to 1.87 per species on the day-night scenarios | 2 digits overall and 1 per species, with 0.3 to spare |
| balanced | `bin/ts1_scenario --threads N --method rodas3 --option rtol=2e-3 --option atol=1e-20 --option quantize=1 --option controller=h211b --option safety=0.9 SCENARIO` | 2 | 3.17 / isoprene 2.17 (polluted-10min); 3.17 to 3.75 and 2.17 to 2.89 | 3 digits overall and 2 per species on most scenarios, with 0.17 to spare on the disclosed ones; declared for level 2 because that margin is below the 0.2 the rule asks |
| L3 | `bin/ts1_scenario --threads N --method extrap --option order=3 --option sequence=shifted --option rtol=1e-3 --option atol=1e-20 --option quantize=1 --option controller=micm --option safety=0.9 SCENARIO` | 3 | 3.61 / isoprene 2.49 (clean-10min); 3.61 to 4.12 and 2.49 to 3.00 | 3 digits overall and 2 per species, with 0.5 to spare |
| recommended | `bin/ts1_scenario --threads N --method extrap --option order=4 --option sequence=shifted --option rtol=5e-4 --option atol=1e-22 --option quantize=1 --option controller=micm --option safety=0.9 SCENARIO` | 3 | 3.62 / OClO 2.69 (diurnal-10min); 3.62 to 4.52 and 2.69 to 3.40 | 3 digits overall and 2 per species, with 0.6 to spare; 2.5 per species on most scenarios |
| accurate | `bin/ts1_scenario --threads N --method extrap --option order=4 --option sequence=shifted --option rtol=2e-4 --option atol=1e-20 --option quantize=1 --option controller=micm --option safety=0.9 SCENARIO` | 3 | 3.87 / BrO 3.00 (diurnal-5min); 3.87 to 4.81 and 3.00 to 3.66 | 3 digits overall and 2 per species with 0.9 to spare; about 4 overall and 3 per species on most scenarios |

Times (sums over the six disclosed scenarios of the median of three runs at
1,000 cells and 10 threads, S13): L2 1.16 s, balanced 2.21 s, L3 2.65 s,
recommended 3.04 s, accurate 3.45 s. MICM's three-stage Rosenbrock at the
tolerances that reach each level on the same scenarios, and the speedups,
are in RESULTS.md.

All five share the carried step quantized to the interval (`quantize=1`),
32-cell groups, the group RMS error norm, and photolysis rate constants
recomputed only when the custom parameters change. The option names are
explained in DESIGN.md (Tier 3, "Options"); `controller=micm`
is MICM's elementary step-size controller and `controller=h211b` the
digital-filter controller H211b with b = 4.

## The recommended point

**recommended**: extrapolated linearly implicit Euler of order 4 with the
shifted substep sequence (2, 3, 4, 5), relative tolerance 5e-4, absolute
tolerance 1e-22 mol m-3. It is the level-3 point with the largest margin for
the worst species at the least cost: 2.69 digits for OClO on diurnal-10min
and 2.85 or better everywhere else, 3.62 digits overall at the lowest, for
15 % more time than L3. By ACCURACY.md's argument, three digits overall and
two for every species keep the solver's error one to two orders of magnitude
below the other errors of an atmospheric model; this point keeps that with
a margin that the spread among the disclosed scenarios (0.1 to 0.2 digit)
does not consume, and it reaches four digits on the polluted and 30-minute
scenarios. Its cost is 7 % of MICM's three-stage method at the tolerances
that give MICM the same level (RESULTS.md).

For a model that needs only two digits, **L2** is the choice at a third of
the recommended point's time; **balanced** sits between them and will
usually deliver three digits. **accurate** is for reference runs: nearly
four digits everywhere, three per species, at 1.13 times the recommended
point's time.
