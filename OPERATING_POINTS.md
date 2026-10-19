# Operating points

Up to five configurations of `bin/ts1_scenario` for the experimenter to
evaluate, each with the accuracy it is expected to reach on every scenario
(ACCURACY.md gives the target and the margins, RESULTS.md the measurements
behind the expectations). The experimenter builds with `CC=clang` and times
10,000-cell scenarios at all cores; every point below runs any scenario and
any thread count.

*(Filled from the scans: a table with, per point, its name, the full
command line, the level it is declared for, the lowest overall and
per-species digits it reached over the six disclosed scenarios, the margin
to its level, and its time relative to MICM's three-stage Rosenbrock at the
tolerance that reaches the same level.)*

| Point | Command | Declared level | Lowest overall / worst species (disclosed) | Expected on an undisclosed scenario |
|---|---|---|---|---|
| | | | | |

## The recommended point

*(Which point, and why: the level-3 point whose per-species margin is
largest at the least cost, per ACCURACY.md's argument that three digits
overall and two per species keep the solver's error one to two orders
below the other errors of an atmospheric model.)*
