# Accuracy target

The intended atmospheric-model setting targets **at least three overall
significant digits (0.1% relative RMS error), and at least two digits for every
scored species**, with a measured margin on the disclosed cases. A faster
two-digit setting is also useful when the model's error budget permits it.
These are chemistry integration targets, not guarantees of atmospheric forecast
accuracy or of individual-value relative errors.

Sandu et al. (1997) motivate efficient low-to-moderate-accuracy sparse one-step
solvers for repeatedly restarted atmospheric chemistry. More recently,
Dreger et al. (2025) report gas-phase relative tolerances of `1e-2` and aerosol
tolerances of `1e-3` in a global model. Their error analysis also shows why a
mean alone can conceal individual-species deviations. We choose a tighter
three-digit recommendation to leave room in the error budget for splitting,
transport and other approximations; we measure the resulting global chemistry
error instead of equating local tolerance with global accuracy.
[Sources and exact uses](REFERENCES.md).

## Metric and what it proves

All scenario accuracy values come from the supplied `score.py`. For reference
concentrations at least `1e-14 mol m-3`, define
`r = abs(y - y_ref) / y_ref`. Overall digits are `-log10(sqrt(mean(r*r)))`
over all included species, cells and snapshots. Each species has its own RMS
and digits. Level L requires overall digits at least L and **every** species
at least L-1. The nonfinite/negativity gate checks all output concentrations,
including those below the accuracy threshold.

The metric assesses the specified chemistry with photolysis held constant in
each external interval. It does not measure the physical error of that forcing
approximation, transport coupling, or operator splitting.

The threshold corresponds to about 6,022 molecules cm-3. An atol of `1e-18`
mol m-3 is about 0.602 molecules cm-3; `1e-16` is about 60.2 molecules cm-3.
Absolute tolerances therefore materially affect low-abundance chemistry.
The default legacy TS1 atol (`1e-3 mol m-3`) was preserved for compatibility,
not selected as the new accuracy target. Its short-step benchmark accuracy
does not carry over to long chemistry intervals, as the disclosed failures show.

An RMS bound is not a per-value bound. The score covers the prescribed snapshots;
errors between those reporting times can be larger. Species below the threshold are excluded
from digit scoring, but still participate in every RHS and Jacobian evaluation.
No reactions or low-abundance species are removed. The supplied references use
a converged six-stage method at rtol `1e-10`, atol `1e-23`; reported reference
convergence exceeds 10 overall digits, comfortably beyond the targets here.

## Independent checks and positivity

All added C formulas agree with a dense Python prototype on nonlinear chemistry
with differing cells and changed rates, with identical operation counts.
Seven analytical stiff-chain tests check concentrations, linear conservation and
nonnegativity; step-halving tests recover the expected second, third or fourth
orders. The prototype implements its own dense LU and reuses factors explicitly.
These small-case checks establish implementation fidelity; scenario scoring
establishes application accuracy.

New methods reject nonfinite trials and concentrations more negative than
`positivity * atol`. Smaller negatives are projected to zero and the change is
charged to the error estimate before accepting. All accepted states are
nonnegative. Projection is not exactly conservative; its counts, sum and maximum
concentration correction are measured. It should not be mistaken for a proof
of chemical-element conservation or universal positivity of the underlying
Rosenbrock/SDIRK formula.

Independent stress cases vary seeds, latitude, season, starting time, pollution,
and duration (48 hours). Their references are made at two tight tolerances and
checked for convergence before scoring candidates. The undisclosed scenarios
cannot be inspected; finite accepted states are enforced by the algorithm,
while accuracy on unseen conditions remains an empirical generalization claim.
An interval that cannot meet the error tests fails explicitly rather than
silently accepting an invalid state.

## Measured margins of the declared points

The final clang measurements cover all cells of six disclosed cases and three
independent cases; GCC passes the same levels. Minima across those nine cases
are:

| Point | Overall digits | Worst-species digits | Limiting species and case |
| --- | ---: | ---: | --- |
| economy | 2.266 | 1.271 | ISOP, diurnal-30min |
| fast3 | 3.086 | 2.057 | ISOP, tropical-clean-48h |
| balanced | 3.433 | 2.490 | ISOP, tropical-clean-48h |
| conservative | 3.867 | 3.065 | BRO, diurnal-5min |

Balanced's lowest overall accuracy occurs in the polluted case. Its 0.433-digit
overall and 0.490-digit species margins correspond to roughly 2.7 and 3.1 times
smaller RMS errors than the target limits. Fast3's species margin is only
0.057 digits, which is why it is not the atmospheric-model recommendation.
Conservative does not meet four overall digits everywhere.

For balanced, the largest individual projection is `9.993e-20 mol m-3` and
the largest trajectory-total projection averaged per cell is `8.184e-20 mol m-3`.
Those increments are small relative to the scoring threshold, but that alone
does not bound indirect kinetic effects; the trajectory comparisons include
those effects. All nine final cases have minimum concentration zero and no
nonfinite trial rejections for every declared point. The full projection
measurements and per-scenario scores are in [RESULTS.md](RESULTS.md).

For larger
10,000-cell timing inputs, any accuracy scored only on the supplied 1,000-cell
prefix is labeled as such; it is not presented as a full 10,000-cell reference
comparison.
