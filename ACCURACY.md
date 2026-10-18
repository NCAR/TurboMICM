# Accuracy target

We recommend a chemistry integration error budget of **0.1% RMS relative error
(three significant digits)** in the supplied metric, with **at least 1% RMS
accuracy for every scored species**. The lower-cost alternative targets 1%
overall and 10% per species. The latter is a throughput option, not our default
recommendation for an atmospheric model.

This is an engineering budget for numerical integration, not a claim that
atmospheric concentrations are known to 0.1%, or that every use of the mechanism
has that accuracy requirement. A host model should repeat the convergence study
for its scientific diagnostics, emissions, transport coupling and time horizon.
No species or reaction has been removed or lumped.

## Basis

Sandu et al. (1997), pp. 3465–3469, focus their atmospheric chemistry
work–precision comparisons on roughly 1% accuracy. They find Rosenbrock methods
competitive in the low to moderate accuracy regime, and stress that local
relative tolerance and measured global accuracy differ. This motivates testing
1% as a useful lower-cost target, then budgeting another decimal digit for the
recommended configuration. See [REFERENCES.md](REFERENCES.md) for sources and
exactly what was used. The extra digit is our design choice, not a universal
threshold prescribed by that paper.

The supplied default MICM three-stage configuration uses atol = 1e-3 mol/m3.
That value is enormous compared with many chemically active TS1 species and
the scoring threshold of 1e-14 mol/m3. A small nominal rtol does not rescue this
absolute tolerance. Disclosed reference logs already show default failures and
poor long-run accuracy. Our experiments therefore set both tolerances explicitly
and measure global error rather than treating rtol as the result.

## What is measured

All reported accuracy values come from the unmodified supplied
`../scenarios/score.py`, not a substitute norm. It pools squared relative errors
for every reference concentration >= 1e-14 mol/m3, across cells, species and
recorded times. Digits are -log10(RMS relative error). It also reports each
species separately. Passing level L requires overall >= L and every species
>= L-1. Nonfinite values or concentrations below -1e-20 mol/m3 fail the gate.
This pooled metric differs from the paper's worst-species aggregate; we report
both overall and the worst species and do not conflate the definitions.

The disclosed references were generated at rtol = 1e-10 and atol = 1e-23, with
six-stage Rosenbrock, and checked against a second tight calculation. The supplied
convergence figures are 10.18–10.79 digits overall and 8.86–9.77 digits for the
worst species. Their numerical uncertainty is comfortably smaller than our
2- and 3-digit targets.

The solver's local norm uses atol_j + rtol max(|y_j|, |y_new,j|). The default
experimental norm is the maximum of the species RMS in each cell, so increasing
the number of cells cannot dilute one difficult cell's error. An optional
component maximum is stricter; group RMS is retained for comparison. Negativity
of a candidate is also normalized by these weights and can reject a step. After
acceptance, residual negative values are projected to zero. This projection is
not exact mass conservation; the small-case tests check invariants where no
projection is needed, while the scenario scores measure its net concentration
error. No post-processing adjustment toward reference values is performed.

## Limits of the evidence

Snapshots include day, night and transitions, but they do not prove a maximum
error bound at every intermediate instant. Species below the reference threshold
are unscored for relative accuracy; they are still integrated and must pass the
finite/nonnegative gate. The worst-species margin protects against an abundant
species dominating pooled RMS, but is not a guarantee for every individual cell.
Disclosed scenarios alone cannot prove performance or accuracy on undisclosed
conditions. Operating-point margins, additional independent conditions, and
convergence sweeps are documented in RESULTS.md before a point is recommended.

## Measured operating-point margins

The final disclosed-case minima are:

| Point | Overall digits | Worst-species digits | Overall RMS error | Worst-species RMS error |
|---|---:|---:|---:|---:|
| L2 | 2.541 | 1.445 | 0.288% | 3.59% |
| L2_margin | 2.799 | 1.689 | 0.159% | 2.05% |
| L3 | 3.341 | 2.322 | 0.0456% | 0.476% |
| recommended | 3.705 | 2.734 | 0.0197% | 0.185% |
| L4 | 4.657 | 3.614 | 0.00220% | 0.0243% |

Thus the recommendation has about five times less overall and species RMS
error than its chosen 0.1% / 1% budget. The fastest L3 point also meets that
budget, with less margin from its looser absolute tolerance. The L2 points
serve the ranking's lower-cost accuracy level; their species budget is looser
than the practical worst-species 1% emphasis of Sandu et al., which is another
reason not to recommend them as the atmospheric-model default.

All five points also pass their declared levels on three independent 48-hour
runs. For the recommendation those extra cases have at least 3.923 digits
overall and 3.194 for every scored species. The references agree to roughly
ten digits overall under a separate tolerance refinement. This strengthens the
case for the chosen margin; it does not turn the sampled error metric into a
uniform bound on all cells, times or future atmospheric conditions.
