# Accuracy policy

For an atmospheric-model configuration, the target is **at least three digits
of overall relative RMS accuracy and at least two digits for every scored
species**, measured over the recorded states of each disclosed scenario. This corresponds to
0.1% overall RMS and at most 1% RMS for each species, with the supplied
concentration threshold. The two-digit tier is a useful lower-cost alternative,
not the default recommendation. A point barely clearing a threshold has less
margin for different meteorology and longer integrations.

This is a numerical-error budget for the prescribed chemistry problem, not a
claim that an atmospheric model predicts observations to 0.1%. Sandu et al.
(1997, REFERENCES.md sources 1–2) studied efficient atmospheric box-model
integration around the 1% region. That supports examining modest accuracy
rather than demanding roundoff agreement from every new method. Here the
recommended budget is stricter overall and also guards individual species.
Transport/splitting errors, mechanism uncertainties, emissions and uncertain
initial conditions remain separate. These experiments contain no transport and
do not establish a universal tolerance for every scientific application.

## What is measured

The authoritative metric is the supplied `../scenarios/score.py`. For reference
values at least `1e-14 mol m-3`, it forms relative errors over every cell,
snapshot and species; digits are minus the base-10 logarithm of their RMS.
The per-species requirement prevents abundant groups of easy species from
hiding a poor trace-species result. The final-time-only Tier 1 scenario uses
its final state; the other disclosed scenarios include every supplied snapshot.
Every output value, including those below the scoring threshold, must also be
finite and at least `-1e-20 mol m-3`. The runner checks finiteness/nonnegativity
at every external step, including steps that are not snapshots.

Local `rtol` is **not** the reported global accuracy. Errors can accumulate,
stiff transients can reduce effective order, clipping can change invariants,
and adaptive step sequences can cause nonmonotone work–precision curves.
Operating points are selected by measured global and species accuracy, not by
renaming their local tolerance as an accuracy guarantee.

The supplied reference solutions were independently checked for convergence
at substantially greater accuracy than the two ranking levels. Additional
48-hour cases use the supplied MICM six-stage solver at `(rtol,atol)` of
`(1e-9,1e-22)` and `(1e-10,1e-23)`; the latter defines the reference only after
checking their agreement. The recipes, hashes and scorer results are retained
in `tools/stress.py`, `results/stress_inputs.json` and
`scans/stress_references.csv`. These additional cases exercise summer, winter,
longer duration, new seeds, different starting hours and both clean and polluted
conditions. They are not the experimenter's undisclosed test set.

## Error weights and floors

New methods use `atol=1e-18 mol m-3` uniformly, without selecting special TS1
species. This is four orders below the scoring threshold. At that threshold,
the recommended `rtol=0.002` contributes `2e-17` to the error scale, so the
absolute floor does not dominate the smallest officially scored concentrations.
Below the threshold, species are still integrated with the same complete
mechanism. Their errors are bounded locally by the absolute tolerance but are
not guaranteed a particular relative accuracy by the official score.

The original TS1 default absolute tolerance is `1e-3 mol m-3`. It is retained
for compatibility, but is far too loose to control many atmospheric
concentrations. The explicit-tolerance MICM comparison therefore matters more
than speed against a failed default run. Tightening the original solver's
weights resolves its failures on the disclosed scenarios.

A controlled floor scan on 64 cells changed only `atol` from `1e-20` through
`1e-15`: worst overall accuracy for fourth-order extrapolation at `rtol=0.002`
fell from 3.588 to 3.238 digits. The chosen `1e-18` gave 3.453. Replacing the
maximum cell RMS norm by an all-cell RMS reduced that minimum to 3.177; the
component maximum increased it to 4.799 at substantially higher cost. We retain
the maximum cell RMS so a larger group cannot dilute a difficult cell by adding
easy cells. These are screening measurements; final points use all disclosed
cells and the official level checks.

## Positivity and conservation

Every new method rejects nonfinite trials. It includes the magnitude of a
negative candidate in the normalized error before accepting and then projects
accepted negative values to zero. The candidate's error contribution is the
maximum of estimated integration error and projection magnitude. This controls
projection locally but is not a proof of exact elemental conservation or a
uniform global bound. Intermediate-stage clipping is disabled in the declared
points. No reactions/species are removed or changed.

On a conservative linear transfer with a stiffness interval of `1e8`, fixed
fourth-order extrapolation and its automatic-order variant produce a final
product excess of approximately `1.65e-5` at `rtol=1e-4`. The accumulated
projection exactly accounts for that excess to the test's roundoff bound.
The other tested methods have roundoff-sized conservation errors on that case.
All eleven variants meet the test's 0.1% solution criterion. The separate smooth
linear-system test verifies near-roundoff conservation and expected convergence
orders without projection, and the nonlinear prototype comparisons also check
conservation. This distinguishes discretization/projection behavior from an
unaccounted implementation error.

For the recommended point on the 1,000-cell disclosed scan, total projected
concentration is zero in the Tier 1 scenario and at most `8.90e-14 mol m-3`
when summed over all species, all 1,000 cells and all accepted steps of a day.
This sum is a diagnostic, not an elemental mass. The scored reference agreement
provides the practical accuracy evidence. Users needing exact conserved
families should assess this numerical projection against their own invariant
budget; it is not silently described as conservative.

## Limits and recommendation

The recommended setting is fourth-order Euler extrapolation with fresh
Jacobians, proportional control, maximum cell RMS, warm step-size starts,
`rtol=0.002`, and `atol=1e-18`. The full disclosed sweep found minima of 3.448
digits overall and 2.550 for the worst species, giving margin beyond the chosen
criterion. OPERATING_POINTS.md and RESULTS.md report final repeated runs,
compiler checks, and the additional 48-hour results, including any operating
point that does not qualify at its intended level on a stress case.

The official score excludes relative error below its concentration threshold
and samples only the specified times. Neither a passing aggregate score nor a
higher formal order proves uniformly accurate instantaneous reaction rates,
all-time peaks, arbitrary new mechanisms, or all possible forcing histories.
The code remains mechanism-generic, and new mechanisms are supported by the
same kernels and tests, but scientific accuracy should be re-established when
the application changes. The experimenter's undisclosed scenarios are not
available here; expected accuracy is supported by disclosed and additional
scenario measurements, not a claim of having tested unavailable data.
