# Accuracy: what is good enough for an atmospheric model

This file argues a target, then shows (once the scans are in) what each
operating point delivers against it on the disclosed scenarios. The
measurements use the experimenter's yardstick, `scenarios/score.py`:
significant digits of accuracy (SDA) against MICM's six-stage Rosenbrock at
relative 1e-10 and absolute 1e-23, overall and per species, over species at
or above 1e-14 mol m-3, with a gate of no NaN, no infinity, and no value
below -1e-20 mol m-3.

## The target, from the literature (REFERENCES.md)

1. **One percent is the level the field designed its solvers for.** The
   benchmarks that established Rosenbrock methods for atmospheric chemistry
   focused on "modest accuracies, say relative accuracies near 1 %", with
   higher accuracy "redundant for the actual practice of air pollution
   modelling" [S97a]; SDA 2 means 1 % in their measure, which is the
   experimenter's. Verwer and coauthors aimed at "two figures at most"
   [VBLS96]. Zhang et al. ask for 0.1 % so that the solver's share of a 1 %
   overall numerical error is small [Z11].
2. **Production models run between one and three digits.** GEOS-Chem's
   KPP Rosenbrock solver has been run at tolerances giving SDA of about 1
   [H07] and at relative tolerances from 5e-4 to 5e-2 over its versions
   [GC], where 1e-3 gave about 2.6 digits on ozone [E09]; EMAC runs its
   gas phase at relative 1e-2 [D25]. CAM-chem, the model that runs TS1,
   integrates it with a first-order backward Euler at a fixed 30-minute step
   with no local error estimate at all, converged to a relative change of
   1e-3 (1e-4 for some species) [CAM6]; replacing that solver by a
   second-order Rosenbrock method at the same step moved summer surface
   ozone over the eastern United States by up to 5 ppb [Su17], which puts
   its own error well above the one-percent level.
3. **Other errors are larger.** Rate constants carry 5 to 30 % uncertainty
   at 298 K and more in the cold upper troposphere, with 10 % on photolysis
   rates, giving 10 to 16 % uncertainty in ozone and OH [NE17]; operator
   splitting with 10 to 60 minute operators changes NOx and secondary
   inorganic aerosol by 5 to 35 % and ozone by up to 2 % [PMK16], with
   strategies differing by up to 10 % [SZY16]. Mechanism reductions
   accepted into GEOS-Chem change medians by 0.2 to 2 % and NO2 by up to
   10 % in places [S20, L23].
4. **Per-species accuracy matters more than the literature's norms admit.**
   The benchmarks take the mean or the maximum over species of a time-RMS
   error with species below 1e6 molecules cm-3 excluded [S97b, E09]; this
   experiment scores every species above 6e3 molecules cm-3 (1e-14 mol m-3)
   and asks for L - 1 digits in each, so radicals and trace species that the
   benchmarks dropped count here. That is the right test for a model whose
   products include OH, HO2, NO3 and the halogen radicals: a solver can
   reach three digits on average by getting ozone and CO right while
   mis-stepping a radical by 20 %.

**The target adopted:**
- **Two digits overall and one per species (level 2)** is the floor for a
  production run: the level the solvers were designed for, below every
  other error source by an order of magnitude or more, and already more
  than CAM-chem's current solver delivers.
- **Three digits overall and two per species (level 3)** is the recommended
  operating point for an atmospheric model: it keeps the solver's error one
  to two orders below the 1 % overall numerical target [Z11], leaves margin
  for the undisclosed conditions a global model meets (other latitudes,
  seasons, pollution levels, step lengths), and costs, as RESULTS.md shows,
  only a fraction more than level 2 with a third-order method. One digit per
  species at level 2 means a 10 % error in some species in some cells is
  acceptable; that is the regime of the mechanism reductions GEOS-Chem has
  accepted, not of a reference run.
- **Beyond three digits** is for reference solutions and solver tests, not
  for a model: a fourth digit costs a higher-order method or a much smaller
  step for no benefit that survives the rate-constant uncertainty.

## Margins

A declared point must reach its level on scenarios it has never seen. The
disclosed set has six scenarios; the undisclosed ones vary seeds, latitudes,
seasons, step lengths, and pollution. The scans measure, for each point, the
lowest overall and per-species digits over the six scenarios and the worst
species; a point is declared for a level only if its lowest per-species
score clears L - 1 by a margin, since the per-species threshold is the one a
new scenario is most likely to test (isoprene and the nitrogen reservoirs
sit near it in these runs). The margin asked is 0.2 digit (a factor 1.6 in
error), the spread seen between the disclosed scenarios themselves.

## What the measurements show

*(Filled from the scans: a table of the operating points with their lowest
overall and per-species digits on the six scenarios, the worst species, the
gate, and the margin to each level; then the day-night structure of the
error, where the error is made, and how the error estimate relates to the
true error, which sets the safety factor.)*
