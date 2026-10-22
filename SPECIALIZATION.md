# Mechanism specialization

The generator consumes CAMP JSON or YAML, including `camp-files`, using the Tier
1 reader. It is independent of species names and reaction counts. TS1 is an input
configuration, not a branch in the generator. Its 209 dynamic species and 547
reactions produce a 2,381-entry filled matrix and 6,986 elimination updates.

Generation fixes species ordering, sparse coordinates, fill, derivative jobs,
reaction arithmetic, stoichiometric yields, rate families and their fixed
coefficients, alpha shifts, elimination operations, and triangular solve rows.
Fused cell loops keep scalar reaction and elimination accumulation order. The
state's temperature, pressure, air density, concentrations, custom parameters,
step length and tolerances remain data. Generated formulas read them each time
rates are updated. No Python, parser, or third-party runtime is linked into C.

The integrator shares MICM's Rosenbrock parameter tables, embedded error estimate,
weighted error norm, and acceptance/rejection controller with the generic solver.
The generated path binds chemistry and linear algebra directly and knows the
species count. Optional fusion handles the three stages and their vector updates;
the generic method loop remains available. Reciprocal pivots avoid recalculating
the same diagonal divisions in the lower factor and all three back solves.

The alternatives and measured effects are recorded below after validation.
