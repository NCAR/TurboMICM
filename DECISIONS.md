# Decisions

The configuration reader will resolve `camp-files` relative to the referring
file and accept both JSON and YAML. Generated C is checked from the same
mechanism representation used by the runtime library, so a benchmark does not
need a separate, hand-maintained chemistry implementation.

TS1 needs several extensions to the example format. A species can have
`parameterization: air density` (the third body is not a state variable),
`initial concentration [mol m-3]`, and a diffusion coefficient. A
`USER_DEFINED` reaction names a rate parameter and a `scaling factor`.
The top-level `initial conditions` object supplies temperature, pressure and
named rate-parameter defaults. Surface reactions retain a `gas-phase reactant`
alongside the ordinary reactant and product maps. Stoichiometric multiplicities
use the examples' `qty` and `yield` fields. These extensions encode data, not
special cases in the C solver.

The example configurations use molecule cm-3 rate units, as shown by their
hand-built equivalents. The generator converts rate coefficients to mol m-3
using reaction order. TS1 is already in molar units, so its configuration sets
`rate units: mol m-3`. Named user parameters are supplied in the library's
molar units, including emission parameters.

Generated mechanisms may use a symmetric minimum-fill ordering computed from
their reaction graph. The benchmark enables it; an independent fixture checks
ordered and original forms by species name. Hand-built mechanisms retain the
caller-specified order. The ordering changes storage and LU fill, not the
chemical equations.

Absent an explicit air density, rate updates calculate the ideal-gas density
from the current temperature and pressure. A C caller can override the density
per cell. The benchmark rejects nonpositive time steps and uses group size 64
by default; both group size and scheduling can be changed at run time.
