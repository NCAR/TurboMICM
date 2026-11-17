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
