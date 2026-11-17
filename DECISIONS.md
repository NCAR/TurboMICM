# Decisions

The initial generic kernel represents integer reaction orders by repeated species
indices, and product stoichiometry by floating-point yields. This keeps the
Jacobian well-defined at zero concentration without dividing by a concentration.
The generator expands a branched reaction into alkoxy and nitrate reactions.
CAMP rate coefficients are converted from molecule/cm3 to mol/m3 using the
reaction order; user-defined parameters are already in the solver units. Surface
parameter labels use the gas-phase species name, as in the tutorial built by hand.
Unsupported types are rejected explicitly.

TS1 configuration extensions: `third body` marks a species whose concentration
is air density rather than a state variable; `rate coefficient units: SI` marks
coefficients already expressed in mole-based units; `USER_DEFINED` accepts a
`scaling factor`; and `surface parameter label` distinguishes the parameter
label from the gas species. Species may specify `initial concentration [mol m-3]`;
the manifest may provide custom parameter defaults, temperature, and pressure.
These retain the values in the generated C++ source without converting them twice.

Additional configuration extensions use `REVERSIBLE` (`A`, `C`, `k_r`),
`TAYLOR_SERIES` (Arrhenius parameters plus a `coefficients` array of length
1–16), and `TERNARY_CHEMICAL_ACTIVATION_JPL19` (k0, kinf, and kint A/B/C/D
fields plus Fc and N). Their coefficient unit conversion follows reaction
order, with the low-pressure JPL19 term having one additional concentration
order. `rate coefficient units: SI` disables conversion.
