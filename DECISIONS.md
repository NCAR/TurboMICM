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

A species may use a `parameterization` object with `constant`, `temperature`,
`pressure`, and `air density` coefficients. Its concentration is the linear
combination of those conditions and it is excluded from state. Each occurrence
as a reactant multiplies the reaction rate; a parameterized product has no state
contribution. As in the reference, a reaction may have at most four such
reactants. The `CALLBACK` rate type binds a named C callback at runtime.

The C state interface uses explicit error codes, clone, move, and destruction
in place of exceptions and C++ special members. A bulk named update checks all
names and lengths before applying any column, avoiding partial updates after an
invalid input. Direct row-major arrays remain available for integrations and
callers that do not need name lookup.

Mechanism generation preserves a separate species catalog with used flags, so
queries include parameterized species even though those species do not occupy
state slots. The catalog retains configuration order when state variables are
reordered; this avoids making name queries depend on the solver permutation.

External algebraic models use caller-owned C callbacks and explicit dependency
indices. The Jacobian callback adds the negative derivative to an existing
row, matching the existing constraint interface and permitting the same sparse
solver path to handle built-in and external constraints.

The C matrix interface uses an explicit runtime group width, including a
row-major option, in place of compile-time C++ matrix policies. The caller
owns the resulting matrix and can access its padded storage; checked row and
column functions operate only on logical cells.
