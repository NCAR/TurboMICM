# Decisions

The initial generic kernel represents integer reaction orders by repeated species
indices, and product stoichiometry by floating-point yields. This keeps the
Jacobian well-defined at zero concentration without dividing by a concentration.
The generator expands a branched reaction into alkoxy and nitrate reactions.
CAMP rate coefficients are converted from molecule/cm3 to mol/m3 using the
reaction order; user-defined parameters are already in the solver units. Surface
parameter labels use the gas-phase species name, as in the tutorial built by hand.
Unsupported types are rejected explicitly.
