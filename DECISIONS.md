# Decisions

The initial generic kernel represents integer reaction orders by repeated species
indices, and product stoichiometry by floating-point yields. This keeps the
Jacobian well-defined at zero concentration without dividing by a concentration.
The initial generator accepts photolysis reactions only and rejects other types
explicitly; this is an incomplete implementation, not a format limitation.
