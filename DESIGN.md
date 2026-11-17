# Design (in progress)

`tools/config.py` reads CAMP configuration manifests and their JSON or YAML data.
`tools/generate.py` emits static C tables for tutorial reaction types,
including rate parameters and labels for caller-supplied values. `src/rates.c`
calculates rate constants from conditions and custom parameters. `src/chem.c` evaluates mass-action
forcing and its dense Jacobian. The derivative is formed by products excluding
one factor, so zero concentrations are handled without division. The interface
separates the mechanism from the kernel and leaves the choice of solver open.

`src/rosenbrock.c` provides a single-cell autonomous ODE solver using callback
functions for forcing and Jacobian, configurable Rosenbrock coefficients, and
reference step-size defaults. It uses dense no-pivot Mozart elimination. The solver clamps negative ODE values after each public solve, as the reference
CPU solver does. The current implementation allocates scratch space per solve; multi-cell scheduling,
sparse matrices, other solver families, and constraints are not yet implemented.

This is a starting point, not the final architecture.
