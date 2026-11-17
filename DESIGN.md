# Design (in progress)

`tools/config.py` reads CAMP configuration manifests and their JSON or YAML data.
`tools/generate.py` currently emits static C tables for photolysis reactions,
with rate constants supplied by the caller. `src/chem.c` evaluates mass-action
forcing and its dense Jacobian. The derivative is formed by products excluding
one factor, so zero concentrations are handled without division. The interface
separates the mechanism from the kernel and leaves the choice of solver open.

This is a starting point, not the final architecture: no solver, threading,
tuning parameters, or linear solver have yet been implemented.
