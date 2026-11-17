# Design (work in progress)

`tools/mechanism.py` reads and validates the configuration, then
`tools/generate.py` emits immutable C tables. `chem.h` also exposes the table
format for mechanisms constructed directly in C. The tables list repeated
reactants explicitly, preserving the product rule even at zero concentration.
Rate constants, forcing, and Jacobian evaluation are separate from the solver
interface to allow replacement of the integration method and linear solver.

The current C routines use a dense row-major Jacobian. No threading or tuning
parameters have been introduced yet. Solver architecture and tuning details
will be documented when implemented.
