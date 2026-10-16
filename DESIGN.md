# Design

The implementation separates mechanism construction, kinetic evaluations, sparse
linear algebra, and integration. Model finalization computes the Jacobian pattern,
symbolic LU fill, and execution plans once. States own concentrations, conditions,
custom parameters, tolerances, and all integration workspaces. A finalized model
can be shared across threads; a state cannot be solved concurrently.

Implementation and tuning details will be recorded alongside validation as the
port is completed.
