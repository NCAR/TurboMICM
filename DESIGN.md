# Design

The implementation is organized into four layers:

1. `tools/generate.py` reads YAML or JSON using PyYAML and normalizes the
   configuration.
2. The generator writes a C descriptor containing species, reaction
   stoichiometry, rate parameters, and default state values.
3. `src/chem.c` builds solver contexts from descriptors and implements rate
   evaluation, Jacobian assembly, the linear solve, and Rosenbrock integration.
4. `src/ts1_bench.c` is a thin CSV and timing driver over the same library API.

Mechanism descriptors are immutable and shared. Per-cell state and mutable
solver workspaces are owned by a group context. This gives the solver and linear
algebra clear replacement boundaries without tying them to a mechanism.

Tuning parameters are `--threads` at run time and the cell group size at build or
run time; current defaults and verified values will be recorded here as the
implementation is completed. To add a Rosenbrock parameter set, add a parameter
descriptor and select it through the solver interface. To change error control,
replace the controller callback. To change the linear solver, implement the
factor/solve interface. To add a rate type, extend the descriptor enum, the
configuration normalization, and the generic C rate evaluator.
