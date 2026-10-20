# Design

## Current architecture

Configuration parsing runs in Python at build or development time. The reader
resolves the configuration's referenced files and normalizes species and
reaction records. The TS1 importer is separate because the permitted benchmark
input is MICM's generated mechanism builder. Its output records species order,
third-body flags, molecular weights, diffusion coefficients, stoichiometry,
rate types, and rate parameters.

`tools/generate_c.py` converts the normalized JSON into immutable C tables. The
C API is mechanism-generic: mechanisms are arrays of species and reactions;
each reaction refers to its reactants and products by species index. The forcing
routine applies mass action and accumulates stoichiometric production and loss.
Rate calculation is selected by a rate-kind enum, so rate implementations can be
replaced independently of mechanism storage.

The generated C tables now build into the CPU executable. Solver state grouping,
the integrator, and benchmark execution are not implemented yet. There are
therefore no tuned group-size, vector-width, or thread-count defaults at this
stage.

## Planned extension points

Solver methods and parameter sets should be data-driven structures consumed by
an integrator interface. The controller should have its own interface for error
norms and step-size decisions. Linear algebra should be selected through a
solver interface, with matrix structure stored independently. New rate types
should add a parameter record and rate evaluator and be reachable from the same
configuration-to-C generator used by tests.
