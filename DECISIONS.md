# Decisions

The configuration reader will resolve `camp-files` relative to the referring
file and accept both JSON and YAML. Generated C is checked from the same
mechanism representation used by the runtime library, so a benchmark does not
need a separate, hand-maintained chemistry implementation.

TS1 needs several extensions to the example format. A species can have
`parameterization: air density` (the third body is not a state variable),
`initial concentration [mol m-3]`, and a diffusion coefficient. A
`USER_DEFINED` reaction names a rate parameter and a `scaling factor`.
The top-level `initial conditions` object supplies temperature, pressure and
named rate-parameter defaults. Surface reactions retain a `gas-phase reactant`
alongside the ordinary reactant and product maps. Stoichiometric multiplicities
use the examples' `qty` and `yield` fields. These extensions encode data, not
special cases in the C solver.

A parameterized reactant can alternatively use a `parameterization` map with
`constant`, `temperature`, `pressure`, and `air density` coefficients. Its
concentration is the corresponding linear combination of the cell conditions;
multiple parameterized reactants multiply their concentrations into the rate.

The example configurations use molecule cm-3 rate units, as shown by their
hand-built equivalents. The generator converts rate coefficients to mol m-3
using reaction order. TS1 is already in molar units, so its configuration sets
`rate units: mol m-3`. Named user parameters are supplied in the library's
molar units, including emission parameters.

For constraints supplied in configuration, a top-level `constraints` list
supports `LINEAR` entries with a species-to-coefficient `terms` map and a
`constant`, and `EQUILIBRIUM` entries with reactant and product stoichiometry,
an `equilibrium constant`, and optional Van't Hoff parameters. Both name the
`algebraic species` whose row is replaced. This extension permits the same
ordering and C generation path to handle constrained mechanisms.

Generated mechanisms may use a symmetric minimum-fill ordering computed from
their reaction graph. The benchmark enables it; an independent fixture checks
ordered and original forms by species name. Hand-built mechanisms retain the
caller-specified order. The ordering changes storage and LU fill, not the
chemical equations.

Absent an explicit air density, rate updates calculate the ideal-gas density
from the current temperature and pressure. A C caller can override the density
per cell. The benchmark rejects nonpositive time steps and uses group size 64
by default; both group size and scheduling can be changed at run time.

Randomized linear-algebra cases use a fixed C random sequence to generate the
same matrix sizes and distributions as the reference tests. This avoids relying
on implementation-specific C++ random-engine sequences while checking the same
factorization and solve residuals at the reference tolerances.

The C API keeps relative tolerances in the selected solver parameters and
per-species absolute tolerances in the state. C state copies use an explicit
deep clone, while transferring ownership of a state uses an ordinary pointer
assignment. Named setters return errors instead of throwing exceptions.
The version accessor reports the compatibility version in the reference
source's version header, rather than a separate version for this C port.

External models declare their possible Jacobian entries explicitly and may
provide residual and gradient callbacks for algebraic constraints. The C API
uses separate model and constraint descriptors in place of C++ templates.
Optional model species and parameters are appended by a composition step;
duplicate names are rejected. Model callbacks and their contexts are borrowed,
so callers control their lifetime and can use independent contexts for different
symbolic plans.

Species metadata uses explicit tagged properties and deep copies in C. Phase
and system objects are borrowed views; an export operation creates the state
species array, preserving property values needed by the runtime. Strict
unused-species validation is an opt-in check of native reactions, because an
external callback can also use a species without a native reaction.
