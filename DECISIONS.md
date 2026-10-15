# Decisions

- C11 is the runtime language. Python is confined to generation and testing in
  the local prefix environment. The C library owns copies of mechanism metadata;
  states borrow an immutable mechanism and own their numerical storage.
- Stable rate-type sorting follows MICM's type order. Species initially retain
  configuration order, excluding parameterized species. A parameterized reactant
  multiplies its rate constant once per occurrence.
- States store each variable contiguously across cells. This supports portable
  OpenMP SIMD and will allow independent groups with separate adaptive steps.
- The configuration reader accepts the tutorial CAMP YAML and JSON structure.
  Default rate-coefficient units are molecule cm-3, converted to mol m-3 using
  1e-6 times Avogadro's constant. Custom rates retain caller units, matching the
  by-hand tutorials. Surface diffusion and molecular weight use their stated SI
  units. Branched reactions expand into alkoxy and nitrate reactions.
- Format extensions: a mechanism or reaction may declare `rate units: mol m-3`;
  species may declare `parameterization: air density`; reactants may be an
  ordered list of names and products an ordered list of `{species, yield}`
  entries, retaining repeated terms and arithmetic order; a mechanism may carry
  `initial state`, a mapping of temperature, pressure, concentrations and custom
  rate labels to values. These extensions represent TS1's existing data without
  converting coefficients twice or losing order. Additional analytic rate types
  use the parameter names of the corresponding C structures.
- `tools/translate_ts1.py` extracts data from the supplied C++ declaration forms,
  asserts the full species/reaction counts, and writes the retained JSON input.
  It neither compiles nor executes source code. The generic configuration reader
  then generates C builders through the same public API used by hand-built tests.
- A missing test implementation is a FAIL, never an N/A. N/A will be used only
  for individually reviewed cases with no meaningful C behavior.
- The supplied task prompt is an input, excluded from the submitted source tree.
