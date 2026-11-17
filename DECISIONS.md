# Decisions

- Configuration files use the tutorial `camp-files`/`camp-data` structure and
  preserve species and reaction ordering. JSON and YAML use the same reader.
  The C generator emits immutable tables, so the runtime does not need Python.
- TS1 adds `third body` on a species, `USER_DEFINED` with a `scaling factor`,
  and `initial conditions` on the root configuration. A third body contributes
  air density to reaction rates but has no state variable. The translated
  surface reactions retain both the phase reactant and the explicit reactant
  and product lists. These extensions express details in the C++ source which
  are absent from the examples.
- Repeated reactants are represented by a positive integer `qty`; products use
  `yield`. The importer checks the source counts and fails on unrecognized
  builder fields instead of silently dropping rate data.
