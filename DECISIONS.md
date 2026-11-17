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
- Tutorial rate constants are interpreted in molecule cm-3 units, as the
  by-hand examples' conversion factors specify. `rate units: mol m-3` marks
  constants already converted to SI; TS1 uses this extension because its C++
  builder already contains converted coefficients. Reaction order determines
  the conversion, including the distinct low- and high-pressure falloff units.
- A cell is currently the unit of solver state and error control. This makes
  independently varying conditions straightforward and allows parallel work
  without shared mutable solver state, at the cost of repeated scalar work.
- The optional `ChemConditions` air-density override supports mechanisms and
  tests that prescribe density independently of temperature and pressure; when
  it is absent, the ideal-gas density is calculated as in the benchmark.
- The JPL-19 falloff, reversible, and Taylor-series rate types use explicit
  configuration type names `TERNARY_JPL19`, `REVERSIBLE`, and `TAYLOR`.
  Taylor coefficients are a list of one to sixteen numbers.
