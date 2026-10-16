# Decisions

- The public library is called `chem`. Its implementation is ISO C11, with doubles
  throughout. Python is restricted to build-time configuration and test tooling.
- Numerical formulas, operation ordering where significant, solver parameters,
  error control, and LU algorithms follow the supplied reference source.
- C++ ownership and template tests will exercise equivalent C allocation, copying,
  destruction, layout, and error-return behavior. Only cases with no meaningful C
  equivalent will be marked N/A, individually and with reasons.
- The supplied experiment prompt remains an untracked local input. It is excluded
  because it contains experiment identity information.
- Reference paths default to sibling directories, overridable in Make variables.
  Builds and runs of reference executables will write only inside this repository.

## Configuration conventions and extensions

- Legacy configuration rate coefficients use molecule/cm3 concentration units.
  The reader converts them to mol/m3 according to reaction order. TS1's supplied
  C++ constants are already in mol/m3; an explicit reaction field
  `concentration units: mol m-3` suppresses conversion. The other accepted value
  is `molecule cm-3` (the legacy default).
- Reactants and products can use the original mappings, including `qty` and
  `yield`. Ordered arrays are also accepted: reactants are species-name strings,
  products are `{species: NAME, yield: VALUE}` objects. This extension preserves
  the exact order and duplicated terms of the C++ TS1 declarations.
- Species may specify `parameterization: air density` for a third body, plus
  `initial concentration [mol m-3]`. The document can supply `initial conditions`
  (temperature and pressure) and `initial rate parameters` (label/value mapping).
  These retain all of TS1 InitState, including its original six photolysis rates
  of 1/s; the benchmark input overrides those values as required.
- `USER_DEFINED` names the general custom-rate type. `PHOTOLYSIS`, `EMISSION`,
  and `FIRST_ORDER_LOSS` retain the legacy convenience forms. The original
  `MUSICA name` key is a data-format label only; no such package is imported.
- Tutorial surface labels differ between file examples (`C surface`) and the
  hand-written programs (`C`). Equivalence tests must map the same physical
  radius/number parameters by their respective labels.
- Species ordering uses the reference's diagonal Markowitz selection and stable
  tie-breaking. Automatic generated kernels use the same symbolic ordering as
  the runtime builder. `preserve_order` supports mechanisms that explicitly need
  their declared variable order.
- `make env` disables conda environment registration and keeps its package cache
  local. This prevents environment creation from updating the global environment
  registry on subsequent runs.
