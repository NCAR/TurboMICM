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
