# Work remaining

The full task is not complete. The last verified checkpoint has 520 passing,
156 failing, and four N/A cases, with both clang and gcc-15. The full provided
TS1 fixture passes both compilers. A full-size clang TS1 run is bit-identical to
MICM for all 2,110,000 concentrations, and locally takes 8.90 s versus 25.38 s.

Completed since the previous checkpoint: dense and sparse matrix storage,
views/callbacks/reducers, state and system operations, constraint collections,
solver construction and option overrides, linear/controller replacement hooks,
and finite-difference Jacobian validation. An external process/constraint
callback interface is implemented; its numerical integration cases remain.

Remaining areas:

- Rate-store and process-set numerical cases, grouped views.
- External models and constraints, aerosol and DAE integration.
- Analytical chemistry integration, Chapman, regression, terminator and OpenMP.
- Portable build audit, malformed-input tests, more TS1 seeds/group sizes and
  final architecture/extensibility documentation.
