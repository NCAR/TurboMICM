# Work remaining

The full task is not complete. `make test-cpu` deliberately reports every missing
case as FAIL. The next implementation areas are:

- Dense/sparse matrix operations, views, callback validation and reducers.
- Species properties, phase/system metadata, state copy/set behavior, builder and
  rate-store cases, process-set numerical tests.
- Constraint-set variants, external process interfaces and constraint integration.
- Analytical chemistry integration, Chapman, aerosol, regression and remaining
  solver cases, retaining all original assertions and tolerances.
- Portable build verification, malformed-input tests, additional TS1 seeds/group
  sizes, performance comparison with matching MICM flags, final documentation.

The last verified state has 140/680 cases passing with both clang and gcc-15.
The complete provided TS1 fixture passes both compilers. No case is yet N/A.
