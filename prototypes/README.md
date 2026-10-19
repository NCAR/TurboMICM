# Prototypes

The Tier 3 methods in Python, over the C library's own forcing, Jacobian,
factorization, and solves (through ctypes), so that a method's accuracy and
operation counts can be studied before and beside its C implementation, and
so that the C implementation can be checked against an independent statement
of the same algorithm. Speed is never measured here.

    make proto                                   # build/libtmc_proto.{dylib,so}
    make micm-scenario                           # bin/micm_scenario, for small references
    .conda/bin/python prototypes/make_mini.py    # prototypes/mini: 4 cells, 12 daylight steps, a reference
    .conda/bin/python prototypes/scenario.py prototypes/mini/scenario.txt --method rodas3 --option rtol=1e-3 --option atol=1e-18
    .conda/bin/python prototypes/agreement.py    # every C method against its prototype

| File | Contents |
|---|---|
| `shim.c` | the library's pieces exported for ctypes, in row-major (cell, variable) arrays |
| `tmclib.py` | the Python side: `Library`, `State`, the parameter sets, MICM's loop constants |
| `methods.py` | `RowSolver` and `ExtrapolationSolver`, mirrors of `src/row.c` and `src/extrapolation.c`, with the same options |
| `scenario.py` | runs a method on a scenario's first cells, scores it as `score.py` does over those cells, counts operations, appends scan rows |
| `agreement.py` | the test: each method in C and in Python on the mini scenario, states to roundoff and equal counts |
| `make_mini.py` | makes the mini scenario and its MICM reference |

## What the prototypes share with C, and what they do not

- The forcing, the Jacobian, the sparse LU, and the solves are the library's,
  so the two implementations differ only in the vector arithmetic of the
  stage combinations, the error norm's summation order, and the Python
  loop's own bookkeeping. Agreement to 1e-13 to 1e-11 relative is the norm
  on the mini scenario; the extrapolation with the H211b controller reaches
  1.3e-9, since two fractional powers in the controller and four rows of
  Aitken-Neville differences amplify the last-bit differences, which the
  stiff system then carries into its most sensitive species. The test's
  bar is 1e-8 relative with equal operation counts; a step-acceptance
  decision at the margin (error within roundoff of 1) could flip and
  separate the two runs by more, which the test would report.
- The parameter sets come from the library (`proto_tableau`), so the
  coefficients exist in one place.
- `scenario.py` scores only the cells it ran, with score.py's definitions
  (species at or above 1e-14 mol m-3, root mean square of relative errors,
  the gate), against the full reference's rows for those cells.
