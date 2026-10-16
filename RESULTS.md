# Verification in progress

These are intermediate results, not completion of the requested correctness gate.
The full 680-case suite is not yet ported.

With Clang and `OPTFLAGS=-O3 -DNDEBUG`:

- 32 rate-constant and lambda-override cases passed their ported assertions.
- The C TS1 run on the supplied 100-cell seed-1 input, 100 timed 30-second steps
  after one warm-up, passed all 21,100 values against the supplied reference.
  The largest relative difference was 9.7e-13; the gate is rtol=1e-10, atol=1e-30.
- The same check passed with one thread/groups of 32, four threads/groups of 7,
  and the generic table-driven chemistry evaluator.
- One preliminary generated-kernel run with one thread/groups of 32 took
  405.04 ms. This is a development observation, not a controlled benchmark or a
  ranking measurement. Timed statistics were 3,600 internal/accepted steps,
  zero rejected, 7,200 function calls, 3,600 Jacobians/decompositions, 10,800 solves.

The supplied C++ reference executable was built using `make micm-ref` with the
same optimization flags. Fresh comparative timing, the 10,000-cell run, complete
case coverage, and the portability audit remain outstanding.

With GCC 15 and the same flags, all 32 currently mapped cases passed. TS1 also
passed the supplied 100-cell gate with four threads and groups of seven; the
largest relative difference was 9.69e-13.

The subsequent constraint checkpoint passes 86 mapped cases with both compilers:
32 rate/lambda cases and 54 linear, equilibrium, set, and initialization cases.
These include rollback on failure, damped Newton cold starts, row scaling,
multi-cell indexing, analytic residuals/Jacobians, and finite-difference checks.
The TS1 benchmark has not yet been rerun after that checkpoint.

The matrix checkpoint passes all 420 implemented cases with Clang and GCC 15
under `-O3 -DNDEBUG`: 117 dense/vector, 190 sparse, 27 grouped-view/category,
and the previous 86 chemistry/constraint cases. These cover the original lane
widths, CSR/CSC indices, exact printed matrix output, native integer arrays,
mutable/const callback views, dimension errors, temporary storage, and reductions.
There are 260 cases still unported, and no N/A declarations.
