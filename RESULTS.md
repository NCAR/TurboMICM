# Results

With `-O3 -DNDEBUG`, the reader tests, TS1 finite-difference Jacobian check,
single-cell Rosenbrock and backward Euler decay checks, and a generated TS1
step with both LU methods passed with Apple clang and gcc-15. All 676 applicable
named reference cases pass with both compilers; four template-specific C++
cases are N/A. The report contains one result for each of the 680 cases.

For the supplied 100-cell seed-1 input, 100 timed 30-second steps after warm-up,
the clang build agreed with the supplied MICM final state in all 21,100 values
at the requested tolerance (largest relative difference 4.18e-14). The GCC build also passed after the species-ordering and clamp changes (largest
relative difference 6.04e-14).
An input omitting temperature, pressure, and all photolysis columns also
matched the reference after three 1e-6-second steps; all 211 reported values
passed, with a largest relative difference of 2.88e-16.

A single 10-thread run on 10,000 cells, 100 timed 30-second steps, seed 2026,
with clang `-O3 -DNDEBUG` measured 30,496 ms for this implementation and
29,617 ms for the same-flag MICM reference. All 2,110,000 output values passed
the requested tolerance; the largest relative difference was 1.77e-13. These are
single exploratory samples, not ranking measurements.
