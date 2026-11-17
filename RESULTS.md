# Results (work in progress)

With `-O3 -DNDEBUG`, the reader tests, TS1 finite-difference Jacobian check,
single-cell Rosenbrock and backward Euler decay checks, and a generated TS1
step with both LU methods passed with Apple clang and gcc-15. Seventy-seven
of the 680 named reference cases have been individually ported; the remaining
603 are reported as `FAIL` by `make test-cpu` rather than treated as passes.

For the supplied 100-cell seed-1 input, 100 timed 30-second steps after warm-up,
the clang build agreed with the supplied MICM final state in all 21,100 values
at the requested tolerance (largest relative difference 5.68e-14). The GCC build also passed after the species-ordering and clamp changes (largest
relative difference 4.26e-14).

A single 10-thread run on 10,000 cells, 100 timed 30-second steps, seed 2026,
with clang `-O3 -DNDEBUG` measured 30,021 ms for this implementation after
sparse-index optimizations. An earlier C run matched all 2,110,000 values from
the same-flag MICM reference at the requested tolerance (largest relative
difference 1.94e-13); the optimized output was bit-for-bit identical to that
earlier run. The reference sample took 28,015 ms. These are single exploratory
samples, not ranking measurements. The full CPU compatibility gate remains open.
