# Results (work in progress)

With `-O3 -DNDEBUG`, the reader tests, TS1 finite-difference Jacobian check,
single-cell Rosenbrock and backward Euler decay checks, and a generated TS1
step with both LU methods passed with Apple clang and gcc-15. Thirty-six
of the 680 named reference cases have been individually ported; the remaining
644 are reported as `FAIL` by `make test-cpu` rather than treated as passes.

For the supplied 100-cell seed-1 input, 100 timed 30-second steps after warm-up,
the clang build agreed with the supplied MICM final state in all 21,100 values
at the requested tolerance (largest relative difference 5.68e-14). The GCC build also passed after the species-ordering and clamp changes (largest
relative difference 4.26e-14).

A single 10-thread run with matching clang `-O3 -DNDEBUG` flags measured 868 ms
for this implementation and 208 ms for the MICM reference on that 100-cell
input. These exploratory timings are not a ranking measurement. The required
full CPU compatibility gate and representative large benchmark remain open.
