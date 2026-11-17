# Results

With `-O3 -DNDEBUG`, both Apple clang and gcc-15 pass `make check`,
`make test-cpu`, and `make check-ts1`. The CPU case report contains 675 PASS
and 5 N/A entries. The N/A cases concern C++ template concepts and category
selection; the C interfaces use explicit view descriptors instead.

For the 100-cell, 100-step TS1 comparison, the largest relative differences
were 4.18e-14 with clang and 6.04e-14 with gcc-15, against the requested
1e-10 relative and 1e-30 absolute tolerance. Each comparison covered 21,100
values, including the temperature and pressure columns.

A single clang run of the full 10,000-cell, 100-step, 30-second workload with
seed 2026 and 10 threads took 36,697 ms, compared with 26,479 ms for the MICM
reference built with the same optimization flags. These are single runs, not
median ranking measurements. The C benchmark advances cells independently;
the reference uses one group per thread, so their summed statistics count
different units. Both performed nine accepted internal steps per 30-second
step, with no rejected steps. All 2,110,000 output values passed the tolerance
check; the largest relative difference was 1.77e-13.

The tutorial executables are compared with output captured from the reference.
Timing output and the order of independent OpenMP thread announcements are
normalized because those fields are nondeterministic. The vectorized matrix
tutorial's reference solve loop starts in the converged state and therefore
updates rate parameters without advancing concentrations.
