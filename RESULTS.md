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

An additional one-cell input omitted temperature and pressure. With the six
unphysical photolysis defaults set to zero, the generated condition defaults
agree with the reference within 3.44e-15 relative after two timed steps. When
those six rates are also omitted and therefore retain their reference defaults
of 1 s-1, the comparison is less close: 44 of 211 values exceed the gate
tolerance after two timed steps, with a largest relative difference of
1.47e-5. This stiff input is outside the supplied ranking inputs, which set
those six rates to zero; its discrepancy remains under investigation.
A warm-up-only check of that same input passes at step lengths of 1e-12,
1e-6, and 1e-3 seconds, but exceeds the tolerance at 1 and 30 seconds. This
locates the divergence during integration at larger steps; it does not yet
identify whether arithmetic conditioning or an implementation difference is
responsible.
At a 30-second warm-up, setting just one of the six photolysis defaults to
1 s-1 while setting the other five to zero produces out-of-tolerance values
for `jo2_a`, `jh2o_b`, `jh2o_c`, and `jch4_b`; the corresponding `jch4_a` and
`jhf` trials pass. This narrows the stress case without establishing its cause.
The reference solution itself is very sensitive on that stress input. With
`OH=1e-10` and the default photolysis rates, increasing the initial H2O value
from 0.746 to the next representable double (0.74600000000000011) changes 128
of 211 reference output values beyond the gate tolerance after the 30-second
warm-up; the largest relative change is 2.84e-5. A one-ULP O2 increase changes
104 values beyond tolerance, with a largest relative change of 4.63e-6. This
shows that roundoff-sized input changes can be amplified to the same scale as
the observed port/reference discrepancy, although it does not rule out an
implementation difference.
