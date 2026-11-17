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
seed 2026 and 10 threads took 25,408 ms, compared with 26,576 ms for the MICM
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

Precomputing each sparse row and column's diagonal position reduced a
single-thread 1,000-cell, 30-step clang run from 6,263 ms to 4,128 ms in a
repeat measurement. The before and after outputs were identical in an
independent run, and the full 10,000-cell comparison above passed. These are
single-run timings; they are not the experimenter's ranking measurements.

All three tutorial configurations are compared against hand-built C
mechanisms. The checks cover generated rate constants, forcing, Jacobians,
and integration trajectories; the two rate-constant tutorials also reproduce
the reference programs' printed output.
The 100-cell TS1 output also passed with 1 thread, 3 threads, and the default
thread count. The case inventory contains 680 unique reference IDs: 675 are
mapped to passing executables and 5 have explicit N/A reasons, with no missing
or extra IDs.
As a diagnostic, generating the reactions in the reference's rate-type order
and using its exact state permutation still left 137 of 211 warm-up values
outside the tolerance on the all-default photolysis stress input (gcc-15 C
port versus clang reference). This experiment did not change the production
ordering and indicates that ordering alone does not explain the divergence.
Three additional 100-cell inputs from the supplied generator, using seeds 42,
12345, and 753, also passed against separately run references. The largest
relative differences were 2.74e-14, 2.54e-14, and 7.66e-14, respectively.
These used gcc-15 for the C port, with three threads, and the clang reference
with one thread.

Static, dynamic, guided, and default OpenMP schedules produced identical
per-cell output in the scheduling check.
