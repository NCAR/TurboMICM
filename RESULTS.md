# Results

With `-O3 -DNDEBUG`, both Apple clang and gcc-15 pass `make check`,
`make test-cpu`, and `make check-ts1`. The CPU report contains 675 PASS and
5 N/A entries. The N/A cases concern C++ template concepts and category
selection; the C interfaces use explicit view descriptors. The reference list
contains 680 unique IDs, with no missing or extra mappings.

For the supplied 100-cell, 100-step TS1 comparison, the largest relative
differences were 1.76e-14 with clang and 3.53e-14 with gcc-15, against the
requested 1e-10 relative and 1e-30 absolute tolerance. Each comparison covered
21,100 values, including temperature and pressure. Three additional 100-cell
inputs from the supplied generator, using seeds 42, 12345, and 753, also passed
against separately run references; the largest relative differences were
3.56e-14, 2.20e-14, and 4.72e-14. These used gcc-15 for the C port with three
threads and the clang reference with one thread.

A single clang run of the full 10,000-cell, 100-step, 30-second workload with
seed 2026 and 10 threads took 25,689 ms, compared with 26,324 ms for the MICM
reference built with the same optimization flags. These are single runs, not
median ranking measurements. The C benchmark advances cells independently;
the reference uses one group per thread, so their summed statistics count
different units. Both performed nine accepted internal steps per 30-second
step, with no rejected steps. All 2,110,000 output values passed the tolerance
check; the largest relative difference was 1.25e-13.

The 100-cell output passed with 1 thread, 3 threads, and the default thread
count. Static, dynamic, guided, and default OpenMP schedules produced
identical per-cell output in the scheduling check. The alternate diagonal
Markowitz generator ordering also passed the gcc-15 suite and TS1 comparison.

All three tutorial configurations are compared against hand-built C
mechanisms. The checks cover generated rate constants, forcing, Jacobians,
and integration trajectories. Tutorial executables are also compared with
reference output; timing output and the order of independent OpenMP thread
announcements are normalized because those fields are nondeterministic. The
vectorized matrix tutorial's reference solve loop starts converged and thus
updates rate parameters without advancing concentrations.

A one-cell input that omits temperature and pressure passes the comparison
when the six unphysical photolysis rates are set to zero. When those rates are
also omitted and retain their 1 s-1 defaults, the same-compiler clang port and
reference give identical output after a 30-second warm-up and after two timed
steps. The gcc-15 port differs from the clang reference by up to 3.73e-5
relative after the warm-up on this stiff input; this input is outside the
supplied benchmark gate, whose generator sets the six rates to zero. The
reference is itself highly sensitive there: changing the initial H2O value by
one representable double changes its output by up to 2.84e-5 relative. The
cross-compiler difference is therefore a reproducibility limitation on that
unphysical input, not evidence that the gate tolerance holds for it.

Precomputing sparse row and column diagonal positions reduced a single-thread
1,000-cell, 30-step clang run from 6,263 ms to 4,128 ms in a repeat measurement
using the alternate ordering. The before and after outputs were identical in
an independent run. These are single-run timings, not ranking measurements.
