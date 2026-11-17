# Results

The work is incomplete. At the latest full case audit, `make test-cpu` emitted
680 case lines under both clang and gcc-15: 326 PASS, 350 FAIL because the case
has not yet been ported, and 4 N/A for C++ template-concept/SFINAE checks with
no C equivalent. The analytic, configuration, constraint, and Doolittle core
checks also passed under both compilers. A failing case report exits nonzero.
The analytical cases include 999 one-second steps over three cells for each
of six rate-law families, stiff variants, both time integrators, and the
surface-reaction case. The six-stage coefficients were checked against the
reference after the branched-rate trajectory exposed transcription errors.
The Robertson, E5, HIRES, and Oregonator trajectories also match the
reference checks, including backward Euler's 18,000 substeps per Oregonator
output interval.

With `-O3 -DNDEBUG`, TS1, seed 2026, 10,000 cells, 100 timed 30-second steps,
64-cell groups, and 10 threads, one local measurement was:

| Implementation | elapsed_ms |
| --- | ---: |
| C benchmark, clang | 8,430.915 |
| MICM reference, clang++ | 23,612.070 |

Both reported 141,300 internal steps, 141,300 accepted, no rejected steps,
282,600 function calls, 141,300 Jacobian updates and decompositions, and
423,900 solves. The supplied comparator passed all 2,110,000 output values;
the largest relative difference was 2.51e-13. These are individual local
measurements, not the experimenter's median ranking runs. The timing was
recorded before the later solver callback hooks and should be repeated after
subsequent changes.

For the supplied 100-cell seed-1 input, 100 steps, four threads, and groups of
32, the generated ordered mechanism passed the same comparison with clang
(largest relative difference 4.18e-14) and gcc-15 (6.04e-14).
