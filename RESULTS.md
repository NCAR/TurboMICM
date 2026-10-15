# Validation and measurements

Work in progress. This is not yet a complete port of the required case suite.

Initial clang build, `-O3 -DNDEBUG`, 100 cells, 100 steps, 30 seconds, seed 1, one thread, group size 32:

- Gate: PASS, all 21,100 output values; largest relative difference from the supplied reference 9.69e-13.
- Elapsed: 6501.244 ms.
- Timed solver totals: 3600 internal/accepted steps, no rejections, 7200 forcing calls, 3600 Jacobian updates/decompositions, 10800 solves. Four independent cell groups account for the totals.

Further compiler, case-suite, and performance checks remain required.
