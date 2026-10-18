# Sources and their use

No external implementation was downloaded or copied. New integration and
experiment code was written here. The retained Tier 1 code and coefficient
tables have their existing MICM attribution.

1. Sandu, Verwer, Van Loon, Carmichael, Potra, Dabdub and Seinfeld (1997),
   *Benchmarking stiff ODE solvers for atmospheric chemistry problems I:
   implicit vs explicit*, Atmospheric Environment 31, 3151–3166.
   DOI: https://doi.org/10.1016/S1352-2310(97)00059-9.
   [Author repository](https://authors.library.caltech.edu/records/wbx07-tym26).
   Used the reported competitiveness of sparse implicit methods around 1% error
   to motivate the method-family scan; not a universal accuracy guarantee.
2. Sandu, Verwer, Blom, Spee, Carmichael and Potra (1997),
   *Benchmarking stiff ODE solvers for atmospheric chemistry problems II:
   Rosenbrock solvers*, Atmospheric Environment 31, 3459–3472.
   [Paper](https://ir.cwi.nl/pub/125/0125D.pdf).
   Used the embedded Rosenbrock family comparison and atmospheric box-model
   work–precision methodology. Coefficient values themselves are the already
   retained MICM tables in `src/parameters.c`.
3. Verwer, Spee, Blom and Hundsdorfer (1999), *A second-order Rosenbrock method
   applied to photochemical dispersion problems*, SIAM J. Scientific Computing
   20, 1456–1480. [Paper](https://ir.cwi.nl/pub/1252/1252D.pdf).
   Used ROS2's second-order consistency with an approximate Jacobian to motivate
   controlled Jacobian/factor reuse experiments; the discussion of internal
   positivity motivated a separate stage-clipping experiment. Clipping is an
   explicitly measured numerical projection and does not guarantee conservation.

Read-only local sources: `../reference/README.md`, `../reference/micm_cases.txt`,
`../reference/compare_final.py`, `../scenarios/README.md`,
`../scenarios/score.py`, and the supplied scenario/reference CSV files. These
define the compatibility gate, forcing convention, and external accuracy metric.
The supplied `micm_scenario.cpp` documents the reference's default tolerances.
