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
4. Alexander (1977), *Diagonally Implicit Runge–Kutta Methods for Stiff O.D.E.'s*,
   SIAM J. Numerical Analysis 14, 1006–1021,
   https://doi.org/10.1137/0714068, and John Butcher's
   [Implicit Runge–Kutta tutorial](https://www.math.auckland.ac.nz/~butcher/ODE-book-2008/Tutorials/IRK.pdf).
   Used the stiffly accurate two-/three-stage SDIRK family and third-order
   diagonal root. The embedded weights were obtained here from the quadrature
   moment equations, with the last embedded weight fixed to zero. Newton
   iteration, stopping rules, and error/reuse experiments are local implementations.


5. Ernst Hairer and Christian Lubich, *Numerical solution of ordinary
   differential equations*, author-hosted overview, sections 3.3 and 4.5.
   [Author PDF](https://www.unige.ch/~hairer/preprints/pcam-ode.pdf).
   Used in the final documentation review to relate the implemented construction
   to Richardson extrapolation and linearly implicit Euler. The substep sequence,
   fixed-Jacobian policy, error estimator and controller are specified and tested
   here; no external implementation was consulted or copied.
6. Hairer and Lubich (1988), *Extrapolation at stiff differential equations*,
   Numerische Mathematik 52, 377–400, DOI: 10.1007/BF01462235.
   [Author repository abstract](https://archive-ouverte.unige.ch/unige:12472).
   Only the public abstract was consulted; full text requires access. Its account
   of stiffness-dependent perturbations in Euler extrapolation error expansions
   supports the distinction between formal nonstiff order and measured stiff
   accuracy. It is not used as a proof for this implementation.

Read-only local sources: `../reference/README.md`, `../reference/micm_cases.txt`,
`../reference/compare_final.py`, `../scenarios/README.md`,
`../scenarios/score.py`, and the supplied scenario/reference CSV files. These
define the compatibility gate, forcing convention, and external accuracy metric.
The supplied `micm_scenario.cpp` documents the reference's default tolerances.
