# Sources and use

No external implementation was downloaded or incorporated. New integrators and
experiment tooling were written here. Existing Tier 1 MICM attribution is in
NOTICE; MICM source e15898ce and its parameter tables remain the baseline.

1. A. Sandu, J. G. Verwer, J. G. Blom, E. J. Spee, G. R. Carmichael and F. A.
   Potra (1997), *Benchmarking stiff ODE solvers for atmospheric chemistry
   problems II: Rosenbrock solvers*, Atmospheric Environment 31, 3459–3472.
   [Paper](https://ir.cwi.nl/pub/125/0125D.pdf).
   Used the work–precision experiment design, thresholded significant-digit
   accuracy concept, emphasis on roughly 1% chemistry accuracy, and the warning
   that local tolerance is not achieved global accuracy. The supplied score.py
   uses pooled RMS plus species requirements; the paper's aggregate uses the
   worst species. We retain the supplied metric rather than equating them.
2. KPP numerical methods documentation,
   [Rosenbrock methods](https://kpp.readthedocs.io/en/stable/num_methods/rosenbrock-methods.html)
   and [integrator options](https://kpp.readthedocs.io/en/stable/num_methods/controlling_the_integrator.html),
   accessed 2026-10-18. Used the stage equations, method families, orders,
   stability classifications, and distinction between autonomous and explicitly
   time-dependent forcing. Existing C coefficient tables are retained from MICM.
3. E. Hairer and C. Lubich (1988), *Extrapolation at stiff differential
   equations*, Numerische Mathematik 52, 377–400,
   [abstract](https://archive-ouverte.unige.ch/unige:12472).
   Motivates testing extrapolation of linearly implicit Euler and treating stiff
   order effects cautiously. Only the abstract was accessible; no claims about
   its detailed algorithm are made. Our fixed 1,2,3-substep Richardson tableau
   is derived directly and is not an implementation of SEULEX.
4. E. Hairer and C. Lubich, *Numerical solution of ordinary differential equations*,
   [overview](https://www.unige.ch/~hairer/preprints/pcam-ode.pdf).
   Used the description of linearly implicit Euler as one Newton iteration,
   L-stability, and stiff order reduction to frame the experiments.
5. SUNDIALS CVODE user documentation,
   [linear solver setup and Jacobian reuse](https://sundials.readthedocs.io/en/v7.5.0/cvode/Usage/),
   accessed 2026-10-18. Used the distinction between reusing Jacobians and
   reusing factorizations, refresh after failure, and bounded reuse age. Our
   implementation only reuses a factorization with exactly unchanged step and
   Jacobian. Reusing an old Jacobian in an ordinary Rosenbrock table does not
   establish W-method order conditions; it is an experimental approximation
   assessed separately, not a claim of retained formal order.
6. G. Söderlind (2002), *Automatic control and adaptive time-stepping*,
   Numerical Algorithms 31, 281–310,
   [author's paper](https://citeseerx.ist.psu.edu/document?doi=f5fd153e99b16dd4b5348bcbde0f2121b7096db8&repid=rep1&type=pdf).
   Used the indexed excerpt of equation (4.21)'s PI.3.4 control structure: current error exponent
   -0.7/q and previous error exponent 0.4/q, with an I-controller safety net
   after rejection. Our multiplicative safety 0.9 and clamps are explicitly
   configurable choices, not coefficients copied from an implementation.
7. KPP 3.5 numerical methods documentation,
   [RODAS3.1 coefficient table](https://kpp.readthedocs.io/en/latest/num_methods/rosenbrock-methods.html#rodas-3-1),
   accessed 2026-10-18. Implemented the published mathematical coefficients as
   the additional `rodas31` method. The documentation attributes them to Long et al. (2026, submitted); that
   manuscript was not used.
   No KPP implementation code was read or copied. Independent linear and
   nonlinear analytic convergence tests verify order three with these numbers.
8. R. Dreger, T. Kirfel, A. Pozzer, S. Rosanka, R. Sander and D. Taraborrelli
   (2025), *Optimized step size control within the Rosenbrock solvers for stiff
   chemical ordinary differential equation systems in KPP version 2.2.3_rs4*,
   Geoscientific Model Development 18, 4273–4291,
   [paper](https://gmd.copernicus.org/articles/18/4273/2025/).
   Used equation (6)'s H211b controller, independent sweeps of controller
   parameters, and the observation that embedded errors can overestimate actual
   integration errors. We retain the acceptance test; only the proposed next
   step changes. No pseudocode or implementation was copied. Tests include
   configurable b and k; the paper motivates exploring b=1 and k near 2 for
   a third-order method rather than assuming an error exponent is universal.
9. A. Sandu et al. (1997), *Benchmarking stiff ODE solvers for atmospheric
   chemistry problems I: implicit vs explicit*, Atmospheric Environment 31,
   3151–3166, [primary repository abstract](https://authors.library.caltech.edu/records/wbx07-tym26).
   Used the abstract's motivation for sparse implicit solvers around two digits
   of accuracy. Detailed method comparisons rely on part II above, not on an
   inaccessible full text of part I.
