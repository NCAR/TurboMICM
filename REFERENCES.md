# References

Every source used in Tier 3, with how it was read and what was taken from it.
No code was copied from any source; the methods were implemented from the
descriptions and tables cited here. "Read" means the author read the text
named; "via review" means the source was read and reported by a research
assistant (`research/`), and the statement used from it is the one quoted in
that report; where both apply, the author checked the point used.

## What the literature settled

- **The accuracy target.** The solver benchmarks for atmospheric chemistry
  were built around "modest accuracies, say relative accuracies near 1 %"
  [S97a], the two-digit level of this experiment; Zhang et al. argue for
  0.1 % so that numerical error stays one to two orders below the 1 %
  overall target [Z11]. Operational models run between the two: GEOS-Chem
  at relative tolerances from 5e-4 to 5e-2 [GC], CAM-chem's own solver
  without an error estimate at all [CAM6]. Rate-constant uncertainties of
  5 to 30 % [NE17] and operator-splitting errors of up to 10 % [SZY16, PMK16]
  bound what a solver's error needs to be. ACCURACY.md builds on these.
- **The method family.** Third-order Rosenbrock methods (ROS3, RODAS3) were
  the most cost-effective to about three digits in the 1997 benchmarks, the
  fourth-order RODAS beyond, and SEULEX-type extrapolation two to three
  times slower at one to two digits [S97b]; KPP's integrators are "selected
  to be very efficient in the low to medium accuracy regime" [KPP]. ROS2 is a
  W-method that ran accurately with fixed 10 to 20 minute steps [VSBH99].
- **Restarts and the step.** With frequent restarts "all solvers spend most
  of their time in the start-up phase" [S97b], and KPP advises using the
  last proposed step as the next start [KPP]: the warm start.
- **Controllers.** KPP's embedded estimate overestimates the local error;
  a safety factor of about 1.4 and the H211b filter cut function evaluations
  by 24 to 43 % in box and global runs [D25].
- **Factorization reuse.** W-methods keep their order with any matrix in
  place of the Jacobian [SW79]; ROS34PW2 is a third-order one [RA05, B16].
  Time-lagged Jacobians were tested in 1980 on test problems, where accuracy
  fell when the Jacobian was frozen and the gain depended on the Jacobian and
  factorization dominating the cost [VS83]; no chemistry study since 1999
  measured it, which this entry's scans do.

## Methods and tableaus

- **[HW96]** E. Hairer and G. Wanner, *Solving Ordinary Differential
  Equations II: Stiff and Differential-Algebraic Problems*, 2nd ed.,
  Springer, 1996. Standard text; the transformed Rosenbrock form with the
  stage vectors u = Gamma k (section IV.7, which MICM's parameter sets
  follow and `from_standard_form` in `src/parameters.c` reproduces), the
  predictive step-size controller as RODAS and RADAU5 use it (IV.8), and the
  extrapolation of the linearly implicit Euler method with the Aitken-Neville
  scheme in powers of h (IV.9). Known to the author; not re-read for this
  work.
- **[S97b]** A. Sandu, J. G. Verwer, J. G. Blom, E. J. Spee, G. R. Carmichael,
  F. A. Potra, "Benchmarking stiff ODE solvers for atmospheric chemistry
  problems II: Rosenbrock solvers", *Atmos. Environ.* 31(20), 3459–3472,
  1997, doi:10.1016/S1352-2310(97)83212-8. Via review (read in full there).
  Taken: the origin of ROS3 and RODAS3 (MICM's three-stage and four-stage
  DAE sets); the benchmark's conclusions (Rosenbrock most cost-effective to
  four digits; RODAS3 and ROS3 best to three; RODAS beyond); the cost of
  restarts; relative tolerances of 0.1 and above make Rosenbrock solutions
  drift; absolute tolerance 0.01 molecules cm-3; the step rule and its
  constants.
- **[VSBH99]** J. G. Verwer, E. J. Spee, J. G. Blom, W. Hundsdorfer, "A
  second-order Rosenbrock method applied to photochemical dispersion
  problems", *SIAM J. Sci. Comput.* 20(4), 1456–1480, 1999,
  doi:10.1137/S1064827597326651. Via review (read there). Taken: ROS2 with
  gamma = 1 + 1/sqrt 2 is second-order consistent for any matrix in place of
  the Jacobian (a W-method), stable with large fixed steps (10 to 20
  minutes), its stages clipped at zero without notable loss of accuracy;
  the first-order embedding is the intermediate stage. The standard-form
  coefficients used in `prototypes/tableau_check.py` (alpha_21 = 1,
  gamma_21 = -2 gamma, b = (1/2, 1/2)) were checked by the author against
  MICM's transformed set.
- **[RA05]** J. Rang and L. Angermann, "New Rosenbrock W-methods of order 3
  for partial differential algebraic equations of index 1", *BIT Numer.
  Math.* 45(4), 761–787, 2005, doi:10.1007/s10543-005-0035-y. Abstract and
  secondary descriptions only (the paper is paywalled). The source of
  ROS34PW2: four stages, order 3 with an embedded order 2, L-stable, order
  kept with an approximate Jacobian.
- **[B16]** D. S. Blom, P. Birken, H. Bijl, F. Kessels, A. Meister, A. H. van
  Zuijlen, "A comparison of Rosenbrock and ESDIRK methods combined with
  iterative solvers for unsteady compressible flows", *Adv. Comput. Math.*
  42(6), 1401–1426, 2016, doi:10.1007/s10444-016-9468-x. Read by the author
  (the authors' manuscript): Appendix A, Table 7, gives ROS34PW2's
  coefficients in the standard form (gamma, alpha_ij, gamma_ij, b, bhat),
  which `tmc_rosenbrock_ros34pw2` converts; the paper also found ROS34PW2
  more robust than the ESDIRK methods it compared.
- **[PETSc]** PETSc manual page TSROSWRA34PW2, petsc.org/main/manualpages/TS,
  read by the author 2026-10-18: ROS34PW2 is strongly A-stable with
  R(infinity) = 0, the embedded method has R(infinity) = 0.48, and only an
  approximate Jacobian is needed.
- **[J19]** T. Jax, *A rooted-tree based derivation of ROW-type methods with
  non-exact Jacobian entries for index-one DAEs*, PhD thesis, University of
  Wuppertal, 2019, doi:10.25926/e36a-vq51. Via review (tables and theorems
  read there). Taken: the W-orders of ROS2 (2), ROS34PW2 (3), ROS3P (1) and
  others, so which tableaus may reuse a factorization without losing order.
- **[SW79]** T. Steihaug and A. Wolfbrandt, "An attempt to avoid exact
  Jacobian and nonlinear equations in the numerical solution of stiff
  differential equations", *Math. Comp.* 33(146), 521–534, 1979. Abstract
  only. The origin of W-methods: a Jacobian kept over several steps.
- **[VS83]** J. G. Verwer and S. Scholz, "Rosenbrock methods and time-lagged
  Jacobian matrices", *Beitr. Numer. Math.* 11, 173–183, 1983 (CWI preprint
  NW 82/80, 1980). Via review. Taken: with a frozen Jacobian accuracy fell on
  test problems, and the gain is worth having only where the Jacobian and
  its factorization dominate the cost, which the cost model in DESIGN.md
  shows they do here (46 % of a step).
- **[TS14]** P. Tranquilli and A. Sandu, "Rosenbrock-Krylov methods for large
  systems of differential equations", *SIAM J. Sci. Comput.* 36(3),
  A1313–A1338, 2014. Via review. Taken: non-W Rosenbrock methods lose order
  with an inexact Jacobian (ROS4 and RODAS4 fell from 4 to about 3 in their
  tests, a W-method kept its order).
- **[D85]** P. Deuflhard, "Recent progress in extrapolation methods for
  ordinary differential equations", *SIAM Review* 27(4), 505–535, 1985. Known
  to the author; the extrapolation of the linearly implicit Euler method, the
  harmonic step-number sequence, and order control by work per unit step.
- **[U22]** Utkarsh, C. Elrod, Y. Ma, K. Althaus, C. Rackauckas,
  "Parallelizing explicit and implicit extrapolation methods for ordinary
  differential equations", arXiv:2207.08135, 2022. Via review. Taken: the
  sequences in use (harmonic, Romberg, Bulirsch) and their work; serial
  extrapolation is not competitive with Rosenbrock methods on standard stiff
  benchmarks, so the extrapolation arm is measured, not assumed.
- **[KC16]** C. A. Kennedy and M. H. Carpenter, "Diagonally implicit
  Runge–Kutta methods for ordinary differential equations. A review",
  NASA/TM-2016-219173, 2016. Via review. Taken: the stage-order limit of 2
  and the Newton-iteration cost that put SDIRK methods behind Rosenbrock
  methods at this accuracy; SDIRK was not pursued.
- **[C22]** N. J. Curtis, K. E. Niemeyer, C.-J. Sung, "Accelerating
  reactive-flow simulations using vectorized chemistry integration",
  arXiv:2205.05681, 2022. Via review. Taken: vectorized ROS3, ROS4, RODAS3
  and RODAS4 perform alike at tolerances 1e-4 to 1e-6 on combustion
  kinetics, the fourth-order pair pulling ahead only at tighter tolerances.
- **[SR97]** L. F. Shampine and M. W. Reichelt, "The MATLAB ODE Suite",
  *SIAM J. Sci. Comput.* 18(1), 1–22, 1997. Via review. Taken: the
  Rosenbrock code ode23s is "especially effective at crude tolerances, when
  a one-step method has advantages".

## Step-size control, norms, and reuse rules

- **[G94]** K. Gustafsson, "Control-theoretic techniques for stepsize
  selection in implicit Runge–Kutta methods", *ACM TOMS* 20(4), 496–517,
  1994. Known to the author through [HW96]; the predictive controller
  h_new = h (h / h_prev) (err_prev / err^2)^(1/k) times the safety factor,
  taken as the smaller of it and the elementary factor, as RADAU5 and RODAS
  do (`controller=predictive`).
- **[So03]** G. Söderlind, "Digital filters in adaptive time-stepping",
  *ACM TOMS* 29(1), 1–26, 2003, and G. Söderlind, "Automatic control and
  adaptive time-stepping", *Numer. Algorithms* 31, 281–310, 2002. Via
  review (the formulas from the author's 2018 lecture slides, read there).
  Taken: the H211b filter h_new / h = (1/err)^(1/(bk)) (1/err_prev)^(1/(bk))
  (h/h_prev)^(-1/b), stable for b >= 1 (`controller=h211b`, `h211_b`,
  `h211_k`).
- **[D25]** R. Dreger, T. Kirfel, A. Pozzer, S. Rosanka, R. Sander,
  D. Taraborrelli, "Optimized step size control within the Rosenbrock
  solvers for stiff chemical ordinary differential equation systems in KPP
  version 2.2.3_rs4", *Geosci. Model Dev.* 18, 4273–4291, 2025,
  doi:10.5194/gmd-18-4273-2025. Read by the author (the article's HTML) and
  via review. Taken: KPP's first-order controller and its constants
  (rejection factor 0.1, start step 1e-5 s); the embedded estimate nearly
  always overestimates the local error, so a safety factor of about 1.4
  (1.3 to 1.5) pays; H211b with b = 1 and k = 1.7 cut function evaluations
  by up to 43 % (gas phase), 27 % (cloud) and 13 % (aerosol) in EMAC and
  saved 11.4 % of a year's core hours; the tolerances in use (relative
  1e-2, absolute 10 molecules cm-3 for gas-phase chemistry); the controller
  memory is initialized neutrally at each solve. These set the `safety`,
  `h211_b` and `h211_k` options and the controller scan.
- **[KPP]** KPP documentation, kpp.readthedocs.io (versions 3.3.0 to
  3.5.0): "Numerical methods", "Integrator options", "Integrator status".
  Read by the author (3.3.0's methods page) and via review. Taken: the
  Rosenbrock methods' stages, function calls and orders (ROS-2 2/2/2(1),
  ROS-3 3/2/3(2), ROS-4 4/3/4(3), RODAS-3 4/3/3(2), RODAS-4 6/5/4(3)); the
  accuracy regime they target (relative errors 1e-2 to 1e-5); the advice
  "for multiple restarts, use Hnew as Hstart in the subsequent run".
- **[SS06]** A. Sandu and R. Sander, "Technical note: Simulating chemical
  systems in Fortran90 and Matlab with the Kinetic PreProcessor KPP-2.1",
  *Atmos. Chem. Phys.* 6, 187–195, 2006. Via review. Taken: KPP's
  integrator set and its sparse-LU design; its benchmark timings put Rodas3
  first and Ros2 near last on a 212-species mechanism.
- **[CV]** SUNDIALS CVODE documentation v7.8.0, "Mathematical
  considerations". Via review (verified against the page there). Taken: the
  BDF codes' reuse heuristics (the Newton matrix is refreshed after 20 steps
  or a 30 % change of the step-dependent factor, the Jacobian after 50
  steps), the model for the age and change limits of `reuse_age` and
  `reuse_change`.
- **[Sh05]** L. F. Shampine, "Error estimation and control for ODEs",
  *J. Sci. Comput.* 25, 3–16, 2005. Via review. Taken: the error-per-step
  convention and the conservatism of step changes; the scaled RMS norm as
  the standard choice.
- **[L09]** J. C. Linford, J. Michalakes, M. Vachharajani, A. Sandu,
  "Multi-core acceleration of chemical kinetics for simulation and
  prediction", *Proc. SC09*, 2009 (and the first author's 2009 slides). Via
  review. Taken: a vectorized Rosenbrock solver takes the largest error over
  the lanes it integrates together, the `norm=cellmax` option.

## Accuracy requirements and solver comparisons

- **[S97a]** A. Sandu, J. G. Verwer, M. van Loon, G. R. Carmichael,
  F. A. Potra, D. Dabdub, J. H. Seinfeld, "Benchmarking stiff ODE solvers
  for atmospheric chemistry problems I: implicit vs explicit", *Atmos.
  Environ.* 31(19), 3151–3166, 1997. Via review (read in full there,
  equations checked on page images). Taken: the definition of significant
  digits of accuracy (the root mean square over time of the relative error,
  species below a threshold excluded; SDA 2 is 1 %); "a numerical comparison
  should focus on modest accuracies, say relative accuracies near 1 %.
  Higher accuracy levels are redundant for the actual practice of air
  pollution modelling".
- **[VBLS96]** J. G. Verwer, J. G. Blom, M. van Loon, E. J. Spee, "A
  comparison of stiff ODE solvers for atmospheric chemistry problems",
  *Atmos. Environ.* 30(1), 49–58, 1996. Via review. Taken: the same 1 %
  focus, and that sparsity pays from about 60 species on.
- **[E09]** P. Eller, K. Singh, A. Sandu, K. Bowman, D. K. Henze, M. Lee,
  "Implementation and evaluation of an array of chemical solvers in the
  Global Chemical Transport Model GEOS-Chem", *Geosci. Model Dev.* 2, 89–96,
  2009, doi:10.5194/gmd-2-89-2009. Via review. Taken: KPP's Rodas3 at
  relative tolerances 1e-1 to 1e-3 reached about 1.1 to 2.6 digits on ozone
  and needed about half the time of SMVGEAR II at equal accuracy; nine
  function calls and three steps per 15-minute advection step and cell;
  "we normally use tolerances of 1e-3 or lower for GEOS-Chem".
- **[H07]** D. K. Henze, A. Hakami, J. H. Seinfeld, "Development of the
  adjoint of GEOS-Chem", *Atmos. Chem. Phys.* 7, 2413–2433, 2007. Via
  review. Taken: GEOS-Chem's Rosenbrock solver was run "at tolerance levels
  that yield an SDA of about 1.0".
- **[Z11]** H. Zhang, J. C. Linford, A. Sandu, R. Sander, "Chemical mechanism
  solvers in air quality models", *Atmosphere* 2(3), 510–532, 2011,
  doi:10.3390/atmos2030510. Via review. Taken: the target of 0.1 % ("about
  3 accurate digits") so that numerical errors sit one to two orders below
  the overall 1 % accuracy; the taxonomy of model errors.
- **[GC]** GEOS-Chem Classic documentation (geos-chem.readthedocs.io, 14.6
  to 14.8) and the v11-02 wiki page. Via review. Taken: the relative
  tolerance moved from 2e-1 to 5e-4 across versions for robustness, with an
  absolute tolerance of 1e-2 molecules cm-3; the default solver is Rodas3
  (now RODAS3.1 with auto-reduction).
- **[S20, S22, L23]** L. Shen et al., *Geosci. Model Dev.* 13, 2475–2486,
  2020; L. Shen et al., *Geosci. Model Dev.* 15, 1677–1687, 2022; H. Lin et
  al., *J. Adv. Model. Earth Syst.* 15, e2022MS003293, 2023. Via review.
  Taken: the errors accepted for mechanism reduction in GEOS-Chem (median
  0.2 to 2 %, NO2 up to 6 to 10 % in places), as the scale of change a
  production model tolerates.
- **[D25]** as above: the gas-phase tolerances in use.

## The target model and the size of other errors

- **[Em20]** L. K. Emmons et al., "The Chemistry Mechanism in the Community
  Earth System Model Version 2 (CESM2)", *J. Adv. Model. Earth Syst.* 12,
  e2019MS001882, 2020. Via review. Taken: the lineage and size of TS1 (T1
  with stratospheric chemistry; the counts in TS1's MICM configuration are
  209 solver species and 547 reactions).
- **[CAM6]** CAM6 Scientific Guide, section 6.4 (the chemistry solver). Via
  review (read verbatim there). Taken: CAM-chem integrates TS1 with a
  first-order backward Euler and Newton iteration at a fixed 30-minute step,
  converged to a relative change of 1e-3 (1e-4 for some species), halving
  the step on non-convergence up to five times, without a local error
  estimate; long-lived species by forward Euler.
- **[Sch18]** M. G. Schultz et al., "The chemistry–climate model
  ECHAM6.3-HAM2.3-MOZ1.0", *Geosci. Model Dev.* 11, 1695–1723, 2018. Via
  review. Taken: the same preprocessor-generated solver's iteration counts
  (2 to 6, most at sunrise and sunset) and thresholds.
- **[Su17]** J. Sun, J. S. Fu, J. Drake, J.-F. Lamarque, S. Tilmes, F. Vitt,
  "Improvement of the prediction of surface ozone concentration over
  conterminous U.S. by a computationally efficient second-order Rosenbrock
  solver in CAM4-Chem", *J. Adv. Model. Earth Syst.* 9, 482–500, 2017. Via
  review (abstract). Taken: replacing the first-order solver by ROS-2 at the
  same 1800 s step changed summer surface ozone by up to 5.2 ppb at about
  half the cost: the operational solver's own error is of that order.
- **[NE17]** B. Newsome and M. Evans, "Impact of uncertainties in inorganic
  chemical rate constants on tropospheric composition and ozone radiative
  forcing", *Atmos. Chem. Phys.* 17, 14333–14352, 2017. Via review. Taken:
  rate-constant uncertainties of 5 to 30 % at 298 K, 10 % for photolysis
  rates, and the resulting 10 to 16 % uncertainties in ozone and OH.
- **[SZY16]** M. Santillana, L. Zhang, R. Yantosca, "Estimating numerical
  errors due to operator splitting in global atmospheric chemistry models:
  transport and chemistry", *J. Comput. Phys.* 305, 372–386, 2016. Via
  review. Taken: splitting strategies differ by at most about 10 %.
- **[PMK16]** S. Philip, R. V. Martin, C. A. Keller, "Sensitivity of
  chemistry-transport model simulations to the duration of chemical and
  transport operators: a case study with GEOS-Chem v10-01", *Geosci. Model
  Dev.* 9, 1683–1695, 2016. Via review. Taken: operator durations of 10 to
  60 minutes change NOx and secondary inorganic aerosol by 5 to 35 %, ozone
  and CO by under 2 %.

## Positivity

- **[VHB02]** J. G. Verwer, W. Hundsdorfer, J. G. Blom, "Numerical time
  integration for air pollution models", *Surv. Math. Ind.* 10, 107–174,
  2002. Via review. Taken: no Rosenbrock method of order above one is
  positive for all steps; with local error control, negative values lead to
  rejections; clipping is the practice, and splitting restarts create
  transients that a one-step method handles.
- **[Sa01]** A. Sandu, "Positive numerical integration methods for chemical
  kinetic systems", *J. Comput. Phys.* 170, 589–602, 2001, and D. Torlo,
  P. Öffner, H. Ranocha, "Issues with positivity-preserving Patankar-type
  schemes", *Appl. Numer. Math.* 182, 117–147, 2022. Via review (abstracts
  and sections). Considered and not adopted: projection costs an
  optimization per step, and Patankar-type schemes lose order near zero
  concentrations and oscillate for large steps.

## Reduced precision (considered for a scan)

- **[SKTC26]** K. Sophocleous, T. Kirfel, D. Taraborrelli, T. Christoudias,
  "GPU-accelerated atmospheric chemical kinetics using single precision in
  EMAC v2.55 with MEDINA v2.0", EGUsphere preprint egusphere-2026-2389, 2026
  (not peer reviewed), and K. Sophocleous, T. Christoudias, "Reduced-
  precision chemical kinetics in atmospheric models", *Atmosphere* 13(9),
  1418, 2022. Via review. Taken: KPP's Rosenbrock solvers in single
  precision kept species masses within 0.004 % of double precision over a
  month, with a 1.4x kernel speedup on a GPU whose single-precision peak is
  twice its double; the note that species depleted to the single-precision
  range differ.
- **[CH18]** E. Carson and N. J. Higham, "Accelerating the solution of
  linear systems by iterative refinement in three precisions", *SIAM J. Sci.
  Comput.* 40(2), A817–A847, 2018. Via review (abstract). Taken: the
  condition-number limits under which a low-precision factorization with
  higher-precision residuals converges.

## Not used

Exponential integrators (Hochbruck, Ostermann, Schweitzer 2009; Tokman;
the chemistry tests of Curtis et al. 2017 and Stewart et al. 2023), QSSA
schemes (Jay et al. 1997), TWOSTEP (Verwer 1994), EBI (Hertel et al. 1993),
SMVGEAR (Jacobson and Turco 1994), and multirate and partitioned methods
(Günther and Sandu 2016; Santillana et al. 2010) were reviewed
(`research/R3`, `research/R5`) and not pursued: the first because a direct
sparse factorization is cheap here, the next three because their accuracy
tops out near the two-digit level or they are mechanism-specific, the last
because partitioning changes the system solved, which the rules exclude.
