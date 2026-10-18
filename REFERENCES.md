# References and provenance

No external implementation code was downloaded or copied. Tier 1 source and its
MICM provenance are retained. New numerical code is written from equations.

- Sandu, A., Verwer, J. G., Blom, J. G., Spee, E. J., Carmichael, G. R., and
  Potra, F. A. (1997). *Benchmarking stiff ODE solvers for atmospheric chemistry
  problems II: Rosenbrock solvers*, Atmospheric Environment 31, 3459–3472.
  [Paper](https://ir.cwi.nl/pub/125/0125D.pdf).
  Used the rationale for sparse one-step methods and restart tests, the
  distinction between local tolerance and measured accuracy, and the warning
  that lagged Jacobians can reduce order. The supplied scoring script defines
  the precise accuracy aggregation used here.
- KPP 3.5 documentation, [Rosenbrock methods](https://kpp.readthedocs.io/en/stable/num_methods/rosenbrock-methods.html),
  read 2026-10-18. Used the mathematical RODAS3.1 table, attributed there to
  Long et al. (2026, submitted), and checked stage/order classifications.
  Older tables remain the full-precision Tier 1 tables; the documentation rounds
  several of their coefficients and is not used to replace them.
- Alexander, R. (1977). *Diagonally Implicit Runge–Kutta Methods for Stiff
  O.D.E.'s*, SIAM Journal on Numerical Analysis 14, 1006–1021.
  [Publisher](https://doi.org/10.1137/0714068),
  [INL mathematical description](https://tmap8.inl.gov/docs/PRs/292/site/source/timeintegrators/LStableDirk2.html).
  Used the two-stage, second-order L-stable SDIRK formula. The embedded estimate,
  chord iteration stopping rule, and safeguards here are explicitly described
  design choices, not an imported solver.
- Dreger, R., Kirfel, T., Pozzer, A., Rosanka, S., Sander, R., and
  Taraborrelli, D. (2025). *Optimized step size control within the Rosenbrock
  solvers for stiff chemical ordinary differential equation systems in KPP
  version 2.2.3_rs4*, GMD 18, 4273–4291.
  [Article](https://gmd.copernicus.org/articles/18/4273/2025/).
  Used equation (6) for the H211b controller and the reported gas/aerosol
  tolerances as context for the accuracy target. No code or pseudocode imported.
- Gustafsson, K. (1991). *Control Theoretic Techniques for Stepsize Selection in
  Explicit Runge-Kutta Methods*, ACM TOMS 17, 533–554.
  [Author's publication record](https://portal.research.lu.se/en/publications/control-theoretic-techniques-for-stepsize-selection-in-explicit-r/).
  Motivated a separate PI-controller experiment; our stated PI formula is a
  design choice, not a claim to reproduce that paper's algorithm.
- Read-only experiment inputs: `../reference/README.md`, `micm_cases.txt`,
  `compare_final.py`; `../scenarios/README.md`, `score.py`, `make_scenario.py`,
  scenario files, convergence logs and reference CSVs. These define gates,
  forcing semantics and all scored reference values. Reference C++ programs
  are built unchanged in this repository only.
