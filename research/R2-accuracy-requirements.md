# R2 — How accurate must a chemistry solver be, and how do published solver comparisons measure and conclude?

Literature review, 2026-10-18. [Brackets] are source numbers from Part 2. "≈" marks values I read off a figure or computed myself; "(unverified)" marks anything I could not confirm in a primary text. Quotations are verbatim.

---

## Part 1 — Findings

### 1.1 Short answers

1. **The target accuracy is modest, about 1 % (SDA = 2).** "A numerical comparison should focus on modest accuracies, say relative accuracies near 1%. Higher accuracy levels are redundant for the actual practice of air pollution modelling." [1, §3; near-identical statement in 3, §3.1]. "…the greatest interest lies in high efficiency for very low accuracy (two figures at most)" [2]. Zhang et al. are stricter: overall numerical accuracy "on the order of 1%", numerical errors "one or two orders of magnitude below the overall target", hence a relative tolerance of 0.1 % ("about 3 accurate digits") [7].
2. **Operational settings give about 1–3 digits** (§1.6), at or below ranking level L = 3. Some schemes have no integration-error estimate: EMEP's TWOSTEP uses prescribed sub-steps and a fixed iteration count [16]; UKCA uses backward Euler with only a Newton convergence test [15].
3. **Time-step and splitting error dwarfs the solver tolerance.** Splitting error in GEOS-Chem is bounded by ≈ 10 % [22]. Operator durations of 10–60 min give relative errors of 5–35 % for NOx and secondary inorganic aerosol (SIA) against a finer-duration "truth" (O3 and CO < 2 % below 30 min), versus model–observation errors of ~30 % (NOx), 20–40 % (SIA), ~20 % (O3) [21].
4. **Reduced-chemistry papers accept 1–2 % median error** against the operational solver's own full-mechanism run, with up to 3 % for OH and 6–10 % for NO2 in some regions [11–14].
5. **SDA conventions differ**, and the concentration threshold ranges from 1 to 10⁶ molec cm⁻³ (§1.2).
6. **Rosenbrock methods win at 1–3 digits**; RODAS (4th order) and SEULEX take over above ≈ 3–4 digits [2]. In GEOS-Chem, KPP Rodas3 needs ≈ half the time of SMVGEAR II at equal SDA [5, 6].

### 1.2 The SDA definition (Sandu et al.) and its variants

From Part I eqs. 15–17 and Part II eq. 6, checked against page images [1, 2]. With y the reference, ŷ the numerical solution, and N+1 = 121 hourly samples over 5 days:

```
ER_k    = sqrt( (1/|T_k|) · Σ_{n∈T_k} ((y_k(t_n) − ŷ_k(t_n)) / y_k(t_n))² ),   T_k = { n : |y_k(t_n)| ≥ a }
SDA_avg = −log10( (1/m) Σ_k ER_k )     (Part I eq. 16)
SDA_max = −log10( max_k ER_k )         (Part I eq. 17; the only form in Part II)
```

It is the time-RMS of the *pointwise* relative error, aggregated over species. "SDA = 2 means 1% accuracy in the error measure used" [1, 2]. Threshold *a* drops individual time samples below *a*; ER_k is dropped only if no sample passes. Purpose: "to avoid chemically meaningless large relative errors for concentration values smaller than 1.0 mlc cm⁻³" [1].

| Source | Norm | Threshold a (molec cm⁻³) | Over species | Reference |
|---|---|---|---|---|
| Sandu I [1] | pointwise rel. RMS over time | 1 | mean and max | RADAU5, rtol 1e-12, atol_i = 1e-15·η_i |
| Sandu II [2] | same | 1e6 troposphere, 1e4 stratosphere; a = 1 gave "nearly the same conclusions" | max | RADAU5, rtol 1e-12 |
| Verwer 1996 [3] | ratio of sums √(Σ(ŷ−y)²/Σy²), 56 points | none printed | mean | "sufficiently accurate" |
| Verwer 1999 [4] | as [1] | 1e-4 × time-mean of species k (scale-free) | mean | accurate solution |
| Eller [5], Henze [6] | as [1], over grid cells | 1e6 | max | RTOL 1e-8; ATOL 1e-3 [5], 1e2 [6] |
| Zhang [7] | norm-wise RMS over species, final time | 100 | one norm | Radau-5A, RelTol 1e-10, AbsTol 1e2 |
| Santillana 2010, Shen 2020 [11, 12] | as [1] | 1e6 | median | full-mechanism run |
| Shen 2022, Lin 2023 [13, 14] | as [1] | boxes holding 99 % of species mass (+ C > 10 in [14]) | median / mean | full-mechanism run |

My conversions (1 molec cm⁻³ = 1.661e-18 mol m⁻³): 1e4 → 1.7e-14 and 1e6 → 1.7e-12 mol m⁻³. The project's 1e-14 mol m⁻³ is ≈ 6.0e3 molec cm⁻³: about 170 times below the 1e6 used with GEOS-Chem, close to the stratospheric 1e4 of [2].

### 1.3 Sandu et al. 1997, Parts I and II [1, 2]

**Problems.** Seven box problems, 5 days. Splitting is mimicked by restarting every 3600 s (120 restarts), i.e. a 1800 s transport step with Strang splitting. Photolysis follows time of day (C⁰ in A; C¹, zero at night, elsewhere). Emissions are applied at the start of each split interval, except in F.

| Problem | Mechanism | Species / reactions |
|---|---|---|
| A (I) | EUSMOG | 15 / 15 |
| A (II) | TMk CH4/CO/HOx/NOx, 40-min restart | 18 / 36 |
| B, C | CBM-IV urban, rural | 32 / 70 thermal + 11 photolytic |
| D, E | "AL" (STEM-II) urban, rural | 84 / 142 + 36 photolytic |
| F | NASA stratosphere, no emissions | 34 / 84 + 25 photolytic |
| G | gas–aqueous | 65 / 77 + 11 gas, 39 aqueous, 39 mass-transfer |

Part II adds 15-min restarts for B and D (7.5-min transport step). Rate coefficients are piecewise constant per interval, because time-continuous exponentials "can be as expensive as the sparse matrix factorization" [2].

**Tolerances.** Part I: per solver, "not specified here". Part II: atol = 0.01 molec cm⁻³; rtol ∈ {1, 3e-1, 1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4}. rtol ≥ 0.1 makes Rosenbrock solvers "drift away from the real solution"; the fix is "rtol ≤ 0.01" [2].

**Headline comparisons** (the work–precision figures are log-scale and I did not digitize them; these are the text's statements):
- Part I: sparse implicit solvers are "often more efficient than the dedicated explicit ones, particularly when approximately two or more figures of accuracy are required". Sparse RODAS is "most competitive in the 1% error region" (VODE wins on A and F). TWOSTEP is the best explicit solver and performs "very efficiently" at accuracies coarser than 1 % for tropospheric gas-phase problems [1].
- Aqueous G, estimated 5-day cost of explicit solvers: QSSA ≈ 350 d, extrapolated QSSA ≈ 650 d, CHEMEQ ≈ 15 d, TWOSTEP-Jacobi ≈ 270 d, TWOSTEP-Seidel ≈ 10.3 h. Implicit solvers need < 10 s [1, Table 10].
- Part II: "within 4 digits of accuracy the Rosenbrock solvers clearly provide the most cost-effective solutions." RODAS3 and ROS3 are best to ≈ 3 digits on B–E, then RODAS. On F, VODE and SEULEX are "more than 2 times slower" at 1 %. On G, SEULEX is ≈ 3× slower at 2 digits but best above 4. With 15-min restarts every run costs 3–4× more [2].

### 1.4 Verwer et al. 1996 and 1999 [3, 4]

- **1996.** TWOSTEP (BDF2 with Gauss–Seidel) against VODE (dense black-box; analytic Jacobian; plus sparse LU). Problems: EUSMOG (15 species), CIRK methane (18–19 species, 46 reactions), EMEP (66 species, ~140 reactions); 112 h, 56 two-hourly restarts; rtol = 10⁻ˡ (l = 1…5), atol = 1.0 molec cm⁻³. Sparsity made VODE ≈ 6× faster than the dense version at low accuracy on the 66-species problem. TWOSTEP "outperforms VODE3 convincingly" on the small problems and is comparable on the large one for "the accuracy range of greatest practical interest". For ~20-species models, sparsity "does hardly pay". An erratum recomputed Fig. 3, shifting SDA by 0.6–1.0.
- **1999 (ROS2).** "The accuracy demand for atmospheric dispersion problems is modest so that a low order method is a natural choice." ROS2 (2nd order, 2 stages) with fixed 10- and 20-min steps was "excellent" on CBM-IV and RIVM; TWOSTEP with 2 iterations "generates significant errors" at 20 min on CBM-IV; ROS2 costs ≈ 2× TWOSTEP. In a 1D test, Strang splitting had "no adverse effect on accuracy whatsoever for step sizes smaller than or equal to 20 minutes". Local error control "creates load imbalancing in parallel implementations and interferes with vectorization".

### 1.5 GEOS-Chem lineage and reduced-chemistry papers

- **Henze 2007, App. A [6].** Global one-day chemistry-only runs, 1-h external step, RTOL 1e-1…1e-5, ATOL 1e6…1e-2 molec cm⁻³. "The Rosenbrock solver is nearly twice as efficient as the SMVGEARII solver." "…we choose to run the Rosenbrock solver at tolerance levels that yield an SDA of ∼1.0 as the standard setting for this work." A week-long run changed by < 10 % for most species and by 10–15 % for HNO2, HNO4, IAP, INO2, ISOP, N2O5, NO, NO2, PP, RIP.
- **Eller 2009 [5].** 106 species, 311 reactions, 43,858 cells, 7-day chemistry-only run; RTOL 1e-1…1e-3, ATOL = 1e4·RTOL molec cm⁻³; reference RTOL 1e-8, ATOL 1e-3.
  - Fig. 4, Ox (≈ ±15 %): Rodas3 reaches ≈ 2 digits in ≈ 300 s and ≈ 2.6 digits at RTOL 1e-3 in ≈ 500 s; SMVGEAR II reaches only ≈ 2.0 digits, at ≈ 620 s; Runge–Kutta needs ≈ 1000 s and SDIRK ≈ 1400 s for 2 digits. RTOL (1e-1, 3e-2, 1e-2, 3e-3, 1e-3) → SDA ≈ (1.1, 1.4, 1.9, 2.3, 2.6) for Rodas3.
  - Table 2: Rodas3 needs 9 function calls and 3 steps per advection step per cell; SMVGEAR II needs 39 and 28.
  - Table 1: the tight-tolerance "references" of Rodas3 and SMVGEAR II agree to only 2.1 (NOx), 3.3 (Ox), 0.5 (PAN), 3.7 (CO) digits, while Rodas3 and Rodas4 agree to 7.6–10.4. My inference: SMVGEAR II does not converge to the true solution at its tightest setting; the paper gives no cause.
  - "Since we normally use tolerances of 10⁻³ or lower for GEOS-Chem". Sparse LSODE gave "3–4 accurate digits" in ≈ 300 s but was dropped because it "does not reach the asymptotic regime" at the loose tolerances.
- **Zhang 2011 [7].** MECCA marine mechanism (223 species, 560 reactions) against a Radau-5A reference, RelTol 1e-6…1e-1. Rosenbrock methods "have become popular due to their efficiency at moderate accuracy requirements" (text only; plot not read). Absolute tolerances must be vectors varying with species and altitude.
- **Jacobson 1998 [8].** SMVGEAR II relative tolerance "usually set to 10⁻³"; absolute tolerance fixed at 1e3 (urban) or 1e5 (global) molec cm⁻³, then made adaptive (up to 1e7). Above 1e7, normalized gross errors (mean |N−E|/E, cutoffs 1e3 and 1e6, reference RTOL = ATOL = 1e-9) exceed 1 %; they stayed at 0.10–0.92 %, and time fell 14–67 %.
- **Reduced chemistry [11–14]** (all against the full-mechanism run):

| Source | Error | Speedup |
|---|---|---|
| Santillana 2010 | one year: NOx ≤ 5 %, CO/OH/O3 < 1 %, "same magnitude" as LSODES vs SMVGEAR II | 50 % |
| Shen 2020 | median 0.2–0.8 %; NO2 1–2 % (up to 10 % at high latitudes), "well within other systematic sources of errors" | 32–41 % |
| Shen 2022 | median < 0.4 % boundary layer, 0.8–2 % free troposphere, 2.4–3.2 % stratosphere | 44–50 % |
| Lin 2023 | mean < 1 % key species; OH up to 3 %, NO2 6–10 % | 32 % |

### 1.6 What operational models run

| Model | Method | Tolerance / step | Source |
|---|---|---|---|
| GEOS-Chem native (≤ v11) | SMVGEAR II (variable-order BDF) | "tolerances of 10⁻³ or lower" [5]; SMVGEAR II default RTOL 1e-3, adaptive ATOL 1e5–1e7 molec cm⁻³ (global) [8]; GEOS-Chem's own ATOL unverified | [5, 8] |
| GEOS-Chem v11-02 | KPP Rosenbrock | ATOL 1e-2; RTOL 0.5e-2 (UCX), 2e-1 (non-UCX), a "temporary fix" for "Step size too small" | [10] |
| GEOS-Chem 14.6.2 docs | KPP Rosenbrock | ATOL 1e-2, RTOL 0.5e-3, tightened Jan 2023 "to avoid integrate errors by halogen chemistry" | [10] |
| GEOS-Chem 14.7/14.8 docs | `rosenbrock_autoreduce`, RODAS3.1 (Hg: Rodas3); 20-min step | per-species KPP_AbsTol / KPP_RelTol; defaults 0.01 and 0.05; BrCl example 1e-3 / 5e-3; units not stated (molec cm⁻³ inferred from [5]) | [10, 12, 14] |
| CMAQ 4.4 (2004) | EBI (mechanism-specific, per-species tolerances), ROS3, SMVGEAR | ROS3: RTOL 1e-3, ATOL 1e-7 ppm (≈ 2.5e6 molec cm⁻³ at the surface); SMVGEAR: 1e-3, 1e-9 ppm; "For most model applications, the EBI solvers should give the best performance." Current defaults unverified | [9] |
| WRF-Chem 4.4 | KPP, "only Rosenbrock type" | RTOL/ATOL in `wkpp_constants.F`, no values in the guide; namelist example chemdt 1.5 min | [20] |
| CAM-chem / MOZART | generated explicit and "semi-implicit" solvers | no tolerances in [18]; 1800 s step [19]; a 2nd-order Rosenbrock cut surface-O3 bias by up to 5.2 ppb, refining the first-order step to 180 s did not | [18, 19] |
| UKCA | backward Euler + Newton–Raphson | relative change 1e-4, ≤ 50 iterations then halve the step; 1-h chemistry step, 20-min transport | [15] |
| EMEP MSC-W | TWOSTEP | no tolerance: 5 × 20 s then 7 × ~160 s per 20 min; 3 / 2 / 1 Gauss–Seidel iterations by layer | [16] |
| CHIMERE | TWOSTEP (BDF2 + Gauss–Seidel) | 1 or 2 iterations; 2 "strongly recommended" at twice the cost; step control unverified | [17] |
| EMAC/MECCA, LOTOS-EUROS | not verified (Part 3) | | [23] |

GEOS-Chem's RTOL moved over a factor ≈ 400 (5e-4 to 2e-1) across versions, and the docs tie changes to robustness: loosening can rescue a first step that will not converge [10]. On Eller's Rodas3 map (§1.5), RTOL 1e-3 gives ≈ 2.6 digits; 5e-4 would give ≈ 3 (my extrapolation).

### 1.7 How authors justify "good enough"

- **Error taxonomy** [7]: data errors (meteorology, rate coefficients, emissions, initial and boundary conditions), model errors, numerical errors ("error due to operator splitting is important, though difficult to quantify"). Rate-coefficient uncertainty is named, not quantified.
- **Splitting and step** [21, 22]: ≈ 10 % bound. GEOS-Chem tolerances were set to "keep simulation times reasonable while maintaining acceptable numerical accuracy" [22]. Going from 2°×2.5° to 4°×5° changes error by an order of magnitude [21].
- **Comparison with other errors:** NO2 error "within other systematic sources of errors" [12]; adaptive-solver differences "of the same magnitude" as between two accepted solvers [11]; quasi-Newton vs Newton differences ~1e-7 for O3, below the 1e-4 solver tolerance for OH [15].

### 1.8 Reading across to the project (my interpretation)

- L = 2 and L = 3 match the literature's 1 % target [1, 3] and the 0.1 % of [7]. "L − 1 digits per species" is stricter than the literature's max-over-species [2, 5, 6] or mean [1, 3].
- The project's species threshold (≈ 6e3 molec cm⁻³) is far below the 1e6 used with GEOS-Chem, so radicals and trace species those papers exclude will count.
- The project uses tight-tolerance references like [1, 2, 5–7]; the reduced-chemistry papers compare with an operational run.

---

## Part 2 — Sources

Read status: **full** = main text read; **partial** = relevant sections only; **abstract** = abstract only; **cited** = only seen cited.

1. Sandu, A., Verwer, J. G., van Loon, M., Carmichael, G. R., Potra, F. A., Dabdub, D., Seinfeld, J. H. (1997). Benchmarking stiff ODE solvers for atmospheric chemistry problems I: implicit vs explicit. *Atmos. Environ.* 31(19), 3151–3166 (DOI not verified). **Full** (CWI scan, equations checked on page images). Defines SDA; seven problems, nine solvers, focus on 1 %. RODAS fastest near 1 %; explicit solvers fail on the aqueous problem.
2. Sandu, A., Verwer, J. G., Blom, J. G., Spee, E. J., Carmichael, G. R., Potra, F. A. (1997). Benchmarking stiff ODE solvers for atmospheric chemistry problems II: Rosenbrock solvers. *Atmos. Environ.* 31(20), 3459–3472. doi:10.1016/S1352-2310(97)83212-8. **Full** (scan, page image checked). Introduces RODAS3, ROS3; rtol range, thresholds, 15-min restarts. Rosenbrock most cost-effective to 4 digits.
3. Verwer, J. G., Blom, J. G., van Loon, M., Spee, E. J. (1996). A comparison of stiff ODE solvers for atmospheric chemistry problems. *Atmos. Environ.* 30(1), 49–58, erratum 30(9), 1531. Printed PII 1352-2310(95)00283-9 (a repository record shows a different check digit). **Full** (scan). TWOSTEP against sparse VODE; rtol/atol as in §1.4; "relative accuracies near 1%".
4. Verwer, J. G., Spee, E. J., Blom, J. G., Hundsdorfer, W. (1999). A second-order Rosenbrock method applied to photochemical dispersion problems. *SIAM J. Sci. Comput.* 20(4), 1456–1480. **Full** (CWI report version, ir.cwi.nl/pub/4728, probably MAS-R9717). Fixed 10–20 min steps; scale-free threshold; splitting error negligible up to 20 min in 1D.
5. Eller, P., Singh, K., Sandu, A., Bowman, K., Henze, D. K., Lee, M. (2009). Implementation and evaluation of an array of chemical solvers in the Global Chemical Transport Model GEOS-Chem. *Geosci. Model Dev.* 2, 89–96. doi:10.5194/gmd-2-89-2009. **Full.** KPP Rodas3/Rodas4/Runge–Kutta/SDIRK against SMVGEAR II: tolerances, SDA, cross-solver agreement, CPU.
6. Henze, D. K., Hakami, A., Seinfeld, J. H. (2007). Development of the adjoint of GEOS-Chem. *Atmos. Chem. Phys.* 7, 2413–2433. doi:10.5194/acp-7-2413-2007. **Partial** (Appendix A). Rosenbrock vs SMVGEAR II; SDA ≈ 1 standard setting.
7. Zhang, H., Linford, J. C., Sandu, A., Sander, R. (2011). Chemical mechanism solvers in air quality models. *Atmosphere* 2(3), 510–532. doi:10.3390/atmos2030510. **Full** (accuracy and results sections). Error taxonomy, 0.1 % target, MECCA benchmark, vector absolute tolerances.
8. Jacobson, M. Z. (1998). Technical note: Improvement of SMVGEAR II on vector and scalar machines through absolute error tolerance control. *Atmos. Environ.* 32(4), 791–796. **Full** (scan, OCR and page images). RTOL 1e-3, adaptive ATOL, NGE < 1 %.
9. U.S. EPA/CMAS (2004). CMAQ CHEM_SOLVER_NOTES, 30 Sep 2004 (CMAQ 4.4). https://cmascenter.org/cmaq/documentation/4.4/CHEM_SOLVER_NOTES.txt. **Full.** Solver ranking, default tolerances, EBI's inadequacy for integrated reaction rates.
10. GEOS-Chem Classic documentation ("Update chemical mechanisms with KPP", "Understand what error messages mean": 14.6.2 and latest, geos-chem.readthedocs.io) and the v11-02 wiki page (wiki.seas.harvard.edu). **Partial** (via a fetch tool that returns extracts). Tolerance history, per-species mechanism, integrators.
11. Santillana, M., Le Sager, P., Jacob, D. J., Brenner, M. P. (2010). An adaptive reduction algorithm for efficient chemical calculations in global atmospheric chemistry models. *Atmos. Environ.* 44(35), 4426–4431. doi:10.1016/j.atmosenv.2010.07.044. **Full.** Errors 1–5 % against LSODES; same size as LSODES vs SMVGEAR II.
12. Shen, L., Jacob, D. J., Santillana, M., Wang, X., Chen, W. (2020). An adaptive method for speeding up the numerical integration of chemical mechanisms in atmospheric chemistry models: application to GEOS-Chem version 12.0.0. *Geosci. Model Dev.* 13, 2475–2486. doi:10.5194/gmd-13-2475-2020. **Partial.** RRMS definition; 20-min step; no RTOL/ATOL stated.
13. Shen, L., Jacob, D. J., Santillana, M., Bates, K., Zhuang, J., Chen, W. (2022). A machine-learning-guided adaptive algorithm to reduce the computational cost of integrating kinetics in global atmospheric chemistry models: application to GEOS-Chem versions 12.0.0 and 12.9.1. *Geosci. Model Dev.* 15, 1677–1687. doi:10.5194/gmd-15-1677-2022. **Partial.** Uses a "fourth-order Rosenbrock" method; errors by atmospheric domain.
14. Lin, H., Long, M. S., Sander, R., Sandu, A., Yantosca, R. M., Estrada, L. A., Shen, L., Jacob, D. J. (2023). An adaptive auto-reduction solver for speeding up integration of chemical kinetics in atmospheric chemistry models: implementation and evaluation in the Kinetic Pre-Processor (KPP) version 3.0.0. *J. Adv. Model. Earth Syst.* 15, e2022MS003293. doi:10.1029/2022MS003293. **Partial.** 291 species, 913 reactions, 32 % faster; tolerances left to users, no values.
15. Esentürk, E., Abraham, N. L., Archer-Nicholls, S., Mitsakou, C., Griffiths, P., Archibald, A., Pyle, J. (2018). Quasi-Newton methods for atmospheric chemistry simulations: implementation in UKCA UM vn10.8. *Geosci. Model Dev.* 11, 3089–3108. doi:10.5194/gmd-11-3089-2018. **Partial.** UKCA solver, 1e-4 convergence test, ~1e-7 differences for O3.
16. Simpson, D., et al. (2012). The EMEP MSC-W chemical transport model – technical description. *Atmos. Chem. Phys.* 12, 7825–7865. doi:10.5194/acp-12-7825-2012, with Supplement §S2.3. **Partial.** TWOSTEP sub-steps and iteration counts.
17. Menut, L., et al. (2013). CHIMERE 2013: a model for regional atmospheric composition modelling. *Geosci. Model Dev.* 6, 981–1028. doi:10.5194/gmd-6-981-2013. **Partial** (solver paragraph; PDF text partly garbled). TWOSTEP adaptation, iteration advice.
18. Lamarque, J.-F., et al. (2012). CAM-chem: description and evaluation of interactive atmospheric chemistry in the Community Earth System Model. *Geosci. Model Dev.* 5, 369–411. doi:10.5194/gmd-5-369-2012. **Partial** (one paragraph). Explicit and "semi-implicit" solvers per species; no tolerances. A CESM forum reply points to Kinnison et al. 2007 (*J. Geophys. Res.* 112, D20302), Auxiliary Material §A1.1; not retrievable, and Kinnison's main text has no solver description.
19. Sun, J., Fu, J. S., Drake, J., Lamarque, J.-F., Tilmes, S., Vitt, F. (2017). Improvement of the prediction of surface ozone concentration over conterminous U.S. by a computationally efficient second-order Rosenbrock solver in CAM4-Chem. *J. Adv. Model. Earth Syst.* 9, 482–500. doi:10.1002/2016MS000863. **Abstract only.** ROS-2 vs first-order implicit at 1800 s.
20. WRF-Chem Version 4.4 User's Guide, chapters 4 and 6. https://etrp.wmo.int/pluginfile.php/87070/mod_resource/content/1/Users_guide.WRF-Chemv4.4.pdf. **Partial.** KPP solvers; where tolerances live.
21. Philip, S., Martin, R. V., Keller, C. A. (2016). Sensitivity of chemistry-transport model simulations to the duration of chemical and transport operators: a case study with GEOS-Chem v10-01. *Geosci. Model Dev.* 9, 1683–1695. doi:10.5194/gmd-9-1683-2016. **Partial.** Error vs operator duration, compared with model–observation errors.
22. Santillana, M., Zhang, L., Yantosca, R. (2016). Estimating numerical errors due to operator splitting in global atmospheric chemistry models: transport and chemistry. *J. Comput. Phys.* 305, 372–386. doi:10.1016/j.jcp.2015.10.052 (arXiv:1505.02835). **Partial.** ≈ 10 % upper bound; GEOS-Chem tolerances "internally set".
23. Sander, R., et al. (2019). The community atmospheric chemistry box model CAABA/MECCA-4.0. *Geosci. Model Dev.* 12, 1365–1385. doi:10.5194/gmd-12-1365-2019. **Partial.** States no integrator or tolerance; its only tolerance is the skeletal-reduction criterion (20 % above 1 pmol mol⁻¹ for HCHO, HO2, NO, O3, OH).
24. Hundsdorfer, W., Verwer, J. G. (2003). *Numerical Solution of Time-Dependent Advection-Diffusion-Reaction Equations.* Springer. **Cited only** (in [11, 22]); not accessed.
25. Brasseur, G. P., Jacob, D. J. (2017). *Modeling of Atmospheric Chemistry.* Cambridge University Press. **Cited only** (in [12–14]); chapter contents not verified.

---

## Part 3 — Open questions and leads

**Citation corrections to the brief.**
- "Lin et al. 2023, Fast integration of chemical kinetics … KPP 3.0.0, *GMD*": not found. The KPP 3.0.0 paper is [14] (JAMES).
- "Shen et al. 2020 (ACP 20, machine-learning-guided…)": the machine-learning-guided paper is Shen et al. 2022, *GMD* 15 [13]; the 2020 paper is *GMD* 13 [12].
- "Shen et al. 2022 … adaptive auto-reduction solver": that title is Lin et al. 2023 [14].
- Santillana et al. 2010: authors Santillana, Le Sager, Jacob, Brenner; pp. 4426–4431 [11].
- Sandu 1997 Part I and Part II have different authors (van Loon, Dabdub, Seinfeld only in I; Blom, Spee only in II).

**Unresolved**
1. *GEOS-Chem*: which species override the documented defaults (ATOL 0.01, RTOL 0.05), and the units. The RTOL used in [12–14] is not stated. Lead: the species database in the GEOS-Chem repository.
2. *CAM-chem/MOZART*: Newton/backward-Euler convergence thresholds. Leads: Kinnison 2007 Auxiliary Material §A1.1, and the preprocessor's implicit-solver module (a search summary mentioned 1e-3 and 1e-4 thresholds; unverified).
3. *WRF-Chem*: RTOL/ATOL values in `wkpp_constants.F` [20]; not read.
4. *CMAQ*: current (5.x) default solver and tolerances; only the 2004 values [9] are verified.
5. *EMAC/MECCA*: default integrator and tolerances; [23] does not state them.
6. *LOTOS-EUROS*: TWOSTEP use not verified; Manders et al. 2017 (*GMD* 10, 4145) has no solver description in the text I searched.
7. *Kinetic-data uncertainty vs solver error*: I found no paper that sets a solver tolerance against rate-coefficient uncertainty; [7] lists it as a model error without a number. Leads (not searched): JPL kinetics evaluations, Monte-Carlo rate-constant studies.
8. *SDA-versus-CPU at 2 and 3 digits per solver* sits in the log-scale figures of [1, 2], which I did not digitize.
9. *Hundsdorfer & Verwer 2003, Brasseur & Jacob 2017, Jacobson 2005, Seinfeld & Pandis*: not accessed. I could not identify a stiff-ODE appendix in Seinfeld & Pandis (unverified whether one exists).
10. *SMVGEAR II's 0.5–3.7-digit disagreement* with Rodas3 at tight tolerance [5]: solver structure or configuration?
