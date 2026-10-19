# R4: Step-size control, error norms, Jacobian/LU reuse, restarts, positivity

*Literature review for the Tier 3 integrator design (about 200 species, 550 reactions, 5–30 min operator-split steps restarted every step, cheap sparse LU, 2–3 digits). Written 2026-10-18. No code was downloaded or quoted. [n] is the source number in Part 2. Read status: **FT** full text, **PT** relevant parts, **AB** abstract or record only, **CIT** only seen cited. The fetch tool cannot decode PDFs, so I saved them and read them with pdftotext or as page images.*

## Part 1: Findings

### 1. Step-size controllers

**1.1 Elementary controller (Rosenbrock/KPP form).** The step is accepted if the embedded estimate has ‖est‖ ≤ τ, and the next step is h·(τ/‖est‖)^(1/(p+1)), applied after accepted and rejected steps alike. Shampine: failed steps are expensive, so be conservative, and large changes cannot be justified by the asymptotics [1]. Rosenbrock practice adds a safety factor and limits:
`h_new = h·min(qmax, max(qmin, δ·Err^(−1/(p̂+1))))`, with p̂ the order of the embedded estimator (exponent 1/3 for the 3(2) pairs RODAS3/ROS3), Err a scaled RMS norm (§3), accept if Err ≤ 1 [8],[7]. Values: δ=0.9, qmax=10, qmin=0.1 [8]; KPP‑2.x defaults δ=0.9, qmax=6, qmin=0.2, h_start=1e‑5 s [7]. [7] attributes these to Hairer–Nørsett–Wanner and Hairer–Wanner [32], which I could not read.

**1.2 Gustafsson's predictive controller.** I could not read the 1994 paper [11]. [10] restates it: assume ‖Est_n‖=C_n·h_n^(q+1) with C_{n+1}/C_n ≈ C_n/C_{n−1}. The rule then uses the two previous steps and, "in combination with the standard one", makes RADAU5 fail fewer steps on some stiff problems (H&W II pp. 124–127, as cited; unverified). From that assumption I derived (not copied) `h_{n+1} = h_n·(h_n/h_{n−1})·(τ/r_n)^(1/k)·(r_{n−1}/r_n)^(1/k)`, k=q+1, which is (β1,β2,α)=(2,−1,−1) in the filter notation below. KPP offers it (ICNTRL(11)=0) for RADAU5/Runge–Kutta [19]; its Rosenbrock codes use the elementary rule [7]. RADAU5 defaults per [10]: growth ≤8, safety 0.9; for their own post-rejection rule they give typical values 0.2 and 0.9 (r̃=max{0.2·r, 0.9·r̂}). I found nothing on a predictive option in RODAS.

**1.3 Söderlind's digital filters.** From his slides [2], with c_n=tol/r_n (r_n the error estimate), ρ_n=h_n/h_{n−1}, k=estimator order+1:
`ρ_{n+1} = c_n^(β1/k)·c_{n−1}^(β2/k)·ρ_n^(−α),  h_{n+1} = ρ_{n+1}·h_n`

| controller | β1 | β2 | α |
|---|---|---|---|
| elementary (I, deadbeat) | 1 | 0 | 0 |
| PI3040 (= PI.3.4, my mapping) | 7/10 | −4/10 | 0 |
| PI3333 | 2/3 | −1/3 | 0 |
| PI4020 (= PI.4.2, my mapping) | 3/5 | −1/5 | 0 |
| H211PI | 1/6 | 1/6 | 0 |
| H211b | 1/b | 1/b | 1/b |

H211b is stable for b ≥ 1 (poles at 0 and 1−2/b). b=2 is the deadbeat FIR filter, and larger b suppresses error noise more but responds more slowly [2]. A search snippet of Arévalo & Söderlind (2017, unread) gives b=4 as default, b∈[3,6] (unverified). I found no readable source for H312b coefficients (my recollection, β=(1,2,1)/b, α=(3,1)/b, b=8, is unverified). The framework papers [3],[4] were available only as abstracts.

*Evidence.*
- **Chemistry [7].** Dreger et al. replaced the KPP Rosenbrock controller by H211b (b=1, k=1.7; Ros3/Rodas3, rtol 1e‑2, 1e‑3 for aerosol). EMAC function evaluations fell 43 % (gas), 27 % (cloud) and 13 % (aerosol), and a one-year run used 11.4 % fewer core-hours. The retuned first-order controller (qmax=100 for gas, δ=1.5 for cloud, h_start=1e‑4 for aerosol) gave only 4 %. In the box model, retuned δ≈1.3 and H211b gave similar savings (24–32 %). Larger b "always" cost more evaluations; b=1 (edge of the stable range) caused no problems but "should be kept in mind".
- **Overestimated error [7].** The embedded estimate overestimates the true local error on nearly all accepted steps. At rtol 1e‑2 the median accuracy was 2.8–4.1 digits and the worst species 1.6–2.3 (target 2), so δ≈1.4 pays. Rejections caused by the overestimate occurred in <5 % of 250 inspected steps.
- **Computational stability [2],[5].** DASSL with filter control gave much better computational stability (error ∝ tol) at no extra cost. The slides show the stiff Chemakzo problem and claim log10(C/c)≈0.05 is possible.
- **Rejections [9].** Ranocha et al. (explicit RK, CFD) accept when the limited factor ≥0.81 and cite Söderlind for fewer unnecessary rejections. They use the smooth limiter κ(a)=1+arctan(a−1) and give no rejection counts. I found no rejection counts for stiff kinetics in any source.

**1.4 Rejection policy and growth limits.**
- KPP Rosenbrock: after a rejection the next step may not grow (limit 1.0, not 10), and a rejected first step gives h/10 [8]. After ≥2 consecutive rejections the step is multiplied by 0.1 [7].
- CVODE (BDF): after a failed test (h'/h)^(q+1)‖Δ‖=ε/6. η is capped at 0.2 after two failures, the order or step is reset after three, and it gives up after seven. Growth is capped at 10 (10⁴ on the first step), and changes <1.5× are skipped [18].
- The H211b pseudocode in [7] applies only the ×0.1 factor after ≥2 rejections (no qmin/qmax).

### 2. Initial step size and restarts

**2.1 Initial step.** Hairer–Nørsett–Wanner Vol. I (Sec. II.4 by my recollection; p.169 per [9]) [32]; not read. The algorithm as I know it, with constants confirmed only through a code-derived description seen in a search summary (unverified): d0=‖y0‖, d1=‖f0‖ (scaled norms); h0=0.01·d0/d1 (1e‑6 if either <1e‑5); one Euler step; d2=‖f(t0+h0,y1)−f0‖/h0; h1=(0.01/max(d1,d2))^(1/(p+1)) (or max(1e‑6, 1e‑3·h0) if both tiny); h=min(100·h0, h1). [12] is abstract level only; [34] is only cited.

**2.2 Restarts in chemistry–transport practice.**
- [13]: "any splitting step involves a start up". One-step Rosenbrock methods "can rapidly change step size", important "in view of the many operator-split restarts". Chemistry is often put in the middle of Strang splitting because each integration needs a start-up; splitting discontinuities create stiff transients that cost extra steps.
- [14]: RODAS3 "can require a rather small step size in the initial phase". In operator splitting these transients are insignificant, so they advocate an a priori chosen start step far larger than the smallest time constant, and use ROS2 with fixed 10–20 min steps (variable steps cause load imbalance and impede vectorization).
- [8]: every split interval restarts at a user h_start. For 1 h restarts h_min=0.1 s and h_start=60 s. Later 60 s was found "too large", and the 15 min restarts used h_start=h_min=0.1 s. With 4× as many restarts every solver became 3–4× more expensive, so cost scaled with the number of restarts: "all solvers spend most of their time in the start-up phase". Too-large h_start/h_min with rtol ≥0.1 made solutions drift. Advice: small h_min/h_start and rtol ≤0.01, since growth of up to 10× per step recovers quickly.
- [7]: at every model step (10–20 min) integrators need initial substeps of 1e‑5 s or shorter, and the first model step of a run needs far more substeps than later ones. h_start=1e‑5 was best in the box model (≥1e‑4 crashed, ≤1e‑6 slightly worse); in EMAC, changing it gave no gain for MECCA or SCAV and −3 % for AERCHEM at 1e‑4. A dynamic initial-step option exists but was not evaluated.
- **Carrying the step over:** KPP documents it: "For multiple restarts, use Hnew as Hstart in the subsequent run" (RSTATUS(3)=Hnew, the last predicted, not yet taken step) [19].
  - I found no published measurement of its benefit.
  - I found none on carrying filter memory across restarts. [7] and [9] initialise the memory neutrally (errOld=facOld=1; ε₀=ε₋₁=1).
  - A 2009 workshop slide starts each call at h=0.1·(t_end−t_start) [17]. This is secondary and differs from KPP's 1e‑5.
  - I found no Sandu or Sportisse paper on restarts and step control.

### 3. Error norms and tolerances

- **Norm.** Scaled RMS in every source I could read: Err=√((1/n)Σ(est_j/Tol_j)²), with Tol_j=atol+rtol·|y_new,j| [8] or rtol·max(|y_old|,|y_new|) [7],[9]. CVODE uses a weighted RMS with W_i=1/(rtol·|y_i|+atol_i) [18]. I found no RMS-versus-max study. Shampine says only "a suitable norm" [1]; Söderlind's 2005 report discusses the choice (abstract only) [6].
- **Typical tolerances (molecules cm⁻³).**

| source | rtol | atol |
|---|---|---|
| benchmark [8] | 1…1e‑4 | 0.01 |
| GEOS‑Chem [15] | "normally" ≤1e‑3 | 1e4·rtol in tests |
| EMAC MECCA/SCAV [7],[16] | 1e‑2 | 10 (1 for key short-lived radicals) |
| EMAC AERCHEM [7] | 1e‑3 | 1 |
| "typical" per [7] citing [33] | up to 1e‑2 | up to 100 |

 I found no source on mixing-ratio solvers.
- **Loose control is dangerous.** rtol ≥0.1 gave breakdowns or >100 % error [8]; "forcing RODAS3 to compute in lower accuracy is sometimes problematical" [13].
- **Yardstick.** SDA=−log10(max_k ER_k) [8]. ER_k is the RMS over time of the relative error of species k, restricted to times with y_k ≥ a, a=1e6 molec cm⁻³ (1e4 stratosphere; a=1 gave nearly the same conclusions). The reference is RADAU5 at rtol 1e‑12. [15] uses a=1e6 too.
- **Tolerance proportionality** [1]. Controlling an order‑p estimate while advancing at order p* makes the error scale as tol^((p*+1)/(p+1)); an internal tolerance μ=τ^((p+1)/(p*+1)) fixes it.
- **Many cells.**
  - [17] (workshop slides of the SC09 work): vectorized KPP Rosenbrock uses E=max(E₁…E_vn) over lanes and one shared h.
  - [16]: one cell per GPU thread; step counts differ between cells, giving thread divergence.
  - [14]: variable steps hamper vectorization.
  - [7]: H211b raised work in a few day–night-transition boxes but caused no net imbalance.

### 4. Jacobian and LU reuse

- **BDF heuristics** (CVODE docs, "inherited from VODE" [18]; the VODE and SUNDIALS papers [22],[23] were not read).
  - The Newton matrix M=I−γJ is refreshed at start, after >20 steps, when |γ/γ̄−1|>0.3, and after a non-fatal convergence failure or an error-test failure.
  - J itself is re-evaluated only at start, after >50 steps, after a convergence failure with an outdated matrix when |γ/γ̄−1|<0.2, or after a failure forcing a step cut.
  - KPP RADAU5/RK: the step is held and the LU reused if Qmin<Hnew/Hold<Qmax, and J is not recomputed if the Newton rate is <ThetaMin (defaults not stated) [19].
- **Rosenbrock.** One J, one LU and s solves per step, "most of the CPU time" [8].
  - Stale J breaks the order conditions.
  - [8] lists, without exploring them: time-lagged exact J with fixed h (one LU for several steps, order kept by changing scalar coefficients, Verwer–Scholz); a simpler J, which loses conservation and needs many more order conditions (W-methods, Steihaug–Wolfbrandt 1979); species partitioning.
  - [20] (third-order method, 72 runs, prescribed steps): accuracy falls when J is frozen (strongly for two problems), and gains appear only where J+LU dominate the cost. Combining it with step-size and error control needed more understanding of the nonlinearity (left as future work in 1980).
  - ROS2 (γ=1+1/√2) is second-order consistent for any matrix A, so it is a W-method [13],[14]. It was used for approximate factorization, not time-lagged reuse.
  - I found no later chemistry paper on frozen-Jacobian Rosenbrock/W-methods.
- **Speedups.** No direct measurement found. Indirectly, Rosenbrock with a fresh J each step beat the BDF code VODE at ≤4 digits [8] and was competitive with LSODE and SPRINT [13], and KPP Rosenbrock was about 2× faster than SMVGEAR II at moderate accuracy [15], though confounded by vectorization. SMVGEAR's own LU/Jacobian policy: not found [21].

### 5. Positivity

- **Why it matters.** Negative values can destabilise chemistry (c'=−λc²-type: the Jacobian sign flips), "even resulting in a numerical blow up within a few time steps" [13].
- **Order barrier.** For a class of stiff linear systems, no RK or multistep method of order >1 is positive for all steps (Bolley–Crouzeix 1978, per [13],[25]). "Truly positive Rosenbrock methods for general non-linear chemical kinetics systems do not exist". RODAS3 is not positive. With local error control, lack of positivity "will automatically lead to step rejections", whereas fixed steps may blow up [13].
- **Clipping (negatives to 0).**
  - ROS2 clips both stages. This breaks mass conservation but gave "no notable loss in accuracy", since clipping is rare for γ=1+1/√2; it "does not guarantee stability"; in fixed-step tests clipping usually enlarged the stable step (one mechanism got worse) [14].
  - [25]: clipping is "a poor choice, which cannot guarantee robust and conservative algorithms", and adaptivity reduces negatives only empirically.
  - KPP's switch: leave, zero, warn or stop (ICNTRL(16)) [19].
- **Projection [26] (AB).** After a mass-preserving step, any negative component is projected onto the reaction simplex (a quadratic optimization). "More accurate than clipping and allow larger time steps" for systems unstable outside the positive orthant; best with positivity-favouring schemes (Rosenbrock variants in [31], unread).
- **Bertolazzi [27] (AB).** Implicit schemes: first order unconditionally positive, higher orders apparently for large steps; they conserve linear invariants and use M-matrix inversions instead of Jacobians.
- **Patankar-type.** MPE/MPRK2 are unconditionally positive and conservative [28]; order theory and third-order variants are in [29]. [24] shows issues even on linear 2×2 problems:
  - oscillations about the steady state unless Δt is below a CFL-like bound (e.g. Δt≤2 for MPRK(2,2,1) on the normalised test);
  - order reduction (down to 1) when a component is near 0;
  - zeros must be replaced by ε>0 (1e‑180 on Robertson, 1e‑35 on HIRES);
  - semi-implicit variants blow up on Robertson;
  - the behaviour carries to nonlinear stiff problems.

  [25] tested a ~200-equation mass-action system, where MPRK was the most efficient of the second-order unconditionally positive, conservative schemes.
- **Production.** Clipping in LOTOS-type models [14]. KPP-based models (GEOS‑Chem, EMAC/MECCA) have selectable zeroing [19], and MECCA has a "posdef" Rosenbrock variant (name only [7]). CMAQ, CAM‑chem and GEOS‑Chem practice: not verified.

### 6. Implications (my inference, not from the sources)
- **Controller.** In a restarted one-step Rosenbrock the first step of each interval is effectively elementary control, so memory-based filters (H211b, PI) pay only after a few steps. A carried-over h_new or a cheap start-up estimate is the lever.
- **Tolerance.** Estimates are conservative (§1.3), so tolerances or δ can be tuned to hit 2–3 digits at less work; verify against a tight reference.
- **Positivity.** Clipping plus rejection of bad steps has the best literature support at low order; Patankar-type schemes lose accuracy near zero concentrations.

## Part 2: Sources

1. Shampine, L.F. (2005). Error estimation and control for ODEs. *J. Sci. Comput.* 25(1–2):3–16. **FT** (Feb-2004 preprint, tlakoba.w3.uvm.edu/math337/paper_Shampine_error.pdf). EPS vs EPUS, step formula (8), local extrapolation, RADAU5's order‑5 method with a filtered order‑3 error estimate, internal tolerance.
2. Söderlind, G. (2018). *Adaptive Time-Stepping, Part I: How Control Theory and Digital Filters Enhance IVP Solver Performance.* Lecture slides, Budapest, V2.0 (math.bme.hu/…/Gustaf1.pdf). **FT.** The author's filter recursion, coefficient table, H211b stability and b-dependence, and Chemakzo/DASSL stability plots. Primary source for §1.3 formulas.
3. Söderlind, G. (2003). Digital filters in adaptive time-stepping. *ACM TOMS* 29(1):1–26. DOI 10.1145/641876.641877. **AB.** Step control as a digital filter, giving smoother step sequences. The one full-text copy I found (nacad.ufrj.br) refused connections.
4. Söderlind, G. (2002). Automatic control and adaptive time-stepping. *Numer. Algorithms* 31:281–310. DOI 10.1023/A:1021160023092. **AB.** Linear-feedback step control; coordination with nonlinear solvers in implicit methods.
5. Söderlind, G., Wang, L. (2006). Adaptive time-stepping and computational stability. *J. Comput. Appl. Math.* 185:225–243. DOI 10.1016/j.cam.2005.03.008. **AB.** Standard DASSL against a digital-filter version: small changes give much better computational stability at no extra cost.
6. Söderlind, G. (2006). Time-step selection algorithms: adaptivity, control, and signal processing. *Appl. Numer. Math.* 56(3–4):488–502; and (2005) *A specification of a control system for adaptive time stepping*, Lund technical report. **AB** both. The report covers choosing the error norm.
7. Dreger, R., Kirfel, T., Pozzer, A., Rosanka, S., Sander, R., Taraborrelli, D. (2025). Optimized step size control within the Rosenbrock solvers for stiff chemical ODE systems in KPP 2.2.3_rs4. *Geosci. Model Dev.* 18:4273–4291. DOI 10.5194/gmd-18-4273-2025. **FT.** The most relevant source: KPP defaults, overestimated error, tuned first-order control vs H211b (CAABA, EMAC, 735-species MOM), tolerances, h_start tests.
8. Sandu, A., Verwer, J.G., Blom, J.G., Spee, E.J., Carmichael, G.R., Potra, F.A. (1997). Benchmarking stiff ODE solvers for atmospheric chemistry problems II: Rosenbrock solvers. *Atmos. Environ.* 31(20):3459–3472. **PT** (OCR text and page images of the CWI-hosted published version). Err formula, step rule, restart set-up, tolerances, SDA metric, §3.2 on approximate Jacobians.
9. Ranocha, H., Winters, A.R., Castro, H.G., Dalcin, L., Schlottke-Lakemper, M., Gassner, G.J., Parsani, M. (2023). On error-based step size control for discontinuous Galerkin methods for compressible fluid dynamics. arXiv:2209.07037v2. **PT** (§2.3, Algorithm 2). Explicit RK; PID form, κ limiter, 0.81 acceptance rule, neutral memory at start. Journal version not checked.
10. Hernández-Abreu, D., González-Pinto, S., Montijano, J.I. (2020). Variable step-size control based on two-steps for Radau IIA methods. Manuscript submitted to ACM (venue unconfirmed; riull.ull.es/xmlui/handle/915/39058). **PT.** Restates Gustafsson's assumption and RADAU5's step rules.
11. Gustafsson, K. (1994). Control-theoretic techniques for stepsize selection in implicit Runge–Kutta methods. *ACM TOMS* 20(4):496–517. **CIT** (via [10],[1]; not read).
12. Gladwell, I., Shampine, L.F., Brankin, R.W. (1987). Automatic selection of the initial step size for an ODE solver. *J. Comput. Appl. Math.* 18:175–192. **AB** (third-party index).
13. Verwer, J.G., Hundsdorfer, W., Blom, J.G. (2002). Numerical time integration for air pollution models. *Surv. Math. Ind.* 10:107–174. **PT**, of the 1998 CWI version (MAS‑R9825, ir.cwi.nl/pub/4620/04620D.pdf); the final PDF is a scan without text. Restarts, Rosenbrock/ROS2/RODAS3, positivity theory, W-methods.
14. Verwer, J.G., Spee, E.J., Blom, J.G., Hundsdorfer, W. (1999). A second-order Rosenbrock method applied to photochemical dispersion problems. *SIAM J. Sci. Comput.* 20(4):1456–1480. **PT** (§§1–3). Start-up argument, fixed prescribed steps, ROS2 positivity and clipping tests.
15. Eller, P., Singh, K., Sandu, A., Bowman, K., Henze, D.K., Lee, M. (2009). Implementation and evaluation of an array of chemical solvers in GEOS-Chem. *Geosci. Model Dev.* 2:89–96. **FT.** GEOS‑Chem tolerances, SDA metric, KPP Rosenbrock vs SMVGEAR II.
16. Alvanos, M., Christoudias, T. (2017). GPU-accelerated atmospheric chemical kinetics in the ECHAM/MESSy (EMAC) Earth system model (version 2.52). *Geosci. Model Dev.* 10:3679–3693. **PT.** One cell per GPU thread, thread divergence, EMAC tolerances.
17. Linford, J.C., Michalakes, J., Vachharajani, M., Sandu, A. (2009). Multi-core acceleration of chemical kinetics for simulation and prediction. *Proc. SC09.* **AB** (UCAR record) plus Linford's IMAGe TOY 2009 slides (**PT**) with vectorized Rosenbrock pseudocode (E=max over lanes). Paper not obtained.
18. SUNDIALS CVODE documentation v7.8.0, "Mathematical Considerations" (sundials.readthedocs.io/en/v7.8.0/cvode/Mathematics_link.html). **PT**, verified against the raw page. Matrix and Jacobian refresh rules, error test, step limits, WRMS norm.
19. KPP documentation (kpp.readthedocs.io/en/stable): "Integrator options (ICNTRL and RCNTRL)" and "Integrator status". **PT.** ICNTRL(16), ICNTRL(11), Qmin/Qmax, Hnew-as-Hstart advice; no FacMin/FacMax/Hstart defaults there.
20. Verwer, J.G., Scholz, S. (1983). Rosenbrock methods and time-lagged Jacobian matrices. *Beitr. Numer. Math.* 11:173–183 (read as CWI preprint NW 82/80, 1980, ir.cwi.nl/pub/8986/8986D.pdf). **PT** (abstract, introduction, conclusions).
21. Jacobson, M.Z., Turco, R.P. (1994). SMVGEAR: a sparse-matrix, vectorized Gear code for atmospheric models. *Atmos. Environ.* 28A:273–284. **CIT** (the Stanford copy was a 2-page scan I could not read).
22. Brown, P.N., Byrne, G.D., Hindmarsh, A.C. (1989). VODE: a variable-coefficient ODE solver. *SIAM J. Sci. Stat. Comput.* 10:1038–1051. **CIT.**
23. Hindmarsh, A.C., et al. (2005). SUNDIALS: suite of nonlinear and differential/algebraic equation solvers. *ACM TOMS* 31(3):363–396. **CIT**; [18] used instead.
24. Torlo, D., Öffner, P., Ranocha, H. (2022). Issues with positivity-preserving Patankar-type schemes. *Appl. Numer. Math.* 182:117–147. arXiv:2108.07347. **PT** (abstract, introduction, §7–8).
25. Formaggia, L., Scotti, A. (2011). Positivity and conservation properties of some integration schemes for mass action kinetics. *SIAM J. Numer. Anal.* 49(3):1267–1288 (read as MOX report 16/2010). **PT** (abstract, introduction, conclusions).
26. Sandu, A. (2001). Positive numerical integration methods for chemical kinetic systems. *J. Comput. Phys.* 170(2):589–602. DOI 10.1006/jcph.2001.6750. **AB** (MTU record).
27. Bertolazzi, E. (1996). Positive and conservative schemes for mass action kinetics. *Comput. Math. Appl.* 32(6):29–43. **AB** (secondary summary only).
28. Burchard, H., Deleersnijder, E., Meister, A. (2003). A high-order conservative Patankar-type discretisation for stiff systems of production–destruction equations. *Appl. Numer. Math.* 47(1):1–30. **AB.**
29. Kopecz, S., Meister, A. (2018). On order conditions for modified Patankar–Runge–Kutta schemes. *Appl. Numer. Math.* 123:159–179. arXiv:1702.04589. **AB.** I could not identify a 2019 paper; the BIT 58:691–728 (2018) third-order paper is known only from [24]'s references.
30. Shampine, L.F., Thompson, S., Kierzenka, J.A., Byrne, G.D. (2005). Non-negative solutions of ODEs. *Appl. Math. Comput.* 170:556–569. **CIT** (via a MathWorks reference list; not read).
31. Sandu, A. (2002). Time-stepping methods that favor positivity for atmospheric chemistry modeling. In *Atmospheric Modeling*, IMA Vol. Math. Appl. 130, Springer, pp. 21–37. **CIT** (via [25]).
32. Hairer, E., Nørsett, S.P., Wanner, G. (1993), *Solving ODEs I*; Hairer, E., Wanner, G. (1996), *Solving ODEs II*, Springer. **CIT** (via [7],[9],[10]); not read.
33. Zhang, H., Linford, J.C., Sandu, A., Sander, R. (2011). Chemical mechanism solvers in air quality models. *Atmosphere* 2(3):510–532. **CIT** (via [7]; details from memory, verify).
34. Watts, H.A. (1983). Starting step size for an ODE solver. *J. Comput. Appl. Math.* 9:177–191. **CIT** (via [7]; details from memory, verify).

## Part 3: Open questions and leads

1. **H312b and recommended filters.** No readable H312b coefficients, and no stiff-specific recommendation beyond [2] and [7]. Leads: the full text of [3] and [6]; the KPP‑3 integrator `rosenbrock_h211b_qssa` (KPP docs, branch feature/h211b) as a reference implementation.
2. **H&W text unread.** Verify II.4 (initial-step constants), IV.8 (Gustafsson formula, RADAU5 defaults, any RODAS predictive option) and the rejection factors; read [11], [12] and [34]. The constants I quote are secondary or recalled.
3. **Restart evidence gap.** No published benchmark of carrying h_new or controller memory across operator-split restarts. Open: what GEOS‑Chem, CAM‑chem and CMAQ actually pass as the initial step.
4. **Max versus RMS over cells.** Only slides [17] and statements in [14],[16]. Leads: the SC09 full text, [33] and the norm discussion in [6].
5. **Units.** No source on mixing-ratio tolerances or species-dependent atol beyond [16].
6. **SMVGEAR [21].** Its LU/Jacobian reuse and block-wise step control are unverified.
7. **Frozen or time-lagged Jacobian Rosenbrock/W-methods with error control** for chemistry: no post-1980 paper found. A hit I did not read: sim.mathematik.uni-halle.de/numdiff14/files/abstracts/jaxt.pdf.
8. **Positivity in production codes.** Not verified for GEOS‑Chem, CMAQ or CAM‑chem. Unread leads: [30] (MATLAB's NonNegative option), [31], the full text of [26]. No MPRK study on a 200-species photochemical mechanism found ([25] has a ~200-equation case, summary parts read only).
9. **Minor.** GEOS‑Chem "Forced exit from Rosenbrock, step size too small" GitHub issues (titles only) hint at step-collapse failures; [13]'s final version not checked against the 1998 text.
