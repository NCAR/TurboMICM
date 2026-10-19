# R3. Extrapolation, exponential and SDIRK/ESDIRK methods versus Rosenbrock for stiff gas-phase chemistry at 2–3 digits

Literature review, 2026-10-18. Bracketed numbers refer to Part 2, which says for each source whether I read the full text, the abstract only, or only saw it cited. **[calc]** marks my own arithmetic or numerical check, not a source claim.

## Part 1. Findings

### Bottom line

1. The only head-to-head at atmospheric-chemistry accuracy is Sandu et al. 1997 [10][11]: box models of about 15–84 species, sparse LU, operator-split restarts. Sparse Rosenbrock codes (RODAS3, ROS3, RODAS) were fastest within about 4 digits. SEULEX (extrapolated linearly implicit Euler) trailed them on every problem at 1–2 digits, by more than 2× on the stratospheric problem and about 3× on the aqueous one, and led only beyond about 4 digits on the aqueous problem. SDIRK4 was comparable to RODAS and VODE but trailed RODAS at 2 digits. No exponential method was tested.
2. In the chemistry studies found, exponential integrators were slower than implicit methods at loose tolerance [13][17] and on Robertson [19]. They led only at tight tolerance, against a matrix-free baseline [17].
3. DIRK-type stage order is at most 2, which limits them to low-to-moderate tolerances [21]. One CFD study finds Rosenbrock equal or better than ESDIRK, and more robust [23].
4. Corrections to the brief: Curtis et al. [13] tested no Rosenbrock method (they used Radau IIA, exp4, exprb43 and CVODE). Imren and Haworth's abstract compares SEULEX only with BDF/CVODE [12]. Stone et al. [14] compare ROS4 only with explicit RKF45. Ehrig–Nowak–Deuflhard is a 1997 preprint, proceedings 1998 [8].

### 1. Extrapolation of the linearly implicit Euler method

**Method and control.** A substep solves (B − hA)(u₊ − u) = h f(u), with A the Jacobian at the start of the macro step or a good approximation [6]. Row j uses h_j = H/n_j and gives T_{j,1}. Aitken–Neville gives T_{j,k} = T_{j,k−1} + (T_{j,k−1} − T_{j−1,k−1})/(n_j/n_{j−k+1} − 1), and ε_j = ‖T_{j,j} − T_{j,j−1}‖ is the error estimate [6][8]. The error expansion has all powers of h, so T_{k,k} has order k [9]. Order and step size are chosen together by minimizing work per unit step [2] (review summary only). The Julia analogue uses a window of orders k−1, k, k+1 and a work count of f evaluations, Jacobian evaluations and substitutions [9]; its order limits are 3 to 12, start 5 [28].

**Sequences.** The harmonic sequence 1,2,3,… is used for this method [6][9]; Romberg (1,2,4,8,…) and Bulirsch (1,2,3,4,6,8,12,…) are the alternatives [9]. Harmonic is cheaper: k rows cost k(k+1)/2 substeps against 2^k − 1 for Romberg **[calc]**, and it load-balances best across threads [9]. I could not verify a "2,3,4,…" variant (Part 3).

**Smoothing step.** None for linearly implicit Euler. Gragg's smoothing belongs to the Bader–Deuflhard linearly implicit midpoint rule, whose error expansion has even powers [9]. Imren and Haworth report that the midpoint variant beats SEULEX in CFD [12].

**Stability.** The base entry R(z) = (1 − z/n_j)^(−n_j) is L-stable. Every T_{j,k} is a linear combination of such terms, so R(∞) = 0 for all entries **[calc]**. A-stability is not inherited. A ray scan of the diagonal T_kk **[calc]** found T₁₁ and T₂₂ A-stable. For k = 3…8 with the harmonic sequence it found A(α) with α ≈ 89.85°, 89.78°, 89.77°, 89.78°, 89.80°, 89.81°, with max |R(iy)| between 1.0016 and 1.013. The {2,3,4,…} and Bulirsch sequences give 89.8–89.9°, and Romberg is nearer 90°. So α stays within 0.25° of 90° and stiff decay is kept. On stiff problems the global-error expansion has perturbation terms, so order k may not be reached [4] (abstract only). I could not open Hairer–Wanner's own table [5].

**Cost per step** is in Table 1: one LU per row (I − h_jJ changes with j) plus n_j back-solves. A SEULEX-type controller also needs a row or two beyond the target order (my inference; the minimum order is 3 [28]). Utkarsh et al. count an LU as roughly one back-substitution for n < 200, a dense-BLAS statement [9].

**Reported chemistry performance.**

- *Atmospheric box models* [10]: atol 0.01 molecules cm⁻³, rtol from 1 downward, 40–60 min restarts. SEULEX never failed, yet it was "always significantly more expensive than RODAS3" at low accuracy and slower than the Rosenbrock codes on the 17-species problem. The authors note that "the virtue of extrapolation manifests itself most clearly when high accuracy is required."
- *ROBER, POLLU and others* [9]: at small tolerances (high accuracy) their extrapolation code was about 4× faster than Hairer's SEULEX on ROBER. Unthreaded it was about 2× ahead of Rodas4 and Radau on 3–8-ODE problems, and threaded about 2.5× ahead on the 20-ODE POLLU problem. They call extrapolation "generally inefficient" against Rosenbrock or SDIRK in serial benchmark suites. The wins come from within-method threading, which cell-parallel runs lack.
- *Combustion:* LIMEX was reported 10–15× faster than DASSL in 1-D flames, second-hand [7]. A semi-implicit Euler extrapolation code was modified and tested on combustion problems [3]. Imren–Haworth are in section 2.

### 2. Combustion comparisons

- **Imren and Haworth** [12]: SEULEX against CVODE in homogeneous reactors and engines, with mechanisms from under 50 to over 5,000 species and CPU time roughly linear in size. SEULEX beat CVODE in accuracy and efficiency in all cases. A conference abstract gives "at least a factor of two faster" for 100 or more species and credits operator-splitting restarts. Tolerances and any Rosenbrock comparison are unknown (body unread).
- **Curtis et al.** [13]: Radau IIA, exp4 and exprb43 against CVODE at atol 1e-10 and rtol 1e-6, on H₂ (13 species) and GRI-Mech 3.0 (53 species, 325 reactions). CVODE won on CPU and Radau IIA on GPU. They recommend Rosenbrock and W-methods next.
- **Stone et al.** [14]: ROS4 was 2.2–2.6× faster than RKF45 on GRI-Mech 3.0, with about 10× the RHS calls per step (finite-difference Jacobian) but about 25× fewer steps. Nothing on extrapolation, exponential or ESDIRK methods.

### 3. Exponential integrators

- **Methods** [15]: exprb32 (order 3, 2 stages, embedded exponential Rosenbrock–Euler) and exprb43 (order 4, 3 stages) are "fully explicit," need no linear solves and need a Jacobian every step. Only one Krylov approximation, of F(uₙ), is large. The tests were PDEs. exp4 was faster at moderate tolerance and exprb43 at stringent tolerance. Radau5 with direct solves took over 10× longer than exprb43 on a Schrödinger problem, a regime where direct solves are costly.
- **Chemistry:**
  - [13]: CVODE was 2.18× (Δt = 1e-6 s) and 10.1× (Δt = 1e-4 s) faster than exp4 on GRI-Mech. Observed orders were about 2 (exp4) and 3 (exprb43) with exact Krylov spaces, and 1 with approximate ones. The authors blame Arnoldi cost and thread divergence.
  - [17]: EPI3V against matrix-free CVODE (BDF, SPGMR) on GRI-Mech 3.0, n-butane (154 species, 680 reactions) and n-dodecane (105 species, 420 reactions). CVODE was modestly faster at loose tolerance, EPI3V led only below error ≈ 1e-5, and on n-dodecane EPI3V was about 10× slower (large real spread of the Jacobian spectrum).
  - [19] (preprint): on Robertson, explicit exponential schemes stall near 1e-4 absolute error while backward Euler, trapezoid and Radau reach 1e-8 to 1e-14.
- **EPIRK** [18][20]: EPIRK-W (order 3) and EPIRK-K (order 4) allow inexact Jacobians and were tested on PDEs. I read nothing chemistry-specific from Tokman or Loffeld–Tokman.
- **Cost reasoning [calc], assumed numbers:** a φ-vector product costs m sparse J·v products plus about 2nm² Gram–Schmidt flops. At n = 200, nnz(J) ≈ 2,000 and m = 20 that is ≈ 2.4×10⁵ flops, about 5 sparse LUs if one costs 5×10⁴ flops. The LU-free advantage disappears when sparse LU is cheap.

### 4. SDIRK and ESDIRK

- **Kennedy and Carpenter 2016** [21]: stage order cannot exceed 2, which "fundamentally limits their utility to settings requiring only low to moderate error tolerances." The recommended methods are stiffly accurate, L-stable and stage-order 2, with embedded and dense-output formulas. ESDIRK4(3)6L[2]SA is the suggested default "at moderate error tolerances" (Table 2).
- **Newton:** no iteration counts are given. The Jacobian is "generally not reevaluated at each stage but rather after many steps." The iteration tolerance is about 0.005–0.1 of the integration tolerance. First-stage derivative reuse and stage-value predictors are the main savings.
- **Order reduction:** little on Pareschi–Russo, moderate on Kaps, severe on van der Pol. Stage-order 2 plus L-stability beat stage-order 1 plus algebraic stability.
- **Versus Rosenbrock:** RODAS was fastest in most atmospheric cases at 1%, and SDIRK4 was second at 2 digits on CBM-IV [11]. Kennedy and Carpenter note that Hairer–Wanner found SDIRK4 "disappointing" while Sandu et al. cast SDIRK "in a more favorable light" [21]. In CFD with Jacobian-free Krylov solvers [23], RODASP and ESDIRK4 were best at fine tolerance, RODASP was more robust (ESDIRK4 failed some runs) and easier to implement, and SDIRK2 was competitive only at very coarse tolerance.

**Table 2. Methods from [21]** (s stages, s_I implicit; a_ii = γ):

| Method | (s, s_I) | order | γ | leading error A⁽ᵖ⁺¹⁾ |
|---|---|---|---|---|
| ESDIRK2(1)3L[2]SA (TR-BDF2 family) | (3,2) | 2 | 0.2929 | 0.0572 |
| ESDIRK3(2)4L[2]SA | (4,3) | 3 | ≈0.436 | 0.0366 |
| ESDIRK3(2)5L[2]SA | (5,4) | 3 | 0.2250 | 7.8e-4 |
| ESDIRK4(3)6L[2]SA (default) | (6,5) | 4 | 0.25 | 1.8e-3 |
| SDIRK4(1), Hairer–Wanner | (5,5) | 4 | 0.25 | 2.5e-3 |

### 5. Atmospheric chemistry, 2015–2026

I found no paper applying extrapolation or exponential integrators to gas-phase atmospheric chemistry. The only ESDIRK use found is TChem-atm [25], which offers the TINES TrBDF2 integrator (second-order trapezoid/BDF2, the ESDIRK2 type of [21]) beside CVODE. It chose rtol 1e-6 for 0.01% global accuracy and gives no TrBDF2-versus-CVODE timing. Rosenbrock methods remain the practice:

- KPP's integrators were "selected to be very efficient in the low to medium accuracy regime," relative error about 1e-2 to 1e-5 [27].
- Dreger et al. report typical rtol 1e-2 and atol 1–10 cm⁻³, using Ros3 and Rodas3 [24].
- A preprint reports ROS2 cutting E3SM chemistry cost by about 33% at a 180 s step against the default implicit solver, for 60 and 155 species [26].

### Table 1. Work per macro step, by family

N is the number of Newton iterations per implicit stage. [21] gives no value, so 2–3 is my assumption. Rosenbrock back-solves equal the stage count (one solve per stage) **[calc]**. Extrapolation counts assume f(y₀) is shared across rows **[calc]**.

| Method (order) | Jacobian | LU | back-solves | RHS evals | Source |
|---|---|---|---|---|---|
| ROS2 (2) | 1 | 1 | 2 | 2 | [27] |
| ROS3 (3) | 1 | 1 | 3 | 2 | [27] |
| RODAS3 (3) | 1 | 1 | 4 | 3 | [27] |
| ROS4 (4) | 1 | 1 | 4 | 3 | [27] |
| RODAS4 (4) | 1 | 1 | 6 | 5 | [27] |
| ESDIRK2(1)3L, TR-BDF2 type (2) | reusable | 1 | 2N | 2N | [21] |
| ESDIRK3(2)4L (3) | reusable | 1 | 3N | 3N | [21] |
| ESDIRK4(3)6L (4) | reusable | 1 | 5N | 5N | [21] |
| Extrapolation, harmonic, p = 3 | 1 | 3 | 6 | 4 | [6][9] |
| Extrapolation, harmonic, p = 4 | 1 | 4 | 10 | 7 | [6][9] |
| Extrapolation, harmonic, p = 5 | 1 | 5 | 15 | 11 | [6][9] |
| Extrapolation, 2,3,4,…, p = 4 | 1 | 4 | 14 | 11 | [6][9] |
| Extrapolation, Romberg, p = 4 | 1 | 4 | 15 | 12 | [6][9] |
| exprb32 (3) | 1 | 0 | 0 | 2 | [15]; 1 large + 1 small Krylov φ-product |
| exprb43 (4) | 1 | 0 | 0 | 3 | [13][15]; 3 Krylov φ-products |

At order 3, extrapolation needs about twice the back-solves and RHS evaluations of ROS3 and 3× the LUs **[calc]**. Rosenbrock methods need a fresh Jacobian every step [14]. Extrapolation evaluates one per macro step, and LIMEX allows an approximate or older one [6]. Modified-Newton ESDIRK can reuse one across many steps [21][14].

*Inference for the design:* with a cheap sparse LU and a 2–3 digit target, this literature favors third-order Rosenbrock or low-order ESDIRK/TR-BDF2 with Jacobian reuse. Extrapolation and exponential methods look like low-priority arms.

## Part 2. Sources

Citations are as I verified them. "Full" means I read the extracted text, all or the relevant sections.

1. P. Deuflhard, "Recent progress in extrapolation methods for ordinary differential equations," SIAM Review 27(4), 505–535 (1985), doi:10.1137/1027140. **Only saw it cited** (bibliographic records; the SIAM page returned 403). Sandu et al. [10] credit it with the idea of extrapolating the linearly implicit Euler method for stiff ODEs.
2. P. Deuflhard, "Order and stepsize control in extrapolation methods," Numer. Math. 41, 399–422 (1983), doi:10.1007/BF01418332. **Review summary only** (MaRDI). Order is chosen by minimizing work per unit step, using subdiagonal error estimates.
3. P. Deuflhard, E. Hairer, J. Zugck, "One-step and extrapolation methods for differential-algebraic systems," Numer. Math. 51(5), 501–516 (1987), doi:10.1007/BF01400352. **Abstract only** (UNIGE archive). It derives the perturbed order pattern of the semi-implicit Euler tableau for DAEs. The modified code was more efficient on combustion-chemistry problems.
4. E. Hairer, Ch. Lubich, "Extrapolation at stiff differential equations," Numer. Math. 52(4), 377–400 (1988), doi:10.1007/BF01462235. **Abstract only** (UNIGE record, seen in a search snippet). It gives global-error expansions on stiff problems for implicit Euler, linearly implicit Euler and the linearly implicit midpoint rule, with an explicit dominant perturbation term.
5. E. Hairer, G. Wanner, *Solving Ordinary Differential Equations II*, 2nd ed., Springer (1996), §IV.9–IV.10. **Not read.** Cited by [10] for SEULEX and by OpenFOAM's seulex documentation.
6. U. Maas, U. Nowak, "An adaptive method of lines for the simulation of complex laminar combustion processes," ZIB Preprint SC 96-56 (1996). **Full text of the time-extrapolation section.** It states the linearly implicit Euler step, the harmonic sequence for this method, the Aitken–Neville formula and the error estimate.
7. P. Deuflhard, J. Lang, U. Nowak, "Adaptive algorithms in dynamical process simulation," ZIB Preprint SC 95-16 (1995). **Read section 2.** LIMEX needs one iteration loop fewer than implicit codes. It quotes (from Maas and Warnatz) 10–15× over DASSL in 1-D and an estimated more than 100× in 2-D combustion, unchecked. The LIMEX reference is Deuflhard–Nowak, *Large Scale Scientific Computing*, Birkhäuser (1987), 37–50 (record only).
8. R. Ehrig, U. Nowak, P. Deuflhard, "Massively parallel linearly-implicit extrapolation algorithms as a powerful tool in process simulation," ZIB Report SC 97-43 (1997); Advances in Parallel Computing 12, Elsevier (1998), 517–524 (from a search summary of the ZIB record, unconfirmed). **Read the introduction.** It restates the tableau and ε_j and says cost is dominated by Jacobian evaluation and linear solves.
9. Utkarsh, C. Elrod, Y. Ma, K. Althaus, C. Rackauckas, "Parallelizing explicit and implicit extrapolation methods for ordinary differential equations," arXiv:2207.08135 (2022); the SciML documentation cites it as Elrod et al., HPEC 2022. **Full text.** Sequences, work model, midpoint smoothing, and benchmarks on ROBER, HIRES, POLLU and QSP. It says serial extrapolation is not state of the art in standard stiff benchmarks.
10. A. Sandu, J. G. Verwer, J. G. Blom, E. J. Spee, G. R. Carmichael, F. A. Potra, "Benchmarking stiff ODE solvers for atmospheric chemistry problems II: Rosenbrock solvers," Atmos. Environ. 31(20), 3459–3472 (1997). **Full text** (CWI scan). RODAS, RODAS3, ROS3, ROS4, SEULEX, VODE and TWOSTEP on seven box models.
11. A. Sandu, J. G. Verwer, M. van Loon, G. R. Carmichael, F. A. Potra, D. Dabdub, J. H. Seinfeld, "Benchmarking stiff ODE solvers for atmospheric chemistry problems I: implicit vs explicit," Atmos. Environ. 31(19), 3151–3166 (1997). **Full text** (CWI scan). Sparse SDIRK4, RODAS, VODE and LSODES plus explicit solvers. RODAS was fastest in the 1% region in most cases.
12. A. Imren, D. C. Haworth, "On the merits of extrapolation-based stiff ODE solvers for combustion CFD," Combust. Flame 174, 1–15 (2016), doi:10.1016/j.combustflame.2016.09.018. **Abstract only** (Penn State Pure). I also saw the abstracts of their ESSCI 2016 paper "Stiff ODE solvers for detailed chemistry in CFD" and their 2017 U.S. National Combustion Meeting paper "A midpoint rule based extrapolation solver for combustion CFD."
13. N. J. Curtis, K. E. Niemeyer, C.-J. Sung, "An investigation of GPU-based stiff chemical kinetics integration methods," Combust. Flame 179, 312–324 (2017), doi:10.1016/j.combustflame.2017.02.005; arXiv:1607.03884. **Full text.** Analytical Jacobians from pyJac. Exponential methods were slower than CVODE and Radau IIA on CPU and GPU.
14. C. P. Stone, A. T. Alferman, K. E. Niemeyer, "Accelerating finite-rate chemical kinetics with coprocessors: comparing vectorization methods on GPUs, MICs, and CPUs," Comput. Phys. Commun. 226, 18–29 (2018); arXiv:1608.05794. **Read the abstract, solver description and results.** It notes that Rosenbrock methods need J every step and implicit RK can reuse it.
15. M. Hochbruck, A. Ostermann, J. Schweitzer, "Exponential Rosenbrock-type methods," SIAM J. Numer. Anal. 47(1), 786–803 (2009), doi:10.1137/080717717. **Full text** (KIT repository). Coefficients of exprb32 and exprb43, the Krylov-oriented reformulation, and two PDE experiments against exp4, RKC, Krogstad and Radau5.
16. M. Hochbruck, A. Ostermann, "Exponential integrators," Acta Numerica 19, 209–286 (2010). **Skimmed** (KIT text). It calls exp4 "the first exponential integrator … competitive for certain time-dependent partial differential equations." I found no stiff-kinetics cost comparison.
17. J. Stewart, M. Tokman, F. Bisetti, V. Dallerit, O. Diaz-Ibarra (author names reordered from the arXiv listing "Stewart Jared, Tokman Mayya, …"), "Variable time-stepping exponential integrators for chemical reactors with analytical Jacobians," arXiv:2307.00136 (2023). **Read.** The abstract says EPI3V can "perform similar if not marginally better" than a BDF method on two mechanisms and worse on a third.
18. M. Narayanamurthi, P. Tranquilli, A. Sandu, M. Tokman, "EPIRK-W and EPIRK-K time discretization methods," J. Sci. Comput. (online 2018), doi:10.1007/s10915-018-0761-3; arXiv:1701.06528. **Abstract and test-problem list only** (shallow water, Allen–Cahn). Methods with inexact Jacobians. No chemistry.
19. C. Fronk, L. Petzold, "Performance evaluation of single-step explicit exponential integration methods on stiff ordinary differential equations," arXiv:2411.19374 (2024). **Abstract and Robertson section.** A preprint, not known to be peer-reviewed. Explicit exponential schemes lag implicit ones on Robertson, HIRES and van der Pol.
20. M. Tokman, J. Comput. Phys. 213(2), 748–776 (2006) (EPI); M. Tokman, J. Comput. Phys. 230, 8762–8778 (2011) (EPIRK); J. Loffeld, M. Tokman, J. Comput. Appl. Math. 241, 45–67 (2013); F. Bisetti, Combust. Theory Model. 16(3), 387–418 (2012). **Only saw them cited** (Acta Numerica, Curtis et al., search records). A Loffeld thesis abstract says exponential integrators "perform favorably" against Newton–Krylov implicit solvers, a PDE-style regime.
21. C. A. Kennedy, M. H. Carpenter, "Diagonally implicit Runge–Kutta methods for ordinary differential equations. A review," NASA/TM-2016-219173 (March 2016). **Full text** (NTRS). Source of Table 2, the stage-order statement, the default recommendation and the summary of earlier comparisons.
22. C. A. Kennedy, M. H. Carpenter, "Diagonally implicit Runge–Kutta methods for stiff ODEs," Appl. Numer. Math. 146, 221–244 (2019). **Only saw it cited** (PETSc manual). No abstract retrieved.
23. D. S. Blom, P. Birken, H. Bijl, F. Kessels, A. Meister, A. H. van Zuijlen, "A comparison of Rosenbrock and ESDIRK methods combined with iterative solvers for unsteady compressible flows," Adv. Comput. Math. 42(6), 1401–1426 (2016), doi:10.1007/s10444-016-9468-x. **Read the abstract, introduction and conclusions** (author manuscript). RODASP, ROS34PW2, ESDIRK3, ESDIRK4 and SDIRK2 in CFD. This is not chemistry.
24. R. Dreger, T. Kirfel, A. Pozzer, S. Rosanka, R. Sander, D. Taraborrelli, "Optimized step size control within the Rosenbrock solvers … in KPP version 2.2.3_rs4," Geosci. Model Dev. 18, 4273–4291 (2025), doi:10.5194/gmd-18-4273-2025. **Read most of the HTML.** Ros3 and Rodas3 with first- and second-order step controllers, typical tolerances, and 13–43% fewer function evaluations.
25. O. H. Díaz-Ibarra et al., "TChem-atm (v2.0.0): scalable performance-portable multiphase atmospheric chemistry," Geosci. Model Dev. 19, 1281–1299 (2026). **Read the first 100,000 characters.** TINES TrBDF2 and CVODE, CB05 and E3SM UCI mechanisms, no TrBDF2-versus-CVODE speed numbers.
26. R.-Y. Chien, J. S. Fu, "Using a two-stage Rosenbrock solver to improve surface ozone prediction and increase computational efficiency in … E3SMv1," EGUsphere preprint egusphere-2026-1451 (28 Apr 2026). **Abstract and introduction.** Not peer-reviewed.
27. KPP documentation, "Numerical methods" (v2.6.0), kpp.readthedocs.io. **Read.** Rosenbrock stages and function calls, the SDIRK, Radau and SEULEX list, and the accuracy-regime statement.
28. OrdinaryDiffEq.jl documentation, "Extrapolation" (docs.sciml.ai). **Read.** The implicit-Euler extrapolation solver is "similar to Hairer's SEULEX." Sequences are harmonic (default), Romberg or Bulirsch, and order limits are 3 to 12 with start 5. Nothing on stability.

## Part 3. Open questions and leads

1. **Hairer–Wanner IV.9–IV.10 not read.** I could not verify SEULEX's default sequence, work function, minimum order or A(α) table. My recollection that SEULEX offers sequences starting at 2 (such as 2,3,4,6,8,…) is unverified. Hairer's Fortran header would answer it, but I fetched no code; a library copy of the book or Deuflhard–Bornemann (2002) §4 is the alternative. The IV.10 work–precision diagrams for POLLU (also in the Bari/CWI test-set reports) would give SEULEX against RODAS and RADAU5 at 1e-2 to 1e-4.
2. **Imren and Haworth full text.** I need the tolerances behind "outperforms," the BDF settings, and whether Rosenbrock or Radau was ever included.
3. **Not accessed:**
   - Stone and Bisetti, AIAA 2014-0822 (solver comparison on combustion chemistry).
   - Bisetti 2012; Tokman 2006 and 2011; Loffeld–Tokman 2013; Kennedy–Carpenter 2019.
   - Constantinescu–Sandu 2010 (SIAM J. Sci. Comput. 31, 4452) and ANL TM-306, which study Euler-based extrapolation and linear stability (TM-306 has an appendix on linearly implicit Euler). Both returned 403. They might confirm my numerical A(α) angles.
4. **Unverified citation:** for the brief's "Ehrig, Nowak, Deuflhard 1999" I found only the 1997 preprint, proceedings 1998 [8]. A 1999 chapter I vaguely recall could not be confirmed, so I did not cite it.
5. **Nothing found for 2015–2026** on extrapolation, exponential or ESDIRK integrators compared against Rosenbrock on gas-phase atmospheric chemistry. Such work may sit inside model-development papers or theses I did not reach.
6. **Untested in any source I found:** a Rosenbrock-W or other frozen-Jacobian scheme on chemistry. Curtis et al. suggest W-methods, and Jacobian reuse is the one cost lever that ESDIRK and W-methods have and plain Rosenbrock lacks.
7. **Possible order-reduction remedy:** a 2025 NASA report on singly implicit RK methods with stage order 3–4 (NTRS 20250008379) was seen only in a search listing and not read.
