# Accuracy policy for atmospheric chemistry

The recommended operating point targets at least **three significant digits overall and two digits for every scored species** on every validation scenario: at most 0.1% RMS relative error overall and 1% RMS for each species. A tighter point is provided for sensitivity studies; two cheaper points serve applications that accept the level-two target. These are engineering accuracy criteria for the supplied chemistry operator, not universal tolerances for every atmospheric application.

Sandu et al. (1997) motivate assessing atmospheric stiff solvers at modest accuracy with work–precision comparisons, rather than judging them by nominal local tolerances. Verwer et al. (1999) highlight the special difficulty of short chemistry intervals and restart transients. More recently, Dreger et al. (2025) show that controller changes need long-run chemical error checks as well as runtime measurements. These motivate the scenario-wide, species-resolved validation here. Sources and exactly what was used are listed in `REFERENCES.md`.

One percent is useful as a cheaper operating regime, but it is too generous as the sole scientific default. Bates et al. (2024), for example, report chemical-mechanism changes with sub-percent effects on important atmospheric quantities. A solver error budget should leave room to distinguish such effects. Choosing the level-three criterion, and a point with margin above it, is this study's judgment; it is not a threshold prescribed by those papers.

## What is measured

The unmodified supplied `scenarios/score.py` measures `-log10(sqrt(mean(relative_error^2)))`. It includes every cell and snapshot and every species whose reference is at least `1e-14` mol m-3. Each species receives its own RMS score, so a good overall average cannot conceal an arbitrarily poor species. All species and reactions still participate in the ODE and error control, including those below the reporting threshold.

The scoring threshold is about 6,022 molecules cm-3. The new uniform absolute tolerance, `1e-18` mol m-3, is about 0.60 molecules cm-3 and is four orders below that threshold. This limits absolute-error dilution without specifying favored species. It is a local error scale, not a guarantee of four digits globally. Adaptive histories, stiff transients, repeated forcing changes, and small-negative projection can all change the accumulated error. That is why operating points are selected from measured reference agreement.

The reported references are the supplied tight six-stage MICM solutions, independently checked for roughly ten digits of convergence. The reference uncertainty is far below the targets here. No method changes the piecewise-constant-in-each-step photolysis convention: comparison isolates integration error for exactly that prescribed forcing, not error from representing sunlight by these outer steps.

## Full disclosed validation

These are the lowest observed scores across all six complete 1,000-cell scenarios and every supplied snapshot, from `scans/operating-points-clang.csv`. Repetitions have identical outputs and operation counts; timing varies.

| Point | Relative tolerance | Minimum overall digits | Overall RMS bound | Minimum species digits | Species RMS bound |
| --- | ---: | ---: | ---: | ---: | ---: |
| fast | `1e-2` | 2.041 | 0.91% | 1.123 | 7.54% |
| balanced | `5e-3` | 2.353 | 0.45% | 1.445 | 3.59% |
| accurate | `5e-4` | 3.313 | 0.049% | 2.409 | 0.40% |
| recommended | `3e-4` | 3.599 | 0.026% | 2.624 | 0.24% |
| precise | `1e-4` | 4.289 | 0.0052% | 3.125 | 0.075% |

Bounds are rounded upward. The recommended point provides margin over level three, and the precise point is useful when trace-species sensitivities justify additional cost. The fast point only narrowly clears the disclosed level-two target; the balanced point supplies more margin. Hidden scenarios may change which point qualifies. The highest errors often involve ISOP or bromine chemistry, which is one reason to retain a per-species requirement and multiple operating points.

The inherited default three-stage method is a compatibility baseline, not an accuracy oracle. At its large absolute floor it performs well on the short Tier 1 benchmark but fails or accumulates very large errors on several longer-step scenarios. Even a finite, nonnegative final state can be inaccurate. Its tighter-tolerance baseline succeeds and is included in the work–precision comparisons.

## Nonnegative output and conservation

Every new trial must be finite and pass the error test. A concentration more negative than its allowed local scale forces a smaller retry. Smaller accepted negative concentrations are projected to zero. `largest_clip` records the largest single correction; it is **not** a cumulative conservation-error bound. The largest observed corrections in the disclosed runs are `2.409e-16`, `2.409e-16`, `3.768e-18`, `3.237e-18`, and `1.165e-18` mol m-3 for the five points in table order. At the recommended point, the largest individual correction is about 0.033% of the reporting threshold.

This policy guarantees nonnegative accepted endpoints when a run succeeds, but does not make the underlying Rosenbrock formulas positivity-preserving or guarantee exact elemental conservation after projection. The independent nonlinear small-case test verifies its stoichiometrically weighted invariant to accumulated roundoff, and the stiff transfer test checks its equilibrium within integration error. Full scenario reference agreement and refinement are the relevant checks of projection's effect on TS1. The official gate is checked on all saved values; the native runner also checks every outer endpoint.

## Limits of the conclusion

RMS criteria do not bound the error of every individual cell/species value. Concentrations below the reporting threshold still need scrutiny if an application uses a very sensitive diagnostic. The local tolerances should be tightened for such work; absolute and relative tolerances remain runtime parameters.

Additional 48-hour tests vary latitude, season, seed, and pollution, with independently converged MICM references; their results are recorded separately in `scans/`. The large 10,000-cell extensions use the supplied references for the first 1,000 cells and check finite/nonnegative output on all 10,000. This extends numerical and throughput evidence without claiming reference-relative accuracy for unreferenced cells.

The recommendation is supported for this full TS1 box-chemistry operator, the disclosed forcing intervals, and the reported additional tests. A coupled atmospheric model should validate its own transport/chemistry splitting, emissions, deposition, diagnostics, and longer duration. Nothing here substitutes a reduced mechanism, a surrogate, or agreement with a loose default solver for that validation.
