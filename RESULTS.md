# Results (in progress)

The focused mass-action, generated Robertson mechanism, and configuration reader
checks, analytic decay, and a generated Robertson integration pass with clang
and gcc-15 at the default optimization settings. They do
not constitute the requested MICM case suite or a numerical solver comparison.
No benchmark result is available yet.

Both tutorial rate configurations match independently hand-assembled reaction
and rate tables over ten integration steps. The non-user-defined tutorial also
matches the copied reference program at its printed precision. This is not a
full-precision comparison or a substitute for the requested case suite.

The TS1 translation contains 210 species (209 state variables) and 547 reactions.
Its initial values match the reference defaults CSV exactly, and generated rates,
forcing, and Jacobian are finite at those defaults. Numerical TS1 integration has
not yet been validated against the reference.

A clang `-O3 -DNDEBUG` build on the supplied 100-cell, 100-step, 30-second
input passed the final-state gate: largest relative difference 4.2e-14. With
10 threads, one preliminary run took 15,406 ms; the reference built with the
same flags took 219 ms. These are single runs, not ranking measurements. The
large gap reflects the current dense, independent-cell implementation.

The same 100-cell gate also passed with gcc-15; largest relative difference
was 3.59e-14 and one preliminary timed run took 19,356 ms.

After adding sparse elimination and automatic species reordering, single runs
of the 100-cell check took 424 ms with clang and 543 ms with gcc-15, versus
219 ms for the clang reference. Both passed the final-state comparison (largest
relative differences 4.18e-14 and 6.04e-14, respectively).

At the requested scale (10,000 cells, seed 2026, 100 timed steps of 30 s,
10 threads), one clang run took 42,002 ms; the clang reference with the same
optimization flags took 27,426 ms. All 2,110,000 output values passed the gate;
the largest relative difference was 1.77e-13. These are single measurements,
not the specified median-of-five ranking measurements.

The case reporter currently has 10 passing ported cases and reports the other
670 cases as unported failures on both compilers. Focused checks include the
backward Euler convergence thresholds, exact implicit decay, failed-step
reductions, and dense/sparse agreement; they are not the full case suite.

The latest case report has 30 passing rate and backward-Euler cases and 650
unported failures with each compiler. New reversible, Taylor-series, JPL-19,
branching, falloff, and surface calculations also have direct reference-value
checks; a separate fixture exercises generation of the added rate families.

The latest case report has 32 passing cases and 648 unported failures with both
compilers. All five Rosenbrock methods pass the 999-step analytical photolysis
and stiff-exchange trajectories at the reference tolerances using both dense
and sparse Mozart solves. This does not yet cover Doolittle or constrained DAEs.

The case report now has 50 passing cases and 630 unported failures for both
compilers. Doolittle reconstruction checks cover dense, random sparse, diagonal,
and extreme initial scratch values, and both LU algorithms pass the analytical
photolysis trajectories. The 100-cell TS1 comparison remains passing.

The latest case report has 92 passing cases and 588 unported failures with both
compilers. The additional cases cover linear, equilibrium, and constraint-set
validation, residuals, Jacobians, sparsity, and multi-cell indexing. The DAE
integration cases are still unported.

The current case report has 107 passing cases and 573 unported failures on
both compilers. Additional checks cover consistent initialization, transactional
rollback, damping, row scaling, and single- and coupled-equilibrium DAE solves
with both linear algorithms. The 100-cell TS1 comparison remains passing.

The latest case report has 143 passing cases and 537 unported failures with
both compilers. New checks cover Mozart decomposition, dense and sparse linear
solves, and Markowitz fill reduction. The 100-cell TS1 comparison still passes.

The latest case report has 149 passing cases and 531 unported failures with
both compilers. New checks cover parameterized reactants and callback rates,
including runtime rebinding and generated mechanisms. The 100-cell TS1 check
continues to pass through the generalized multiplier path.

The latest case report has 158 passing cases and 522 unported failures on
both compilers. New checks cover process-set algebraic row masking and
finite-difference Jacobian comparison and sparsity verification. Both
compilers pass the 100-cell TS1 comparison.
