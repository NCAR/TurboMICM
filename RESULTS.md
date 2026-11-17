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

The latest case report has 180 passing cases and 500 unported failures with
both compilers. The additional cases cover multi-cell state storage, ownership,
copy/move behavior, named concentrations and custom rates, and invalid input.
Both compilers continue to pass the 100-cell TS1 comparison.

The latest case report has 202 passing cases and 478 unported failures on
both compilers. Additional analytical trajectories exercise all five
Rosenbrock methods and backward Euler with six rate families, including the
fast equilibrium variants, for three independently varying cells. Rate values
and every time step are checked against the reference formulas and tolerances.
Both compilers still pass the 100-cell TS1 comparison.

The latest case report has 203 passing cases and 477 unported failures on
both compilers. The random process-set case now covers the reference system
sizes, with additional forcing, Jacobian, and structural checks. A fixed
process fixture also checks the reference forcing and Jacobian, but its full
case remains unported because parameterized-species name enumeration is not
yet represented in the C mechanism. TS1 still passes with both compilers.

The latest case report has 205 passing cases and 475 unported failures on
both compilers. The fixed process-set fixture now exercises generated
parameterized-species metadata, forcing, and Jacobian values. The 100-cell TS1
comparison remains passing on both compilers.

The latest case report has 208 passing cases and 472 unported failures on
both compilers. The additional checks cover every published Robertson and
HIRES output time in the reference, using each available Rosenbrock method
and both LU algorithms; HIRES also checks backward Euler. The 100-cell TS1
comparison remains passing with both compilers.

The latest case report has 211 passing cases and 469 unported failures on
both compilers. New checks cover all E5 and Oregonator reference output times
for the Rosenbrock methods and both LU algorithms, plus the 18,000-substep
backward Euler Oregonator checks. Both compilers pass the 100-cell TS1 check.

The latest case report has 216 passing cases and 464 unported failures on
both compilers. Added checks cover the surface reaction trajectory, the
terminator mechanism across changing illumination and cells, and monotonic
integrated reaction rates. Both compilers pass the 100-cell TS1 check.

The latest case report has 223 passing cases and 457 unported failures on
both compilers. New Chapman checks cover 100-step oxygen and spectator
conservation, rate constants, forcing, and analytical Jacobians. Both
compilers continue to pass the 100-cell TS1 check.

The latest case report has 230 passing cases and 450 unported failures on
both compilers. New external constraint checks cover residual/Jacobian finite
differences, conservation, row replacement, and multiple cells. Custom
nonlinear projection and composition with built-in constraints are also
exercised. Both compilers pass the 100-cell TS1 comparison.

The latest case report has 235 passing cases and 445 unported failures on
both compilers. New checks compare kinetic and composed DAE systems, built-in
and callback constraints at each step, multiple equilibria, and temperature-
dependent parameter updates. Both compilers pass the 100-cell TS1 comparison.

The latest case report has 273 passing cases and 407 unported failures on
both compilers. Additional matrix cases cover row-major and interleaved group
layouts, initialization, indexing, conversion, arithmetic, column operations,
and printing. Both compilers pass the 100-cell TS1 comparison.

The latest case report has 305 passing cases and 375 unported failures on
both compilers. New checks cover reusable row callbacks, multiple matrices,
temporary vectors, size errors, and strict and padded reductions for matrix
layouts. Both compilers pass the 100-cell TS1 comparison.

The latest case report has 339 passing cases and 341 unported failures on
both compilers. Additional matrix tests cover reusable signatures, matrix
column checks, vector and row-count mismatches, empty arguments, and mutable
vectors. Both compilers pass the 100-cell TS1 comparison.

The latest case report has 377 passing cases and 303 unported failures on
both compilers. New sparse checks cover CSR and CSC structures, grouped
blocks, duplicate and invalid coordinates, diagonals, copying, and printing.
Both compilers pass the 100-cell TS1 comparison.

The latest case report has 419 passing cases and 261 unported failures on
both compilers. New sparse tests cover reusable block views, row callbacks,
multiple structures, dense and vector interoperability, index lookup, fill,
and copy. Both compilers pass the 100-cell TS1 comparison.

The latest case report has 505 passing cases and 175 unported failures on
both compilers. New checks cover sparse function block/vector dimensions,
reusable functions, empty blocks, and incompatible group widths. Both
compilers pass the 100-cell TS1 comparison.

The latest case report has 526 passing cases and 154 unported failures on
both compilers. New grouped-view checks cover dense and sparse lanes, const
access, mixed views, partial groups, and fill/copy across CSR and CSC layouts.
Both compilers pass the 100-cell TS1 comparison.

The latest case report has 533 passing cases and 147 unported failures on
both compilers. The additional DAE checks cover constrained overshoot across
four Rosenbrock methods, algebraic tolerance sensitivity, and normalized
error with constraints. Both compilers pass the 100-cell TS1 comparison.
