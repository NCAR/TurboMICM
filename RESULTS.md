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
