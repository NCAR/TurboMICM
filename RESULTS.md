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
