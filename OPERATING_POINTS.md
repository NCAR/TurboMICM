# Operating points

The conservative configuration below has passed every disclosed scenario.
Additional throughput candidates are still being measured and are not declared.
All options are runtime parameters; no compiler flag changes are used.

## Conservative — recommended pending throughput study

```
bin/ts1_scenario --method ros3 --option rtol=1e-4 --option atol=1e-18 SCENARIO
```

Other experimental defaults apply: max-cell RMS norm, positivity rejection,
exact forcing reuse on retries, a current Jacobian for every accepted step,
I controller, warm step proposals, no step deadband, 32-cell groups (shrinking
for small inputs), and all available cores. No approximate Jacobian reuse.

Expected accuracy on the disclosed cases: at least **4.47 digits overall** and
**3.35 digits for every species**, hence a margin above both ranking levels.
The initial family scan used `reuse_rhs=0`; forcing reuse gives identical
trajectories and only removes redundant evaluations at unchanged y. The point
will be rerun in the final performance measurements before selection is final.
This conservative margin is useful for an atmospheric model because local
weights alone do not bound accumulated species errors over a full day.
