# Numerical structure study

The selected formula is the stiffly accurate four-stage, third-order Rosenbrock method (`rodas3`), with revised startup, retries, positive endpoint checks, and error control. Five tolerance settings cover the useful work–precision range. `scans/` contains the attempted configurations and individual native runs, including configurations that failed or produced inaccurate finite answers. `REFERENCES.md` gives the mathematical sources and their use.

The study deliberately separates structural comparisons from tolerance adjustment. A common tolerance is not a common achieved accuracy: both are reported. Python prototypes establish arithmetic and operation-count agreement on small systems, while only the native executable supplies times. The final figures and reported 1,000-cell times use five serial repeats, not timings from the small screening inputs.

## Family and order

The first scan uses every disclosed scenario, the first 128 cells, `rtol=1e-4`, `atol=1e-18`, maximum cell RMS control, a 1-second initial proposal, and the same growth/safety bounds. Only the formula or BDF maximum order changes. No species is removed, lumped, or treated by hand. The legacy rows use their original controls and are explicitly separate anchors.

| Formula | Minimum overall digits | Minimum species digits | Accepted group steps, all six scenarios |
| --- | ---: | ---: | ---: |
| ROS2 | 3.684 | 2.676 | 210046 |
| ROS3 | 4.593 | 3.463 | 53737 |
| ROS4 | 4.893 | 3.859 | 53943 |
| RODAS3 | 4.331 | 3.142 | 27226 |
| RODAS4 | 4.704 | 3.646 | 21991 |
| BDF, order ≤1 | 1.977 | 1.117 | 184821 |
| BDF, order ≤2 | 2.996 | 1.894 | 106066 |
| BDF, order ≤3 | 3.288 | 2.303 | 99068 |
| BDF, order ≤5 | 3.293 | 2.304 | 96570 |

RODAS4 takes fewer steps than RODAS3, but six stages and six forcing evaluations per step erase that benefit at the required accuracy. ROS3 is accurate but takes more steps. ROS2 is robust, as anticipated by Verwer et al., but its lower-order embedded estimator is costly here. BDF's history and nonlinear solves do not pay off when photolysis repeatedly restarts the integration; reducing Jacobian refresh frequency helps only part of that cost. This is evidence about this independent BDF implementation and these forcing intervals, not a claim that Rosenbrock always beats production BDF solvers.

The next tolerance scan checks competing families at different achieved accuracies. RODAS3 at `rtol=1e-3` is near the level-three boundary, while `3e-4` gives margin. RODAS4 at `1e-3`, and ROS3 at `3e-3`, also reach roughly that accuracy range, with higher work. The full-cell sweep later confirms the choice. Tightening relative tolerance is generally effective but the error is not strictly monotone at every scenario: adaptive step histories and positivity projection matter. The full sweep is therefore used instead of an interpolation between two tolerance settings.

## Norm and absolute floor

Holding RODAS3 and `rtol=3e-4` fixed, group RMS reduces work but lowers the minimum species accuracy from 2.632 to 2.379 digits in the prefix scan. Maximum component control raises it to 3.716 digits at appreciably greater work. Maximum cell RMS is the default compromise: adding easy cells cannot hide a difficult cell's error. The additional maximum-norm tolerance sweep checks that this conclusion is not merely a difference in nominal tolerance.

At `rtol=1e-3`, raising `atol` from `1e-18` to `1e-16` loses accuracy in low-concentration species. The absolute floor remains `1e-18`: it is four orders below the scoring threshold, rather than the inherited `1e-3` floor that effectively disables chemical error control in these units. Accuracy is still measured against references; that separation of scales alone is not a global-error bound.

## Startup, controller, and retries

Initial proposals from 0.001 to 100 seconds produce similar converged accuracy after retries. Only the first call uses this proposal; subsequent calls retain their adaptive next step. This avoids restarting a nearly equilibrated box at a fixed millionth of every outer interval. Step proposals never cross an outer forcing boundary, and a changed light factor clears controller/multistep history.

The PI policy increases work at common tolerance. A separate H211b scan varies its smoothing and exponent parameters, motivated by Dreger et al. (2025); it supplies no consistent improvement over proportional control here. The literature result therefore does not transfer automatically to this controller, mechanism, and restart policy. Raising the proportional safety factor to 1.3 increases rejections rather than improving efficiency. Growth bounds 2, 6, 10, and 100 do not improve on the default 3 sufficiently to justify changing it. Every controller retains the same error and positivity acceptance tests.

Nonfinite stage results reject the trial and reduce the step instead of immediately aborting the interval. Large negative endpoints also reject the trial. Small accepted negatives are projected with explicit instrumentation; the largest measured corrections for the declared points are reported in `ACCURACY.md`. This policy is tested on a stiff conservative transfer and by refinement against the scenario references. It does not make the formula intrinsically positivity-preserving or exactly conservative after projection.

## Jacobians and factorizations

Increasing RODAS3 Jacobian age from 1 to 2, 4, and 8 lowers the prefix minimum overall digits from 3.611 to 3.607, 3.369, and 3.218, respectively; the corresponding species minima are 2.632, 2.493, 2.244, and 2.090. Arbitrary approximate Jacobians need not preserve its full formal order. Holding nearly unchanged step sizes enables factor reuse but adds enough steps that it does not improve the measured work–precision result. No declared point lags a Jacobian.

ROS2 remains a useful controlled W-method experiment. With a 1.2 step deadband, matrix ages 1, 4, and 16 maintain similar accuracy while reducing linear-algebra work. They still require substantially more steps than RODAS3 at the target accuracies. BDF's approximate matrix reuse is guarded by the nonlinear convergence test and shift-change threshold; increasing the age to 20 helps factorization cost, but its startup/iteration cost remains high. Failed and compiler-overlapped screening measurements remain identifiable in the logs.

## Grouping and validation scope

The prefix grouping scan changes 32 cells to 8, 16, or 64, including partial groups where applicable. The error norm remains per-cell and no group-specific species selection exists. These short inputs expose controller overhead and uneven worker occupancy; they are insufficient to decide throughput on 10,000 cells. A separate larger-input scheduling/grouping scan tests that issue before final throughput measurements.

The final five-point study runs all 1,000 disclosed cells and every supplied snapshot, with five native repetitions. Default MICM3 is also run independently, and so is an accuracy-controlled MICM3 baseline that succeeds on all six scenarios. Default finite output is not treated as accurate: several defaults miss both targets or fail outright. Additional locally generated 48-hour inputs vary seed, season, latitude, and NOx; their reference convergence is checked independently. Neither these checks nor the disclosed cases prove performance on unknown future inputs.

The study keeps double precision and solves the coupled full mechanism. Splitting and mixed precision are not claimed to have been scanned. They would introduce separate coupling and roundoff questions; the family/control/reuse comparisons already establish a useful improvement without those additional assumptions. There is no claim of a global optimum over all possible integrators.
