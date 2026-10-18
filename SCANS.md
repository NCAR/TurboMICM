# Structure scans

Every native scenario run is recorded in `scans/*.csv`, including failed runs.
Matching JSON files give exact command-line configurations. `tools/scan.py`
runs one native process at a time, obtains `elapsed_ms` from that process, and
invokes the **supplied** `../scenarios/score.py --json`. Its full JSON result,
operation counts, compiler label, source and binary hashes, and failure stderr
are retained. Accuracy is never replaced with a home-grown approximation.
Python prepares and scores experiments; it supplies no speed measurements.

The initial scans use all six disclosed 1,000-cell scenarios and ten threads.
The standard adaptive controls are cell-wise maximum RMS, fresh Jacobians,
integral control, warm step prediction, and groups of 32. A structural row
changes one choice relative to its stated parent. Tolerance sweeps hold the
algorithm fixed and vary one accuracy control. They are not evidence that a
tolerance value itself guarantees a number of digits.

## Scan sequence and controls

| Files | Question / controlled changes |
| --- | --- |
| `01-families` | Same rtol `1e-4` and atol `1e-18`: ROS2, ROS3, ROS4, RODAS3, RODAS4, SDIRK2; separate unmodified legacy baseline |
| `02-tolerances` | Per-family rtol `1e-4`, `1e-3`, `1e-2`, fixed atol; adds published RODAS3.1 coefficients |
| `03-structure` | RODAS3, rtol `1e-3`, atol `1e-18`; controllers, restart policy, error norm, Jacobian age, LU holding, absolute floor, projection bound |
| `04-workprecision` | Intermediate relative tolerances for promising methods; atol changed within fixed-rtol series; repeat of the structural baseline after norm-loop reordering |
| `05-baselines` | Native MICM three-stage defaults and a tighter robust baseline, plus the same tight setup in the C port |
| `06-stepcaps` | Per-method/per-tolerance pairs, changing only maximum internal step to 300 or 100 seconds |

In scan 03, `structure-base` is the parent unless stated otherwise.
`restart-large` differs only in starting h from `restart`; `h211b-b1` differs
only in b from `h211b`; `h211b-b1-k2` differs only in k from `h211b-b1`;
`lu15` differs only in LU holding from `jac4`. Each absolute-tolerance row
changes only atol; no species-specific tuning is done.

Early pilot binaries evolved in explicitly mathematical-equivalent ways:
exact RHS reuse on rejected trials, followed by SIMD across cells in the error
norm. The projection correction was also charged to the embedded error before
acceptance. Source hashes distinguish these stages. The norm-loop change
produced byte-identical output CSV hashes on all six baseline scenarios.
Scan 04's first configuration overlapped reference compilation; its exploratory
timing is not used as an idle-machine performance claim. Early single-run times
guide candidate selection, while final operating-point comparisons use repeated
native measurements without concurrent build or test jobs.

## Findings from the controlled scans

At the common tight setting in scan 01:

| Formula | Sum of six elapsed_ms | Minimum overall digits | Minimum species digits |
| --- | ---: | ---: | ---: |
| ROS2 | 58,848 | 3.625 | 2.639 |
| ROS3 | 24,663 | 4.466 | 3.352 |
| ROS4 | 31,492 | 4.597 | 3.718 |
| RODAS3 | 14,738 | 4.150 | 3.008 |
| RODAS4 | 20,388 | 4.611 | 3.589 |
| SDIRK2 | 55,957 | 4.420 | 3.408 |

These are single observations at a common tolerance, **not** an equal-accuracy
ranking. RODAS3 is a strong starting point: fewer stages alone did not win.
The second-order methods' embedded estimates demanded many more accepted steps.
The retained default three-stage method failed three disclosed scenarios and
achieved only 1.215 overall digits on the 5-minute diurnal case.

At rtol `1e-3`, RODAS3's minimum overall accuracy fell to 2.837 digits and its
minimum species accuracy to 1.962. That tolerance cannot be called a universal
three-digit setting. RODAS3.1 reduced cost in several comparisons but likewise
needed its own work–precision curve; at rtol `5e-4`, one species fell below
two digits. Favorable stage coefficients do not remove the need for validation.

The main structural observations at the scan-03 setting were:

- Restarting at 0.01 s instead of carrying h increased total time from 5,523 to
  11,963 ms, while improving the minimum overall accuracy from 2.837 to 3.421.
  Starting restarts at 1 s reduced the total to 9,065 ms. Restart cost matters
  particularly for frequent chemistry calls.
- PI and tested H211b controllers did not improve the combined accuracy/cost
  tradeoff here. Their histories are discarded at discontinuous forcing
  boundaries. This measured result does not contradict improvements reported
  for other chemistry systems and controller parameters.
- RMS pooling across cells was somewhat cheaper but weakened the worst species
  from 1.962 to 1.679 digits. The infinity norm was more accurate and much more
  expensive. Cell-wise maximum RMS was retained for difficult-cell protection.
- Four-step Jacobian lagging reduced updates from 121,563 to 78,976, but raised
  factorizations to 123,912 and weakened the worst species to 1.695 digits.
  Its 5,844 ms total did not beat the 5,523 ms fresh-matrix baseline. Eight-step
  lagging reduced accuracy further. Age four plus a 1.5 LU holding band cut
  factorizations to 106,854, but its 5,819 ms total still did not establish a win.
  These are lagged Rosenbrock experiments, not methods with proved W-order.
- Raising atol from `1e-18` to `1e-16` reduced total time to 4,376 ms but weakened
  the worst species to 1.675 digits. The absolute floor is an accuracy choice,
  especially for radicals and low-abundance species, rather than a free speedup.
- Increasing the permitted tiny projection from 0.1 to 1 times atol did not
  improve reported accuracy and did not establish a time benefit. The more
  conservative default was retained.

## Scope of the conclusions

This is a systematic search over two families, seven formulas, three
controllers, three error norms, restart policies, Jacobian/factorization reuse,
tolerance floors, and internal step limits. It is not a proof of global
optimality. BDF/NDF, exponential, Patankar, partitioned/multirate, and mixed
precision implementations were not evaluated. Keeping the sparse full mechanism
and double precision reduced the number of new failure modes while testing
several materially different integration structures. There is no mechanism
reduction or learned surrogate.

The independent dense prototype agrees with every C method to double-precision
roundoff and matches all counters. Formal order tests, stiff analytical tests,
and independent 48-hour scenarios supplement the disclosed benchmark. Final
selected configurations, repeated measurements, and scope of generalization
are reported in `OPERATING_POINTS.md`, `ACCURACY.md`, and `RESULTS.md`.
