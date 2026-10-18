# Structure scans

Every native research scenario run is recorded in `scans/*.csv`, including failed runs.
Matching JSON files give exact command-line configurations. `tools/scan.py`
runs one native process at a time, obtains `elapsed_ms` from that process, and
invokes the **supplied** `../scenarios/score.py --json`. Its full JSON result,
operation counts, compiler label, source and binary hashes, and failure stderr
are retained. Accuracy is never replaced with a home-grown approximation.
Python prepares and scores experiments; it supplies no speed measurements.

The completed record contains **960 native research runs**. It retains 21 failed
legacy/default-baseline runs; no declared operating-point run failed. Final
10,000-cell measurements comprise 120 C runs (four points, six cases, five
repeats) and 30 supplied MICM baseline runs. Failed runs are never assigned a
speedup or included in timing medians.

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
| `07-caps-tolerance` | Intermediate step caps, then relative/absolute tolerance changes with cap held fixed |
| `08-stress-candidates` | Independent 48-hour scenarios with converged references; four candidate operating points |
| `09-groups` | 10,000 cells; only group size changes (16, 32, 64, 128), scoring the disclosed prefix |
| `10-copy-elision` | Repeat two leading groups after eliminating the unnecessary preserved-Jacobian copy |
| `11-sdirk-projection` | Retest SDIRK on the final engine after correcting its signed projected-error estimate |

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
timing is excluded from the work–precision figure and is not used as an
idle-machine performance claim. Early single-run times
guide candidate selection, while final operating-point comparisons use repeated
native measurements without concurrent build or test jobs.

MICM baseline runs retain the supplied runner's default one group per thread;
new C runs use the recorded group size. The unchanged MICM timer includes its
in-memory snapshot copies, whereas the native C timer excludes them. Baseline
comparisons consequently reflect implementation and grouping as well as the
numerical method. Formula comparisons within a structural scan keep these
execution choices fixed.

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

The step-cap scan gave a useful improvement in the global error/cost tradeoff.
For RODAS3 at rtol `0.003`, atol `1e-18`, a 300-second cap improved the minimum
overall/species digits from 2.484/1.546 to 3.086/2.060, with a six-scenario
total of 4,279 ms. A 200-second cap gave 3.433/2.509 digits at 5,174 ms, and
a 100-second cap gave 3.867/3.065 digits at 7,583 ms. RODAS3.1 remained weaker
on the 5-minute case unless capped more tightly. These results motivate a
step-limited RODAS3 recommendation; the cap limits internal integration steps
and never changes the prescribed photolysis update times.

The three independent 48-hour cases vary latitude, season, initial solar time,
NOx scaling, and random seed. Their two tight reference integrations agree to
10.193–10.412 overall digits and at least 9.494 digits per species. All four
candidate settings passed their intended accuracy levels. Reference generation
and bidirectional convergence scores are recorded in `stress-reference.csv`
and `validation/stress-references.json`.

At 10,000 cells, scan 09 initially favored groups of 16 over 32, 64 and 128.
The fresh-Jacobian path then removed an unnecessary unfactored matrix allocation
and copy. Scan 10 compared groups 16 and 32 again, with two repeats per scenario:
the sums of the six per-scenario median times were 50,199 and 48,628 ms, and
their geometric means were 7,458 and 7,334 ms. Group 32 was faster in four of
six scenarios, but slower on the short Tier 1 case. The final points use 32;
the advantage is modest (1.7% by geometric mean), not a universal group-size
rule. Every full output hash and every operation counter matched the respective
pre-copy-elision run across all 24 comparisons; see
`validation/copy-elision.json`. Preliminary GCC point checks with group 16 are
retained with `-group16` suffixes, separately from the final group-32 checks.

The final review found that SDIRK stored `embedded - candidate`, while the shared
projection correction assumed `candidate - embedded`. A scalar regression with
`y' = -100*y`, `y(0)=1`, `h=1`, and atol `2.25` reproduced a false acceptance:
the projected normalized error should be 1.022886, but the reversed sign gave
0.983723. C and the independent prototype now use the consistent signed error,
and both reject that trial. Both compiler suites were rerun. Scan 11 repeats
SDIRK's original six-case setting on the final engine: 48,798 ms total,
4.420 minimum overall digits and 3.408 minimum species digits. Other engine
optimizations also occurred since scan 01, so the time difference is not an
isolated performance effect of the sign correction. The superseded six SDIRK
observations stay in scan 01 but are excluded from the plot.

The RODAS3 operating-point formulas are unaffected. Earlier point measurements
remain under `pre-sdirk-fix-*`; final accuracy and five-run timing series were
refreshed against the corrected binary. The unchanged MICM baseline binary and
its complete five-run series remain the same reference comparison.
All 192 before/after operating-point output hashes and counter records match
exactly, covering both compilers, the full accuracy cases and the large timing
cases. Each of the 24 final large configurations also has identical output
hashes and counters across its five repeats. The comparison record is
[validation/sdirk-fix-operating-points.json](validation/sdirk-fix-operating-points.json).

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

## Reproduction and retained records

Use the pinned local conda environment (`make env`) and build the native runners
before beginning timing. Every Python command below uses `.conda/bin/python`.
The default scanner output root is `runs/scans`; `--run-root` separates series
that share configuration IDs. `--keep` retains the large native CSV files as
well as the stdout, stderr and scorer JSON. CSV logs append and must have a
fresh filename for a new experiment; they are not silently overwritten.

```sh
make bench bin/micm_scenario CC=clang
.conda/bin/python tools/prepare_experiments.py --large --stress

# Two independent tight reference integrations; keep both sets of outputs.
.conda/bin/python tools/scan.py scans/stress-reference.json --log runs/recheck-stress-reference.csv --scenario-root runs/stress --scenarios winter-night-48h summer-polluted-48h tropical-clean-48h --keep
.conda/bin/python tools/finalize_references.py --log runs/recheck-stress-reference.csv --manifest runs/recheck-stress-references.json
```

`finalize_references.py` reads the two named reference directories under
`runs/scans`, checks convergence in both directions with the supplied scorer,
and installs the tighter reference only inside `runs/stress`. It adds convergence
scores to the selected log and writes the requested manifest. With no arguments
its destinations are the canonical `scans/stress-reference.csv` and
`validation/stress-references.json`. Reference work is completed before any
candidate timings.

The final clang series is generated sequentially with:

```sh
.conda/bin/python tools/scan.py scans/operating-points.json --log scans/final-1000.csv --run-root runs/final-1000 --keep
.conda/bin/python tools/scan.py scans/operating-points.json --log scans/final-stress.csv --run-root runs/final-stress --scenario-root runs/stress --scenarios winter-night-48h summer-polluted-48h tropical-clean-48h --keep
.conda/bin/python tools/scan.py scans/operating-points.json --log scans/final-10000.csv --run-root runs/final-10000 --scenario-root runs/large --reference-root ../scenarios --prefix-cells 1000 --repeat 5
.conda/bin/python tools/scan.py scans/final-baseline.json --log scans/final-baseline-10000.csv --run-root runs/final-baseline --scenario-root runs/large --reference-root ../scenarios --prefix-cells 1000 --repeat 5
```

These are the original filenames; choose new log names when repeating completed
experiments. Full 10,000-cell output hashes are logged even though only the
disclosed 1,000-cell prefix has reference accuracy. Each run's finite-state and
minimum-concentration diagnostics cover all cells. The final plotting/report
scripts read the canonical final logs and fail if required measurements or
accuracy gates are missing. GCC checks use the same point JSON, with
`--compiler gcc-15`, after rebuilding; they are accuracy checks, not the primary
timing ranking. Compiler labels describe the built binary and do not invoke a
compiler themselves.
