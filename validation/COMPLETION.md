# Final evidence audit

Completed on 2026-10-18, after the final numerical change in commit `a307b8a`.
The remaining changes contain reports, plots, input provenance and measurements;
they do not alter native numerical code.

| Requirement | Final evidence |
| --- | --- |
| Retain Tier 1 behavior | Both compilers pass all 680 exact supplied case IDs and both TS1 gates at relative tolerance `1e-10`; [gate record](tier3-gates.json) |
| Preserve the inherited chemistry and solvers | All 24 inherited files under `src/`, `configs/`, and `include/chem_array.h` match commit `0dc7f67` byte for byte; the public header adds the new API |
| Native runtime selection | Seven formulas in the new ODE driver, five retained Rosenbrock choices, documented controls and scenario CLI in [DESIGN.md](../DESIGN.md) |
| Independent implementation checks | 147 C/prototype comparisons per compiler with identical counters, seven stiff analytical cases, linear-chain order checks, failure tests, projected-error regression and scenario-contract checks; [index](README.md) |
| Systematic numerical research | 960 recorded native runs, including 21 failed legacy/default-baseline runs, in [scans/](../scans/); [controlled comparisons](../SCANS.md) |
| Declared configurations | Four points, all passing their target level on six disclosed and three independent full cases with both compilers; [commands](../OPERATING_POINTS.md) |
| Primary-size timing | 120 successful C runs: four points, six 10,000-cell cases, five repeats each; 30 separate MICM baseline attempts; [results and timing scope](../RESULTS.md) |
| Repeatability | Every final large configuration has identical output hashes and counters across its five repeats; all 192 comparisons with measurements before the SDIRK correction also match; [record](sdirk-fix-operating-points.json) |
| Input provenance | SHA-256 hashes for 80 files and exact disclosed-prefix/forcing/settings checks; [manifest](input-manifest.json) |
| Accuracy justification | Literature, per-species constraints, projection measurements, independent reference convergence and generalization limits in [ACCURACY.md](../ACCURACY.md) and [REFERENCES.md](../REFERENCES.md) |
| Portable execution | ISO C11, double precision, standard C/libm and CPU OpenMP; inspected native includes, pragmas, source constructs and linked libraries; no target, floating-point, LTO or profile optimization flags |
| Reproducible reports | [Report generator](../tools/report_results.py), [plot generator](../tools/plot_scans.py), pinned [environment](../environment.yml), and PNG/SVG work–precision figures with the PNGs visually inspected |

Final native-source digest, using `tools/scan.py`'s source-file selection and
hashing procedure:

```text
29d95d1c369cfad7112351a42c8a73bbc2eadf583dd3b1061c1fe32edb9387bd
```

Final clang scenario executable SHA-256:

```text
3fcb362a893621e335749bd571f730b6e5e28d60c1f09b7c00eb02e0e929bc46
```

The six final C series (disclosed, stress, large, both GCC accuracy series and
corrected SDIRK) match the current native-source digest. The clang rows match
the current scenario executable. The current Tier 1 executable matches its
last clang gate record, and the current MICM scenario executable matches all
30 baseline measurements. Executable hashes identify this build; rebuilds need
not produce byte-identical binaries.

The scenario C executable and library link only the system C runtime and
OpenMP runtime on this host. All listed native-source checks passed on the
supplied macOS machine; Linux execution was unavailable. Accuracy for the
10,000-cell timing inputs is checked only on the disclosed 1,000-cell prefix.
All cells still undergo native finite/nonnegative-state checks. No claim is
made that unseen scenarios will necessarily satisfy the measured accuracy
levels, or that machine-wide idleness was independently certified.
