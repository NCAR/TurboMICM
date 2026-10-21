# TurboMICM delivers a 5.11x MICM speedup in CheMPAS-A TS1

## Latest result: the Tier 1 gains reach the coupled model

**Good news for the optimizer team: optimized Fable Tier 1 delivers 5.11x
faster MICM calls and 3.43x faster overall execution in CheMPAS-A TS1.**
That includes the TMT interface's validation, packing and unpacking overhead.
The optimization work is paying off in the coupled application.

The earlier TS1 slowdown measured the large-state compatibility configuration.
Selecting 32-cell groups, width 2, fast kernels and minimum-fill ordering reduced
MICM time from 115.624 s in that configuration to **16.893 s**. The saved MUSICA
baseline took **86.362 s**. Three-stage Rosenbrock, the MICM step controller,
rtol 1e-8 and all absolute tolerances of 1e-14 mol/m3 were retained. Tier 3
profiles and relaxed tolerances were not used.

| Scope | Saved MUSICA seconds | Optimized TMT seconds | MUSICA / TMT |
|---|---:|---:|---:|
| MICM library calls | 86.362159 | 16.893161 | **5.1123** |
| TUV library calls | 40.956201 | 6.224055 | **6.5803** |
| Entire process | 151.582252 | 44.182417 | **3.4308** |

This is the 2026-10-21 run: 209 species, 547 reactions, 112,000 chemistry cells,
50 six-second steps, 18 MPI ranks and one OpenMP thread per rank on an Apple M5
Pro. One model job ran at a time. MUSICA was not rerun; the exact saved inputs
and outputs were reused. These are single observations against that saved
baseline, with generation/compilation excluded from the timed run. Component
times are maximum accumulated rank times and must not be added together.

Per-call MICM instrumentation is now available. The first call's maximum rank
time was **0.506 s**, followed by means of **0.319, 0.335, 0.349 and 0.359 s**
over calls 11–20, 21–30, 31–40 and 41–50. The benefit is present throughout the
run; there was no need to wait for a late performance transition.

Independent group error norms change the integration trajectory. The earlier
strict field-scaled 1e-7 agreement check is not met: maxima are 5.41e-4 after
two steps and 1.78e-5 after 50 steps. We retained that failed check and added a
tolerance-aware comparison after the smoke run, without changing solver
settings. Differences normalized by
`atol + rtol*max(abs(reference),abs(candidate))` in mol/m3 have RMS/max of
0.02898/1.665 after two steps and **0.003466/0.6958 after 50 steps**.
All final timed concentration differences are below one tolerance unit.
Nonchemistry fields are identical, photolysis agrees to 1.15e-10 field-scaled,
and there are no solver failures, negative outputs or clipping. This measures
cross-backend agreement, not independent global integration error.

The updated adapter is now installed for normal CheMPAS-A runs, with
`TMT_MICM_TIER1=1` enabled by `/Users/dfillmor/Data/CheMPAS/env.sh`. An additional
18-rank, two-step installation smoke test confirmed Tier 1 on every rank, and
the TMT chemistry/photolysis plus MUSICA MIEM coexistence check passed.
CheMPAS model Fortran and the executable are unchanged. The only subsequent
CheMPAS build change fixes a macOS library search path.

Code and evidence:

- TMT adapter: `chempas-a-tmt`, commit `69272ea`.
- Reusable runner and full result record: CheMPAS-A-Dev `develop`, `bce8390`;
  [report](../CheMPAS-A-Dev/test_cases/ts1_benchmark/TIER1.md) and
  [runner](../CheMPAS-A-Dev/scripts/benchmark_ts1_tier1.py).
- CheMPAS macOS installation fix: `chempas-a-tmt`, `800c7fcc`.
- Pinned Fable source: `cb0f8a3020f1c0b9ec7e0ad079e79b4f4a8f9bdd`.
- Raw TS1 outputs, private libraries and per-rank CSV traces:
  `/Users/dfillmor/Data/CheMPAS/ts1-tier1-20261021/`.
- Generated cache: `1238ed4958cde14b6e217878d91599f3ea7214b72e15c2c317a471aa7168e833`;
  2,293 LU nonzeros and 938 Jacobian jobs.
- Installed-runtime verification and previous installation backup:
  `/Users/dfillmor/Data/CheMPAS/tmt-runtime-update-20261021/`.

The optimized toy-ring case has not yet been rerun. Its older result remains
below as an open measurement, separate from this successful TS1 result.

## Historical compatibility-mode investigation

Recorded 2026-10-20, America/New_York (run timestamps are UTC, 2026-10-21).
Requested for the optimizer agents after a matched MUSICA/TMT comparison.

**Scope correction, 2026-10-21:** these runs use TMT's default compatibility
path, not the best-performing Fable configuration. Although the pinned source
is Fable Tier 3 `cb0f8a3`, CheMPAS requests `RosenbrockStandardOrder` and never
calls `configure_fable`. TMT selects three-stage Rosenbrock, width 1, exact
kernels and Markowitz ordering; each rank owns a roughly 6,000-cell state.
Published Tier 1 gains used width 2, fast kernels, minimum-fill ordering and
32-cell groups. Tier 3 gains also use selected solver/controller/tolerance
profiles. TMT exposes L2, L3r and fourth-order extrapolation through
`TmtConfigureFable`, but none is active here. Reproduce those settings with
an explicit accuracy comparison before treating these results as a regression
of the optimized solution. Longer runs alone do not activate those settings.

## Original toy-ring result

A connected 200-species synthetic mechanism ran correctly through both backends,
but the TMT MICM library call was slower. The longer comparison used 50 outer
steps of 6 seconds, 18 MPI ranks, one thread per rank, and one model job at a
time on an Apple M5 Pro (18 cores, 24 GiB). The TMT generated-code cache was warm.

| Scope | MUSICA seconds | TMT seconds | MUSICA / TMT |
|---|---:|---:|---:|
| Entire process run | 278.768587 | 303.528700 | 0.9184 |
| MICM library calls | 214.497469 | 275.433860 | 0.7788 |
| TUV library calls | 44.320496 | 6.284185 | 7.0527 |

Whole-process values are wall-clock seconds. Library values are the maximum
accumulated wall time over MPI ranks, measured around the public Fortran
`micm%solve` / `tuvx_solver%run` calls. They include the TMT interface's own
validation/copies but exclude CheMPAS state conversion and model diagnostics.
Different timer maxima can occur on different ranks; they are not additive.
There were 900 global MICM calls and 80,000 global TUV calls in each timed run.
This is one observation per backend, not a repeated statistical benchmark.

TurboMICM took 28.4% longer inside the MICM calls. TurboTUV saved 38.0 seconds;
the overall TMT run still took 24.8 seconds (8.9%) longer.

## The first step matters, but do not overinterpret its attribution

| MPAS integration steps (whole-step timers) | MUSICA seconds | TMT seconds |
|---|---:|---:|
| First step | 62.050984 | 84.878407 |
| Sum of steps 2–50 | 210.176620 | 211.409889 |
| Mean of final ten steps | 3.519248 | 3.262399 |

The first step accounts for 22.83 seconds, about 92% of the overall wall-time
penalty. However, these are whole integration-step measurements, including
photolysis. The MICM timers are cumulative, so the 60.94-second MICM-specific
penalty has **not** been assigned to individual steps. Do not infer from the
nearly equal later whole-step times that later MICM calls are equally fast:
TurboTUV offsets part of the chemistry cost.

A preceding two-step validation pair measured MICM at 63.882947 seconds
(MUSICA) versus 95.904983 seconds (TMT), and TUV at 1.499390 versus 0.258666
seconds. TMT's cold code-generation cost was outside the MICM timer. Those
short-run total wall times are not the warm-cache benchmark.

## Workload and numerical controls

- 100 monomers A_001–A_100 and 100 dimers A_101–A_200, plus inert background O3.
  Dimer i contains monomers i and i+1, with cyclic indexing.
- 100 association, 100 thermal dissociation and 100 photolysis reactions.
- Thermal reference rates span 1e-6 to 1e3 s^-1, uniformly in log space and
  independently shuffled (seed 1729). Association coefficients span 1 to
  1e9 m3 mol^-1 s^-1 at reference concentration 1e-6 mol/m3. B, C and E are zero.
- Relative solver tolerance 1e-8 and per-species absolute tolerance 1e-14
  mol/m3 in both backends. The generated module's manifest records the default
  construction tolerance 1e-6; CheMPAS overrides that to 1e-8 through
  `SetRosenbrockSolverParameters`. Both runtime logs confirm the override.
- CheMPAS selects `RosenbrockStandardOrder`. TMT uses generated width-1, exact
  mode, Markowitz ordering, three-stage Rosenbrock kernels. The generated
  mechanism reports 201 variables, 300 reactions, 995 LU nonzeros and 400
  Jacobian jobs. TMT's constructor checks that generated kernels were accepted.
- 1,600 columns x 70 atmospheric levels = 112,000 chemistry grid cells, batched
  over 18 ranks. Rank zero has 6,160 chemistry grid cells. This is a large
  standard-order batch, not the vector-ordered/small-batch benchmark path.
- 16 upper-atmosphere extension layers, 102 wavelength intervals, 100 photolysis
  outputs from 73 TS1 templates (the first 27 repeat). Host ozone is inert.
  92 channels become positive in this column; eight remain zero.
- All toy species start at q=2.4166666666666668e-8 kg/kg. At reference dry density
  1.2 kg/m3, monomers start at 1e-6 and dimers at 0.5e-6 mol/m3. The first step
  therefore includes a substantial fast chemical adjustment.

The largest timed chemistry difference was 2.8796153766e-13 at field scale;
photolysis was 5.8854233605e-11. This metric is maximum absolute difference
divided by the maximum absolute value of either field over all saved records,
not a pointwise relative error near zero. All nonchemistry output fields were
bitwise identical. Maximum relative drift of any mass-weighted monomer-unit
domain integral was 8.88e-16 for MUSICA and 1.11e-15 for TMT. There were no
negative concentration clipping messages or solver errors.

An initial checker incorrectly required every conserved sum to remain locally
uniform. Nonlinear tracer transport limiters can redistribute such sums while
conserving every tracer globally. The corrected check integrates each monomer
unit using dry-air mass; local departures are retained as diagnostics.

## Reproduce and inspect

Scripts and generated mechanism are tracked in **CheMPAS-A-Dev**, not CheMPAS-A:

- `../CheMPAS-A-Dev/scripts/generate_toy_mechanism.py`
- `../CheMPAS-A-Dev/scripts/init_toy_mechanism.py`
- `../CheMPAS-A-Dev/scripts/benchmark_toy_backends.py`
- `../CheMPAS-A-Dev/micm_configs/toy_ring200/`
- `../CheMPAS-A-Dev/micm_configs/TOY_MECHANISM.md`
- `../CheMPAS-A-Dev/test_cases/toy_ring200/` contains the compact results.

Complete local run artifacts:
`/Users/dfillmor/Data/CheMPAS/toy-ring200-20261020/`.
Start with `summary.json`, `timed-comparison.json`, and each backend's
`timed/{run_status.json,verification.json,log.atmosphere.0000.out,output.nc}`.
The cache under `mechanism-cache/` contains the generated C, normalized input,
compile flags and manifest. Cache key:
`e69c3c0bd0f5025b31e0d140f446d53b8db1e525046f13e77d2f9a2803d9df7f`.

The executables are
`/Users/dfillmor/Data/CheMPAS/tmt-validation-20261020/bin/atmosphere_{musica,tmt}_timed`.
Model code is on CheMPAS-A `chempas-a-tmt`, commit `4a8ad9fb` (instrumentation
was committed in `2c26c5ae`). The TMT interface is commit
`0f9ce1655c6dd24bdb45aaf3266d3860638bc3b9`, linked against pinned Fable
`cb0f8a3020f1c0b9ec7e0ad079e79b4f4a8f9bdd`. MUSICA is pinned to
`1403e3d22717bc87f3bf9d0aa591caf039c92bbc`. Compilers are GNU 16.2, optimized;
generated kernels use `-O3 -ffp-contract=off`, binary64, width 1, exact mode.

Source the local `/Users/dfillmor/Data/CheMPAS/env.sh`, then use the documented
benchmark command with the existing 18-rank mountain-wave fixture:
`/Users/dfillmor/Data/CheMPAS/tmt-validation-20261020/timing18/inputs/full_mountain_wave_abba`.
The environment now enables Tier 1; explicitly set `TMT_MICM_TIER1=0` after
sourcing it when requesting the compatibility configuration. The executables
above resolve the current installed TMT library. For exact historical runtime
reproduction, use the recorded commits or relink a copy against the archived
installation at
`/Users/dfillmor/Data/CheMPAS/tmt-runtime-update-20261021/previous-tmt-install/`.
The exact command/options are preserved in the run's `summary.json`. Use a fresh
output directory for a new measurement. Do not benchmark competing model jobs
concurrently, and exclude cold code generation from the timed comparison.

## Useful next measurements

1. Rerun the toy ring with the same optimized Tier 1 settings, retaining the
   positivity and 100 independent monomer-unit conservation checks. Its earlier
   timing should not be generalized to the optimized configuration.
2. Repeat the TS1 timings to measure variability. The saved MUSICA baseline uses
   one large standard-order state per rank; record that configuration explicitly.
3. Use the new per-call traces to separate startup from later work, then profile
   packing, parameter refresh and generated kernels individually if worthwhile.
   Per-call totals and solver counters are already recorded for Tier 1 TS1.
4. Evaluate Tier 3 as a separate experiment with an explicit accuracy target.
   Its method/controller/tolerance profiles require their own comparison.
5. Keep the width-1 large-state path as a compatibility optimization target.
   The TS1 result now provides a working optimized configuration to compare with.

The toy ring is much sparser than TS1 and has 100 independent conserved monomer
units; the two mechanisms remain distinct profiling cases.

## Original TS1 compatibility result: 2026-10-21

The original coupled TS1 benchmark also showed slower MICM calls through TMT. It used the
same executables, 18-rank grid, 50 six-second steps and numerical tolerances as
above, with 209 integrated species and 547 reactions. The generated cache is warm.

| Scope | MUSICA seconds | TMT seconds | MUSICA / TMT |
|---|---:|---:|---:|
| Entire process | 151.582252 | 144.105712 | 1.0519 |
| MICM library calls | 86.362159 | 115.624007 | 0.7469 |
| TUV library calls | 40.956201 | 6.813316 | 6.0112 |

The combined TMT run is 4.9% shorter because of the photolysis savings, despite
MICM taking 33.9% longer. Whole-step TMT means settle at 2.668, 2.653, 2.651 and
2.662 seconds over windows 11–20, 21–30, 31–40 and 41–50; the corresponding
MUSICA means are 2.837, 2.860, 2.865 and 2.868 seconds. This is not evidence of
a later MICM-only speed advantage: the per-step measurements include TUV and
the MICM timer remains cumulative.

Configuration and reproducible scripts are in
`../CheMPAS-A-Dev/micm_configs/TS1_BENCHMARK.md` and
`../CheMPAS-A-Dev/test_cases/ts1_benchmark/`, committed and pushed on `develop`
as `c79f9afbbc4f32fb00443518313b6c3671a02cc1`. Full local artifacts are under
`/Users/dfillmor/Data/CheMPAS/ts1-20261021-closed/`. Cache key:
`c6620f3d27d6b1e7d32d48103330dc92536d01479e6482317fdceb4fc87bf1de`.
The generated kernel reports 209 variables, 547 reactions, 2,381 LU nonzeros
and 938 Jacobian jobs; width 1, exact mode, Markowitz ordering.

This is a prescribed-rate, no-emissions TS1 case. Of 123 photolysis reactions,
122 use TUV (72 direct templates plus 50 CAM aliases), and NO dissociation
`jno` uses the reference fixed value; `jno_i` ionization is not substituted.
The 13 SURFACE reactions retain zero host aerosol parameters. Seventeen `het*`
rates use reference constants, and `usr_CO_OH`/`usr_DMS_OH` remain zero as in
the standalone reference InitState. H2O is renamed `TS1_H2O` to retain a closed
water tracer instead of binding to dry meteorological qv. Ozone uses a USSA
vertical profile; other species use reference mole fractions. These differences
from a standalone TS1 workload matter when comparing speedups.

All four runs pass validation. The timed maximum field-scaled differences are
5.58e-12 for chemistry and 6.01e-11 for photolysis; other fields are bitwise equal.
There are no solver errors or negative clipping. These historical measurements
motivated the adapter work; the optimized TS1 outcome is now at the top of this
note. The original toy-ring measurement is retained until its Tier 1 rerun.
