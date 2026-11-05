# The paper's configuration in the compiled library, 2026-11-04 and 05

The timing of the retaining solver of "A Fast ODE Solver for Atmospheric
Chemistry: ROS3 with H211b Control and Retained Steps" (the short JAMES
manuscript, `paper_efficient_ode/` in the ODE repository) in the reference
solution's library, beside the library's own restarting loop and MICM's
loop, on the 10,000-cell versions of all eleven Tier 3 scenarios. The
session ran from 23:34 on 4 November to 00:52 on 5 November 2026 on an Apple
M5 Pro laptop (see `environment.md`; not the machine of the ranking), at ten
threads, three rotating rounds per scenario, every run passing the gate.

## Programs (`programs.txt`)

- `micm/C2`, `micm/C3`: MICM's three-stage Rosenbrock in its own driver, at
  the ranking's baseline settings (relative 1e-2 and 3e-4, absolute 1e-15,
  32-cell groups).
- `loop/C2`, `loop/C3`: the library's ROS3 restarting loop (`warm_start=0`,
  `quantize=0`, `controller=micm`) at the paper's least-cost restarting
  points, relative 3e-2 and 3e-4 with absolute 1e-15.
- `retaining/C2`, `retaining/C3`: the paper's retaining solver, ROS3 under
  H211b (b = 4, safety 0.9) with the carried and quantized step and the three
  interval rules (`exact_end=1`, `reject_elementary=1`, `strict_min=1`), at
  relative 1e-2 and 5e-3 with absolute tolerance rtol x 1e-14.

Branch `tier3/fable/paper` (`296beb4` for the Rosenbrock loops; `7a23d26`
adds `--method be` and `substeps=N`, not timed). Both loops reproduce the
archived Python solver's operation counts exactly (`agreement.md`), and the
backward Euler support reproduces its digits (`agree_be.py`).

## Results (`session-summary.md`, `out/<set>/<scenario>/`)

Sums over the eleven scenarios of the median elapsed time (s):

| Program | All | Disclosed | Held out |
|---|---:|---:|---:|
| micm/C2 | 333.14 | 158.28 | 174.86 |
| micm/C3 | 402.24 | 187.28 | 214.96 |
| loop/C2 | 74.47 | 34.55 | 39.92 |
| loop/C3 | 90.85 | 42.17 | 48.67 |
| retaining/C2 | 21.12 | 9.64 | 11.48 |
| retaining/C3 | 30.60 | 13.70 | 16.90 |

| Ratio of summed times | All | Per-scenario range |
|---|---:|---:|
| restarting loop / retaining solver, criterion 2 | 3.53 | 1.49 to 6.74 |
| restarting loop / retaining solver, criterion 3 | 2.97 | 1.28 to 6.72 |
| MICM / retaining solver, criterion 2 | 15.78 | 6.72 to 30.24 |
| MICM / retaining solver, criterion 3 | 13.14 | 5.75 to 29.98 |
| MICM / restarting loop (the kernels), criterion 2 | 4.47 | 4.26 to 4.78 |
| MICM / restarting loop (the kernels), criterion 3 | 4.43 | 4.30 to 4.58 |

The time ratio between the two loops, which share the kernels, equals the
paper's arithmetic ratio (3.6 and 3.0 at the two criteria) and follows the
paper's per-scenario attempt ratios (1.3 on Diurnal 30 min, 6.6 on
Benchmark 30 s). The factor against MICM's loop includes the generated
kernels, worth 4.4 to 4.5 on this machine.

## Files

- `programs.txt`, `time_paper.sh`: the session; `time.log` its log;
  `summarize.py` the tables above from `out/`; `apply_timing.py` the
  script that wrote the numbers into the paper.
- `out/<set>/<scenario>/timings.csv`, `summary.md`, `environment.txt`: every
  run, as `admin/time_entries.py` records them (the `runs/` scratch
  directories are not kept).
- `agree.py`, `agreement.md`: the agreement of the C loops with the Python
  solver; `agree_be.py`, `time_be.sh`: the backward Euler agreement and a
  timing script that was not run (the paper keeps its arithmetic comparison).
- `environment.md`: machine, compiler, library, baseline, scenarios, and
  scoring references.
