# Tier 2's undisclosed inputs

Five held-out day-night scenarios and the TS1 ranking input, which the agents
never saw, used with the five disclosed scenarios to check every entry and to
rank them. They were fixed on 2026-10-22, before any Tier 2 run, in a git
repository with no remote, outside the run root (`~/NSF-NCAR/.tier2-heldout`),
and their hashes were recorded in [`admin/runs.md`](../../admin/runs.md) before
the first launch. They are published now that the last Tier 2 run, luna-high,
has ended.

- **`make_heldout.sh`**, committed there at 09:30:00 (`eb4ae5f`); its SHA-256,
  `c8dfef98d54b09c52910a764cffd7791836278323f73c7b9f6e325543365750b`, is that
  of the copy here. Claude Code (Opus 5.5) chose the parameters, orchestrating
  Tier 2 at the SSE's instruction. Tier 3's held-out set is public, so Tier 2
  has its own.
- **The ranking input:** `reference/make_cells_csv.py 10000 SEED` with the
  63-bit seed in `SEED`, drawn with Python's `secrets` and committed there at
  09:30:31 (`e10fd84`). The cells file's SHA-256,
  `be7ccfdbc8c74a68c1112122e14bb2407009eba50d84485a784af8597e263772`, and that
  of MICM's final state for it, `cb7d0cc7014b598515a70805364e3854d5f318ccc1466ebd90e01b1b5bac62ce`,
  were recorded before the first launch.

The scenarios vary what the agents' scenarios README announced: seeds,
latitudes, seasons, start times, and initial NOx, with steps of 5 to 30
minutes and runs of up to 48 hours, recording a snapshot every 6 hours.

| Scenario | Seed | Steps | Latitude, day of year | Start (local solar time) | Initial NOx | MICM's internal steps at the gate's settings | Reference convergence: overall (worst species) | MICM at the gate's settings, against the reference | MICM's default |
|---|---:|---|---|---|---|---:|---|---:|---|
| spring45-20min | 301 | 72 of 20 min (24 h) | 45° N, 110 | 06:00 | x1 | 60,435 | 10.20 (OClO 9.28) | 4.84 | 2.67 digits |
| tropics-polluted-15min | 302 | 96 of 15 min (24 h) | 5° N, 200 | 06:00 | x3 | 75,293 | 10.29 (HO2NO2 9.46) | 5.10 | fails |
| arctic-summer-10min | 303 | 144 of 10 min (24 h) | 70° N, 172 (midnight sun) | 06:00 | x0.5 | 103,211 | 10.10 (OClO 9.10) | 5.04 | 3.64 digits |
| south-48h-30min | 304 | 96 of 30 min (48 h) | 35° S, 15 | 06:00 | x2 | 95,163 | 10.21 (BrCl 9.18) | 4.97 | 2.36 digits |
| night-start-5min | 305 | 288 of 5 min (24 h) | 40° N, 300 | 18:00 | x0.1 | 132,855 | 10.44 (BrO 9.52) | 4.73 | 3.61 digits |

MICM's results at the gate's settings (relative 1e-6, absolute 1e-15 mol m⁻³,
32-cell groups) and the reference solutions (the six-stage Rosenbrock at 1e-10
relative and 1e-23 absolute, checked against 1e-9 and 1e-22) were made with
`tier2/admin/make_scenario_data.sh` from 09:40 to 09:57 on 2026-10-22; MICM's
default is its three-stage Rosenbrock with absolute tolerance 1e-3, also with
32-cell groups.

## Files

| File | Contents |
|---|---|
| `make_heldout.sh` | the five scenarios: `make_heldout.sh MAKE_SCENARIO DIR [CELLS]`, with `tier2/scenarios/make_scenario.py` |
| `SEED` | the ranking input's seed |
| `SHA256SUMS` | the checksums of the scenario files, MICM's results, the references, their logs, the ranking input, and MICM's final state for it, relative to `~/NSF-NCAR/.tier2-heldout` |
| `micm-data-2026-10-22.log` | the log of MICM's runs and the references |

To regenerate (about 20 minutes; it needs the run root's `micm/` and
`tier2/scenarios/`):

    mkdir -p ~/tmp/heldout && cd ~/tmp/heldout
    zsh ~/NSF-NCAR/TurboMICM/tier2/heldout/make_heldout.sh ~/NSF-NCAR/ODE/benchmarks/tier2/scenarios/make_scenario.py heldout
    zsh ~/NSF-NCAR/ODE/benchmarks/tier2/admin/make_scenario_data.sh MICM_SCENARIO heldout/*
    python3 ~/NSF-NCAR/ODE/benchmarks/reference/make_cells_csv.py 10000 $(cat ~/NSF-NCAR/TurboMICM/tier2/heldout/SEED) --output ts1/cells.csv

The scenario files and the ranking input should match `SHA256SUMS` exactly;
MICM's results match on this machine and compiler, and the references may
differ in their last digits elsewhere.

Each entry's results on these scenarios are in its
`results/2026-10-22-tier2-<run>-verification/report.md`, beside the disclosed
ones.
