# Tier 3's held-out scenarios

Five scenarios the agents never saw, used with the six disclosed ones to check and
rank the declared operating points. The SSE approved their parameters on
2026-10-18, before any Tier 3 result. `make_heldout.sh` was committed at 01:50 that
day in a separate repository with no remote, outside the run root, and its
SHA-256, recorded in [`admin/runs.md`](../../admin/runs.md), is that of the copy
here: `08016cac9dfb41c1824648e843b2962686de94a1fcc2ddefb58a9445268e195c`. They are
published now that the last Tier 3 run, astra-max, has ended (11:14 on 2026-10-18).

All start at 06:00 local solar time with diurnal photolysis and record a snapshot
every 6 hours, varying what the disclosed README announces: seeds, latitudes,
seasons, and 48-hour runs.

| Scenario | Seed | Initial NOx | Steps | Latitude, day of year | Reference convergence: overall (worst species) | MICM's default |
|---|---:|---|---|---|---|---|
| summer-15min | 201 | x1 | 96 of 15 min (24 h) | 40° N, 172 | 10.16 (OClO 9.14) | 3.14 digits |
| tropics-20min | 202 | x1 | 72 of 20 min (24 h) | 10° N, 80 | 10.26 (OClO 8.98) | fails |
| winter60-10min | 203 | x1 | 144 of 10 min (24 h) | 60° N, 355 | 10.69 (Br 9.74) | fails |
| polluted-48h-30min | 204 | x5 | 96 of 30 min (48 h) | 30° N, 200 | 10.29 (phenol 9.36) | 1.43 digits |
| clean-48h-5min | 205 | x0.2 | 576 of 5 min (48 h) | 50° N, 260 | 10.16 (phenol 9.40) | 3.04 digits |

The references were made like the disclosed ones, with
`tier3/admin/make_references.sh` (MICM's six-stage Rosenbrock at 1e-10 relative and
1e-23 absolute tolerance, checked against 1e-9 and 1e-22), on 2026-10-18 between
sol-xhigh and sol-max, with 10 threads; the log is `references-2026-10-18.log`.

## Files

| File | Contents |
|---|---|
| `make_heldout.sh` | the five scenarios' parameters: `make_heldout.sh MAKE_SCENARIO DIR [CELLS]` |
| `make_data.sh` | makes everything the ranking needs in `~/NSF-NCAR/.tier3-heldout/`: MICM's scenario driver, the held-out scenarios at 1,000 cells with references and at 10,000 cells, and the disclosed scenarios at 10,000 cells |
| `SHA256SUMS` | checksums of the files `make_data.sh` made on 2026-10-18, relative to that folder |
| `references-2026-10-18.log` | the log of the reference runs |

To regenerate (about 30 minutes; it needs the run root's `micm/` and
`tier3/scenarios/`, and must not run while a contestant runs):

    zsh tier3/heldout/make_data.sh
    cd ~/NSF-NCAR/.tier3-heldout && shasum -a 256 -c ~/NSF-NCAR/TurboMICM/tier3/heldout/SHA256SUMS

The scenario files should match exactly; the references and their logs may differ
in the last digits, and the logs in their times, on another machine or compiler.

Each entry's results on these scenarios are in
`results/2026-10-18-tier3-<run>-verification/heldout-report.md`.
