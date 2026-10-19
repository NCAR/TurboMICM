# Tier 3 follow-up: timing Fable's reference solution beside the contest's winners

**State on 2026-10-19 at 08:30.** Tier 3 is published: the contest's ranking
(`results/2026-10-19-tier3-ranking/`, wiki Tier-3-Results), the blind review,
and Fable's reference solution (`tier3/fable/max` at `cb0f8a3`, verified in
`results/2026-10-19-tier3-fable-max-verification/`, wiki
Fable-Tier-3-Reference-Solution). One measurement was deferred by the SSE:
the session that places Fable's five points beside the contest's level
winners under the ranking's protocol. It was started at 08:15 and stopped
after five runs. Everything it needs is in place.

## The session (about 1 h 45 min, idle machine, AC power)

    cd ~/NSF-NCAR/tier3-work
    tmux -L turbomicm new-session -d -s tier3-timing-fable "zsh $HOME/NSF-NCAR/tier3-work/run-timing-fable.sh"

- `run-timing-fable.sh` runs `admin/time_tier3.sh` under `caffeinate` with
  `programs-fable-vs-contest.txt` on `disclosed-10000,heldout-10000`, 3 rounds,
  10 threads, with a one-minute load audit, into `~/NSF-NCAR/tier3-work/timing-fable/`,
  and touches `DONE` at the end (363 runs).
- The programs: the contest's two fastest level-2 points (astra-xhigh/fast,
  sol-max/fast) and two fastest level-3 points (astra-xhigh/accurate,
  sol-xhigh/L3) from the 2026-10-18 re-verification builds in
  `~/NSF-NCAR/tier3-work/verify/<run>/clang/<run>/bin/`; Fable's five points
  from `~/NSF-NCAR/tier3-work/verify/fable-max-2/clang/fable-max/bin/ts1_scenario`
  (built from `cb0f8a3`); MICM's baselines from `~/NSF-NCAR/.tier3-heldout/bin/`.
  If those builds are gone, rebuild them with `admin/verify_tier3.py`
  (`--stages build`) and `tier3/heldout/make_data.sh`.
- Nothing else may run on the machine meanwhile (no builds, scans, or agents).

## Afterwards

1. `python3 admin/rank_tier3.py ~/NSF-NCAR/tier3-work/timing-fable --threads 10 > ranking.md`;
   the entries are the four contest points' owners, `fable-max`, and `micm`.
   Rescore exactly anything it flags.
2. `results/2026-10-19-tier3-fable-vs-contest/` (or the day's date): `programs.txt`,
   `ranking.md`, per scenario `timings.csv`, `summary.md`, `environment.txt`,
   `time_tier3.log`, `load-audit.log`, and a README after
   `results/2026-10-19-tier3-ranking/README.md`'s pattern: procedure, the
   level-2 and level-3 tables, margins, the load audit, and the statement that
   the contest's ranking stands as published (this session has only the
   winners' points, three rounds, and a non-contestant).
3. Wiki: replace the "Beside the contest's winners, provisionally" section of
   Fable-Tier-3-Reference-Solution with the measured table; the sentence in
   Tier-3-Results' reference-solution section; Home's paragraph.
4. `admin/runs.md`, then delete this file and the open item in CLAUDE.md.

## Other open items (unchanged)

Luna's Tier 3 arm (new held-out scenarios first), the x86-64 Linux timing,
the analysis pass on every entry, and the framework work, as listed in
CLAUDE.md.
