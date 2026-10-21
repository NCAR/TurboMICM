# Tier 1 follow-up: timing all nine entries

**State on 2026-10-20 at 18:50.** All nine Tier 1 runs have finished and are
verified. Eight were timed together on 2026-10-20, 12:45 to 14:24
(`results/2026-10-20-luna-xhigh-verification`, wiki Tier-1-Results and
Leaderboard). luna-high, the rerun, finished at 18:20 and was verified at 18:31
(`results/2026-10-20-luna-high-verification`: 662 PASS, 0 FAIL, 18 N/A with
both compilers; gate 5.34e-14 with clang, 1.43e-13 with gcc-15), and is
collected as `runs/luna/high` (`d3c2480`). At 18:33 the SSE deferred its timing
session, all nine entries with MICM's three programs, to this follow-on
session. The wiki says "not yet timed" wherever luna-high's time would go.

## Before starting

- No Codex run in progress, and nothing else heavy: no builds, scans, agents,
  or subagents during the session. The SSE has had an agent editing the wiki's
  slides from another session; it ran no process on this machine, but check.
- AC power (`pmset -g batt`); after a reboot the laptop may be on battery. The
  timing script waits for AC itself.
- The builds, all in `~/NSF-NCAR/tier1-work` (outside the run root):
  - each entry's clang build at `verify/<run>/clang/<run>/bin/ts1_bench_cpu`,
    for sol-high, sol-xhigh, sol-max, astra-high, astra-xhigh, astra-max,
    luna-max, luna-xhigh, and luna-high;
  - MICM's programs at `timing/bin/micm-standard-in-place-mozart` and
    `timing/bin/micm-vector8-in-place-mozart`;
  - the gate files `verify/tier1-gate/cells_2026.csv` (SHA-256 `43febf31…ad8b5`)
    and `verify/tier1-gate/micm_final_2026.csv` (`13329009…6ba7d`).

  All were present at 18:50. If any are gone: `prepare.sh` there rebuilds MICM's
  programs, the gate files, and the six Sol and Astra entries; rebuild each Luna
  entry with `admin/verify_tier3.py <run> runs/luna/<effort>
  ~/NSF-NCAR/tier1-work/verify/<run> --stages build,ts1 --compilers clang`.
  Rebuild before the session, not during it.

## The session (about 2 hours once idle)

    tmux -L turbomicm new-session -d -s tier1-timing \
      "zsh $HOME/NSF-NCAR/tier1-work/time-when-idle.sh luna-high luna-max luna-xhigh luna-high"

- The first argument is the label; the runs after it are added to the six
  Sol and Astra entries that `run-timing.sh` always times. tmux keeps the
  session alive if Claude Code exits; a Claude Code background task would not.
- `time-when-idle.sh` waits for 10 checks in a row, a minute apart, with no
  Codex session record written for 10 minutes, total CPU below one core, and
  AC power, logging to `~/NSF-NCAR/tier1-work/timing-luna-high-wait.log`. Then
  `run-timing.sh` runs `admin/time_entries.py` under `caffeinate`: 5 rotating
  rounds of all twelve programs at 10 threads and then 1 thread, each run
  checked against MICM's final state, with a one-minute load audit, into
  `~/NSF-NCAR/tier1-work/timing-luna-high/`, touching `DONE` at the end.
- The 2026-10-20 session took 99 minutes for eleven programs. luna-high adds
  five runs at each thread count; its single verification runs at 10 threads
  took 51 s (clang) and 57 s (gcc-15), and its one-thread time is unknown.
- Stay idle while it runs. To watch without loading the machine, read
  `timing-luna-high/time_entries.log` now and then.

## Afterwards

1. Check: `DONE` exists, the log ends `exit 0`, and every run says
   `gate=PASS` (120 runs).
2. Copy into `results/2026-10-20-luna-high-verification/` (its README says the
   times will be added there): `timings.csv`, `summary.md`, and
   `environment.txt` and `load-audit.log` with `/Users/fillmore` replaced by
   `~`. Regenerate the boards over the seven sessions:

       S="results/2026-10-15-micm-configurations/inplace-mozart/timings.csv results/2026-10-15-micm-configurations/others/timings.csv results/2026-10-15-astra-verification/timings.csv results/2026-10-16-astra-max-verification/timings.csv results/2026-10-19-luna-max-verification/timings.csv results/2026-10-20-luna-xhigh-verification/timings.csv results/2026-10-20-luna-high-verification/timings.csv"
       python3 admin/leaderboard.py --out DIR10 ${=S}
       python3 admin/leaderboard.py --threads 1 --out DIR1 ${=S}

   and save them as `leaderboard.md` and `leaderboard-1-thread.md`.
   `ENTRY_GROUPING` already has luna-high ("one per thread").
3. Add Timing and Findings sections to that README, after
   `results/2026-10-20-luna-xhigh-verification/README.md`: the wait, the
   procedure, the largest relative differences, the twelve-program table at 10
   and 1 threads with speedups over MICM's default and fastest, and findings:
   - luna-high against MICM's default and the other Luna entries;
   - its scaling from 1 to 10 threads (luna-xhigh's 2.79 was the least);
   - whether the other eight keep their order;
   - agreement with the 2026-10-20 session (within 2.2 % at 10 threads, 0.6 %
     at 1 thread last time);
   - its self-report as one line: its RESULTS.md gives 48.5 s against MICM's
     27.3 s at 10 threads;
   - the load audit.
4. Wiki, after `git pull --rebase` (the SSE's other agent edits the slides,
   `Overview.md`, and `tools/rank_diamonds.py`; leave those alone, and never
   force-push):
   - **Tier-1-Results:** the intro (all nine timed); every matrix cell from the
     new session, luna-high's included; a "### All nine entries, <date>"
     Verification subsection above "### luna-high, 2026-10-20", with the
     table and findings; the Run details TS1 columns for every row; luna-high's
     "not yet timed".
   - **Leaderboard:** replace both tables with the regenerated ones, then run
     `python3 tools/rank_diamonds.py Leaderboard.md` to restore the medal
     marks. Update the text: seven sessions and their list; the largest
     session-to-session differences (2.9 % entries, 3.1 % MICM last time;
     recompute); run counts; other load; the placement paragraph, dropping
     "luna-high ... joins the board after its timing session"; the one-thread
     paragraph; the data links.
   - **Home:** the status line and the speed figures.
   - **Tier-1-Implementations:** the "At a glance" TS1 column for every row,
     and luna-high's "not yet timed".
   - **Tier-3-Results:** the Luna sentence.
5. `admin/runs.md`: a "Timed with all entries" bullet under "luna-high
   (rerun)". `CLAUDE.md`: the Luna item and the timing reference (the newest
   results folder). Then delete this file and CLAUDE.md's pointer to it,
   commit, and push the repo and the wiki.

## Not done in the session that wrote this

- No times for luna-high. The verification's single runs above were not
  under the timing protocol (one run, machine not checked for idleness); do
  not publish them as timings.
- The wiki slides (the SSE's other agent) still show the state before
  luna-xhigh: slide 6 "7 of 9", slide 7 luna-max at 30.91 s, slide 16 Luna
  xhigh and high pending.
- For the analysis pass (CLAUDE.md): luna-xhigh's setup inside its timer and
  its 2.79-fold scaling; luna-high's timed loop should be read the same way.
- Codex's weekly usage window read 10 % at 18:20 (resets 2026-10-27 00:33).
- luna-high's section of Tier-1-Implementations is published (17 citations
  checked). Two findings from it are on Tier-1-Results' Notes and in the
  results folder: its case report prints PASS lines from a fixed list once
  its one C test program passes, and its generator mislabels species for CAMP
  configurations that its Markowitz order permutes (TS1 unaffected). For the
  analysis pass or the review: audit the mapping of `tests/test_core.c`'s
  checks to MICM's case ids.
