# Tier 1 follow-up: timing all twelve entries

**State on 2026-11-17 at 11:45 MST.** All twelve Tier 1 entries are finished
and verified. Eight were timed together on 2026-10-20, 12:45 to 14:24
(`results/2026-10-20-luna-xhigh-verification`, wiki Tier-1-Results and
Leaderboard). Four are not yet timed:

| Entry | Branch, commit | Verification | Grouping |
|---|---|---|---|
| luna-high | `runs/luna/high`, `d3c2480` | `results/2026-10-20-luna-high-verification` | one group per thread |
| terra-max | `runs/terra/max`, `9071941` | `results/2026-11-17-terra-max-verification` | 64 cells |
| terra-xhigh | `runs/terra/xhigh`, `288f274` | `results/2026-11-17-terra-xhigh-verification` | one cell |
| terra-high | `runs/terra/high`, `a8a2ea0` | `results/2026-11-17-terra-high-verification` | one cell |

The timing session for all twelve was armed on 2026-11-17 at 08:49 and
deferred by the SSE at 11:43, before it started, to this follow-on session (the
idle gate had waited three hours; see "The machine" below). The wiki says "not
yet timed" wherever these four entries' times would go.
`admin/leaderboard.py`'s `ENTRY_GROUPING` already has all twelve.

## Before starting

- No Codex run in progress, and nothing else heavy: no builds, scans, agents,
  subagents, video calls (a Google Meet call added about 1.3 cores on
  2026-11-17), or other work on this laptop during the session.
- AC power (`pmset -g batt`); after a reboot the laptop may be on battery. The
  timing script waits for AC itself.
- **The machine.** On 2026-11-17 two macOS daemons, `suggestd`
  (CoreSuggestions, about 65 % of a core) and `maild` (about 41 %, Mail.app not
  open), ran all day and kept the total near 1.1 cores, above the gate's
  one-core limit. Check first with `top -l 2 -s 5 -n 8 -o cpu -stats
  pid,command,cpu`. If they are still there, the SSE chose on 2026-11-17 to
  relax the gate to two cores (`IDLE_CPU_LIMIT=200`) and record the load;
  confirm that choice with the SSE, or restart the daemons (`killall suggestd
  maild`; launchd restarts them) with the SSE's agreement and see whether
  they settle.
- The builds, all in `~/NSF-NCAR/tier1-work` (outside the run root):
  - each entry's clang build at `verify/<run>/clang/<run>/bin/ts1_bench_cpu`,
    for sol-high, sol-xhigh, sol-max, astra-high, astra-xhigh, astra-max,
    luna-max, luna-xhigh, luna-high, terra-max, terra-xhigh, and terra-high
    (the Terra builds of 2026-11-17 used Apple clang 21.0.0, the binary dated
    2026-10-19 that built the Luna entries and MICM's programs);
  - MICM's programs at `timing/bin/micm-standard-in-place-mozart` and
    `timing/bin/micm-vector8-in-place-mozart`;
  - the gate files `verify/tier1-gate/cells_2026.csv` (SHA-256 `43febf31…ad8b5`)
    and `verify/tier1-gate/micm_final_2026.csv` (`13329009…6ba7d`).

  All were present on 2026-11-17 at 08:49, hashes checked. If any are gone:
  `prepare.sh` there rebuilds MICM's programs, the gate files, and the six Sol
  and Astra entries; rebuild each Luna or Terra entry with
  `admin/verify_tier3.py <run> runs/<model>/<effort>
  ~/NSF-NCAR/tier1-work/verify/<run> --stages build,ts1 --compilers clang`.
  Rebuild before the session, not during it.

## The session (about 2.5 to 3 hours once idle)

    tmux -L turbomicm new-session -d -s tier1-timing \
      "IDLE_CPU_LIMIT=100 zsh $HOME/NSF-NCAR/tier1-work/time-when-idle.sh all-twelve luna-max luna-xhigh luna-high terra-max terra-xhigh terra-high"

- The first argument is the label; the runs after it are added to the six Sol
  and Astra entries that `run-timing.sh` always times: 15 programs with MICM's
  three. tmux keeps the session alive if Claude Code exits; a Claude Code
  background task would not. `IDLE_CPU_LIMIT` (percent of one core, default
  100) is the gate's CPU limit; 200 if the SSE confirms it (above).
- `time-when-idle.sh` waits for 10 checks in a row, a minute apart, with no
  Codex session record written for 10 minutes, total CPU below the limit, and
  AC power, logging to `~/NSF-NCAR/tier1-work/timing-all-twelve-wait.log`.
  Then `run-timing.sh` runs `admin/time_entries.py` under `caffeinate`: 5
  rotating rounds of all 15 programs at 10 threads and then 1 thread, each run
  checked against MICM's final state, with a one-minute load audit, into
  `~/NSF-NCAR/tier1-work/timing-all-twelve/`, touching `DONE` at the end.
- Estimate: the 2026-10-20 session took 99 minutes for eleven programs. The
  four new entries' single verification runs at 10 threads took 51 s
  (luna-high), 10.0 s (terra-max), 29.6 s (terra-xhigh), and 24.4 s
  (terra-high); their one-thread times are unknown (luna-xhigh's was 60 s,
  luna-max's 183 s).
- Stay idle while it runs. To watch without loading the machine, read
  `timing-all-twelve/time_entries.log` now and then.

## Afterwards

1. Check: `DONE` exists, the log ends `exit 0`, and every run says
   `gate=PASS` (150 runs: 15 programs, 5 rounds, 2 thread counts).
2. Make `results/<date>-tier1-all-twelve-timing/` with `timings.csv`,
   `summary.md`, and `environment.txt` and `load-audit.log` with
   `/Users/fillmore` replaced by `~`, and add a Timing section to each of the
   four verification READMEs above pointing to it. Regenerate the boards over
   the eight sessions:

       S="results/2026-10-15-micm-configurations/inplace-mozart/timings.csv results/2026-10-15-micm-configurations/others/timings.csv results/2026-10-15-astra-verification/timings.csv results/2026-10-16-astra-max-verification/timings.csv results/2026-10-19-luna-max-verification/timings.csv results/2026-10-20-luna-xhigh-verification/timings.csv results/<date>-tier1-all-twelve-timing/timings.csv"
       python3 admin/leaderboard.py --out DIR10 ${=S}
       python3 admin/leaderboard.py --threads 1 --out DIR1 ${=S}

   and save them as `leaderboard.md` and `leaderboard-1-thread.md` in the new
   folder.
3. Its README, after `results/2026-10-20-luna-xhigh-verification/README.md`:
   the wait (and the gate's limit), the procedure, the largest relative
   differences, the fifteen-program table at 10 and 1 threads with speedups
   over MICM's default and fastest, and findings:
   - the Terra entries against MICM, each other, and the other nine; terra-max's
     64-cell grouping against MICM's 32-cell and default groupings;
     terra-xhigh's and terra-high's one-cell solves; terra-high's allocation
     inside its timer (`src/rosenbrock.c:167-175`);
   - luna-high against MICM's default and the other Luna entries;
   - each new entry's scaling from 1 to 10 threads (luna-xhigh's 2.79 was the
     least);
   - whether the other eight keep their order, and agreement with the
     2026-10-20 session (within 2.2 % at 10 threads, 0.6 % at 1 thread last
     time);
   - the self-reports as one line each: terra-max 9,632 ms against MICM at 64
     cells 23,835; terra-xhigh 30,496 against 29,617; terra-high 25,689
     against 26,324; luna-high 48.5 s against 27.3 s (all at 10 threads);
   - the load audit (the daemons, if present).
4. Wiki, after `git pull --rebase` (the SSE's other agents may edit the
   slides, `Overview.md`, and `tools/rank_diamonds.py`; leave those alone, and
   never force-push):
   - **Tier-1-Results:** the intro; every matrix cell from the new session,
     the four new ones included; a "### All twelve entries, <date>"
     Verification subsection above "### terra-high, 2026-11-17", with the table
     and findings; the Run details TS1 columns for every row; the four "not
     yet timed".
   - **Leaderboard:** replace both tables with the regenerated ones, then run
     `python3 tools/rank_diamonds.py Leaderboard.md` to restore the medal
     marks. Update the text: eight sessions and their list; the largest
     session-to-session differences (recompute); run counts; other load; the
     placement paragraph; the one-thread paragraph; the data links.
   - **Home:** the status line and the speed figures.
   - **Tier-1-Implementations:** the "At a glance" TS1 column for every row,
     and luna-high's "not yet timed" (the Terra entries have no sections yet).
   - **Tier-3-Results:** the Luna sentence.
   - **Entrant-<run>, the nine:** the At a glance rows "TS1, 10 threads",
     "TS1, 1 thread", and "Speedup, 1 to 10 threads" (luna-high's "not yet
     timed"), the summary's time where it gives one, and the Records line
     naming the session the times come from. Tier-1-Judging holds no times;
     Tier-1-Judging-Opus quotes one comparison (astra-xhigh 3.8 % slower than
     sol-max at ten threads, under "What the order shows"): update it if the
     gap changes.
5. `admin/runs.md`: a "Timed with all entries" bullet under "luna-high
   (rerun)" and each Terra section. `CLAUDE.md`: the Terra and Luna items and
   the timing reference (the newest results folder). Then delete this file
   and CLAUDE.md's pointer to it, commit, and push the repo and the wiki.

## Open, beyond the timing

- The Terra entries have no sections in Tier-1-Implementations and no
  `Entrant-<run>` pages (the nine have both; the template is
  `results/2026-10-20-tier1-judging-opus/brief/PAGE_TEMPLATE.md`), and no
  judging; whether Terra enters Tiers 2 and 3 is decided after Tier 1.
- terra-high's case report: each PASS id is the exit status of one of 47 test
  programs (675 ids, 86 from one). A partial audit found 24 vector-ordering LU
  and linear-solver ids passing on single-cell, standard-storage checks, and 5
  more in programs with no grouped-layout code
  (`results/2026-11-17-terra-high-verification`); the rest of the mapping
  (`tools/test_cpu.py:9-293`) is untraced. For the SSE, with the Luna gate
  findings.
- The archived transcripts (all of them, Terra's included) contain the SSE's
  Codex approved-command list, which names a work email address; the SSE has
  been told; redaction is the SSE's call.
- Codex's weekly usage window read 29 % at the end of terra-high.
