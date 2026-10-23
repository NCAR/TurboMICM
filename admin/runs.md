# Run log

Launch record for the MICM-to-C runs. Times are local: MDT for Tier 1 and Tier 3,
EDT for Tier 2, which runs on another machine. Tier 1 comes first; the Tier 2 runs
follow, under [Tier 2](#tier-2), and the Tier 3 runs, under [Tier 3](#tier-3).

Each run is an interactive session started in its run directory,
`benchmarks/<run>/`, and launched with one instruction:

    /goal read PROMPT.md and work until done, commit, push

The run repositories have no remote, so "push" has nowhere to go.

A run is `<model>-<effort>`; repeat instances of the same model and effort, run as
time permits, add a suffix: `sol-high-2`, `sol-high-3`. Each finished run is kept
on the TurboMICM branch `runs/<model>/<effort>`, or `runs/<model>/<effort>-<n>`
for a repeat.

Contestant changes
- 2026-10-16: the contestants are OpenAI's models only (Luna, Sol, Astra), for
  practical reasons given usage plans. Anthropic's Claude Opus and Claude Fable,
  at max reasoning effort, design the experiment, orchestrate the runs, and
  validate, evaluate, and judge the entries, and will later build on the results
  a framework for MUSICA, CheMPAS-A, and CAM-chem. The sonnet-*, opus-*, and
  fable-* run directories in the run root were never started and are unused.
- 2026-10-15, before any Terra run: Terra was replaced by Luna. Terra was available
  only as gpt-5.6-terra, which Codex lists as an "Older balanced model for
  straightforward work", while every other contestant is its family's current
  release. Luna is gpt-6-luna, "Fast and affordable model for easier tasks", with
  high, xhigh, and max effort levels. The runs are luna-high, luna-xhigh, and
  luna-max; the prompts differ from Terra's only in the name. The Claude side
  (Sonnet, Opus, Fable) has no fast-tier counterpart; Claude Haiku would be the
  match if one is added.

Prompt changes
- 2026-10-15, after the three Sol runs finished: RULES gained "and do not look at
  other solutions to this task, on this machine or online", because finished
  solutions now live in ~/NSF-NCAR/TurboMICM and on GitHub. Every prompt in
  prompts/ carries it. The Sol runs received the prompt without it, which the
  experimenter judged immaterial to their solutions.

Launch changes
- 2026-10-16, from luna-high on: at the SSE's request, Claude Code launches the
  runs with `admin/launch_codex_run.sh`. It is the same interactive Codex session,
  command, and goal, typed into a detached tmux session started from an empty
  environment with the SSE's login shell, and checked on screen and in the session
  record. The method was first tested on a throwaway repository. The Sol and
  Astra runs were launched by the SSE.
- 2026-10-19, during luna-max: the SSE authorized Claude Code to stop a Luna run
  it judges to be floundering. There is no time limit; slowness alone is no
  reason. Claude Code checks each run every 20 minutes (commits, the agent's
  failing-case counts, tokens) with a read-only probe of the session record,
  and records any stop here with its reasons.

Reference changes
- 2026-10-15: `reference/micm_ts1_ref.cpp` gained `--group-size G`, so MICM can be
  timed with the same cell grouping as an entry. Its default, one group per
  thread, is unchanged and gives bit-identical results.

Harness changes
- 2026-10-15, between astra-high and astra-xhigh: Codex CLI was upgraded from
  0.162.1 to 0.163.0 (Homebrew cask, 21:14), and the experimenter chose to keep
  it. The Sol runs and astra-high ran on 0.162.1; astra-xhigh and later runs use
  0.163.0. Codex's instructions to gpt-6-astra are the same text in both versions.

Verification
- 2026-10-15: the three Sol entries were rebuilt from their branches with clang
  and gcc-15, their tests rerun and counted against the case list, their builds
  audited against the portability rules, and their TS1 timings replicated with
  `admin/time_entries.py` alongside MICM. The results are in
  `results/2026-10-15-sol-verification/` and on the wiki's Tier-1-Results page.
- 2026-10-15 and 16: astra-xhigh was verified the same way (astra-high's tests
  had been verified on 2026-10-15), and all five entries were timed in one
  session, from 23:49 to 00:55, with MICM's default, the default with 32-cell
  groups, and MICM's fastest configuration, at 10 and 1 threads. All 80 runs
  passed the gate. The results are in `results/2026-10-15-astra-verification/`.
- 2026-10-16: astra-max was verified the same way, and all six entries were
  timed in one session, from 02:58 to 04:09, with the same three MICM programs.
  All 90 runs passed the gate. The results are in
  `results/2026-10-16-astra-max-verification/`.
- 2026-10-16: documenting the implementations showed that astra-high's benchmark
  recomputes every rate constant before each timed step (`src/ts1_bench.c:154`),
  where MICM's driver and the other entries compute them once, before the
  warm-up. Its times include that extra work; the earlier verification had not
  noted it.

Transcripts
- Each finished run's session transcript is archived on the branch `transcripts`,
  in `transcripts/<run>.jsonl.gz`, kept apart from `main` and the run branches
  because a transcript names the model and shows every step. The files are made
  with `admin/archive_transcript.py`, which redacts the account identifiers in
  the session header and gzips the rest unchanged.
- 2026-10-15: sol-high, sol-xhigh, and sol-max archived. Their reasoning is
  encrypted by OpenAI and cannot be read; messages, tool calls with their full
  output, and token counts are complete.
- 2026-10-15: astra-high archived. A scan for credentials found only chance
  matches inside the encrypted reasoning.
- 2026-10-15: astra-xhigh archived. The credential scan's only matches (8, of the
  `sk-` pattern) were inside encrypted reasoning.
- 2026-10-16: astra-max archived. The scan's only matches (6, of the same
  pattern) were inside encrypted reasoning.

| Run | Model | Effort | Harness | Session | Final commit | Status |
|---|---|---|---|---|---|---|
| sol-high | gpt-6.1-sol | high | Codex CLI 0.162.1, interactive | 2026-10-14 21:55-23:32 | `1d944f4`, 23:32 (6 commits) | finished |
| sol-xhigh | gpt-6.1-sol | xhigh | Codex CLI 0.162.1, interactive | 2026-10-15 03:54-05:29 | `171f54a`, 05:29 (6 commits) | finished |
| sol-max | gpt-6.1-sol | max | Codex CLI 0.162.1, interactive | 2026-10-15 07:26-09:13 | `001fb64`, 09:13 (7 commits) | finished |
| astra-high | gpt-6-astra | high | Codex CLI 0.162.1, interactive | 2026-10-15 12:10-14:26 | `4055d7c`, 14:26 (15 commits) | finished |
| astra-xhigh | gpt-6-astra | xhigh | Codex CLI 0.163.0, interactive | 2026-10-15 21:17-23:35 | `fbd365c`, 23:35 (11 commits) | finished |
| astra-max | gpt-6-astra | max | Codex CLI 0.163.0, interactive | 2026-10-16 01:04-02:46 | `0dc7f67`, 02:46 (7 commits) | finished |
| luna-high | gpt-6-luna | high | Codex CLI 0.163.0, interactive, launched by Claude Code | 2026-10-16 06:47-08:11 | `217bf6e`, 08:10 (23 commits) | stopped; to rerun |
| luna-max | gpt-6-luna | max | Codex CLI 0.163.0, interactive, launched by Claude Code | 2026-10-19 19:14-21:05 | `da2d722`, 21:05 (5 commits) | finished |
| luna-xhigh | gpt-6-luna | xhigh | Codex CLI 0.163.0, interactive, launched by Claude Code | 2026-10-20 09:31-12:20 | `ada55d5`, 12:20 (28 commits) | finished |
| luna-high | gpt-6-luna | high | Codex CLI 0.163.0, interactive, launched by Claude Code (rerun) | 2026-10-20 14:27-18:20 | `d3c2480`, 18:20 (128 commits) | finished |

Token counts are as Codex reports them: total = uncached input + output.

## sol-high

- Session start: `codex -c features.memories=false` in `benchmarks/sol-high`; model
  and effort from `~/.codex/config.toml`.
- Told during the session not to push; replied "Understood. I'll commit locally
  and skip the push." No remote was added.
- Tokens: 1,399,332 total; input 1,105,996 (plus 56,121,856 cached); output
  293,336 (97,708 reasoning).
- Codex session `01a13db3-1129-76d1-8797-c3b73194043b` (`codex resume` with that id).
- Final commit `1d944f4`, "Complete CPU case ports and verify both compiler gates".
  README.md, DESIGN.md, DECISIONS.md, RESULTS.md, environment.yml, and Makefile
  are present; the working tree is clean.

## sol-xhigh

- Session start: `codex -c features.memories=false` in `benchmarks/sol-xhigh`. The
  default effort in `~/.codex/config.toml` was changed to xhigh at 03:55, and the
  session record shows xhigh for every turn.
- Tokens: 1,095,606 total; input 859,920 (plus 43,877,504 cached); output 235,686
  (61,432 reasoning).
- Codex session `01a13efc-1abd-7381-a20e-142540e76a25`.
- Final commit `171f54a`, "Validate portable CPU port, select configurable grouping
  and document completed implementation". No remote; the working tree is clean;
  all deliverables are present.
- Self-reported, not yet verified: 675 PASS, 0 FAIL, 5 N/A (C++-only) with both
  compilers; TS1 gate passes; clang benchmark 7.71 s against MICM's 26.48 s (about
  3.4 times faster); about 1 h 33 min of work.

## sol-max

- Session start: `codex -c features.memories=false` in `benchmarks/sol-max`, with
  the default effort in `~/.codex/config.toml` set to max; the session record shows
  max for every turn.
- Tokens: 1,086,171 total; input 833,061 (plus 43,570,560 cached); output 253,110
  (79,758 reasoning).
- Codex session `01a13fbc-d221-7b81-88a8-3d73454dd2b7`.
- Final commit `001fb64`, "Record final compiler validation, accuracy gates, and
  benchmark results". No remote; the working tree is clean; all deliverables are
  present.
- Self-reported, not yet verified: 675 PASS, 5 N/A, no failures with both
  compilers; all TS1 accuracy gates pass; 10,000-cell clang median 6.94 s against
  MICM's 26.01 s (3.75 times faster); about 1 h 47 min of work.

## astra-high

- Session start: 2026-10-15 12:10, `codex -c features.memories=false` in
  `benchmarks/astra-high`. The session opened on gpt-6.1-sol, still the model in
  `~/.codex/config.toml`, and was switched to gpt-6-astra at high effort before
  the goal was given. Its record therefore carries Sol's base instructions, with
  Astra's in a model-switch message; every turn ran on gpt-6-astra at high
  effort.
- Goal: "read PROMPT.md and work until done, commit", with no push.
- Codex session `01a140c1-b679-7082-aef7-62cd2c85941c`.
- Tokens: 1,429,356 total; input 1,148,542 (plus 42,448,384 cached); output
  280,814 (41,030 reasoning).
- Final commit `4055d7c`, 14:26, "Validate complete CPU port and final TS1
  acceptance gates". No remote; the working tree is clean; all deliverables are
  present, plus a WORK_LOG.md.
- Self-reported, not yet verified:
  - all 680 cases pass with both compilers, none N/A: it wrote C equivalents of
    the five C++-only cases;
  - the TS1 gate passes;
  - clang benchmark 9.63 s against MICM's 27.42 s (2.85 times faster);
  - about 2 h 16 min of work.
- Other load during the session: preparation for Tier 2 (since renumbered Tier 3)
  ran from 13:24 to 13:34, using
  at most two threads at nice 10–15. It included two compiles, a few seconds of
  test runs, and the start of the reference computation. Astra ran no timing in
  that window; it was porting the sparse-matrix tests. At 13:34 the reference
  computation was suspended (SIGSTOP) until the run finished, to keep the machine
  idle. No other work from the experiment ran during the Sol sessions.

## astra-xhigh

- Session start: 2026-10-15 21:17,
  `codex -m gpt-6-astra -c model_reasoning_effort=xhigh -c features.memories=false`
  in `benchmarks/astra-xhigh`, on Codex CLI 0.163.0; the session record shows
  gpt-6-astra at xhigh.
- Goal, 21:17: "read PROMPT.md and work until done, commit".
- Codex session `01a142b6-b159-7d52-801a-7ac70fa3fed4`.
- Before it, `codex --yolo` was opened in the directory at 21:14, with memories on,
  and closed before any message; it left no session record.
- The whole run was one turn, from the goal at 21:17:29 to its end at 23:35:50;
  the session's context was compacted four times along the way, at 21:52, 22:21,
  22:50, and 23:17 (the transcript's `compacted` records; this log said twice
  until 2026-10-21). Every turn context
  shows gpt-6-astra at xhigh, and no model switch.
- Tokens: 1,879,105 total; input 1,607,225 (plus 31,004,544 cached); output
  271,880 (48,132 reasoning).
- Final commit `fbd365c`, 23:35, "Validate complete CPU port and TS1 gate; record
  portable benchmark results". No remote; the working tree is clean; all
  deliverables are present, plus WORK.md and LICENSE.
- Self-reported, not yet verified: all 680 cases pass with both compilers; the
  10,000-cell TS1 gate passes with both; clang median 7.31 s against MICM's
  28.08 s at ten threads (3.84 times faster); about 2 h 18 min of work.
- Other load during the session: Claude Code edited the wiki and run log, read
  files, and once ran a shell for a few seconds to test a launch environment. No
  builds, tests, or timing ran.
- After the run, from 23:46 to 23:50, a race between parallel builds in Claude's
  verification script left two self-referencing links inside the run root,
  `micm/micm` and `reference/reference`. They were removed at 23:50, before any
  other run started, and the script now makes its links without following them.

## astra-max

- Session start: 2026-10-16 01:04,
  `codex -m gpt-6-astra -c model_reasoning_effort=max -c features.memories=false`
  in `benchmarks/astra-max`, on Codex CLI 0.163.0; the session record shows
  gpt-6-astra at max, base instructions from gpt-6-astra, and no model switch.
- Goal, 01:05: "read PROMPT.md and work until done, commit".
- Codex session `01a14386-84c0-70a3-867a-6cf69b22b824`.
- The whole run was one turn, from the goal at 01:05:31 to its end at 02:46:29;
  the session's context was compacted three times. Every turn context shows
  gpt-6-astra at max.
- Tokens: 1,022,510 total; input 831,287 (plus 16,987,776 cached); output 191,223
  (31,930 reasoning).
- Final commit `0dc7f67`, 02:46, "Record complete compiler validation and TS1
  performance results". No remote; the working tree is clean; all deliverables
  are present.
- Self-reported, not yet verified: 680 of 680 cases pass with both compilers;
  the TS1 gates pass; clang median 8.24 s, 3.25 times faster than MICM on all
  cores; about 1 h 41 min of work.
- Other load during the session: Claude Code made two brief progress checks and
  edited the run log and wiki. No builds, tests, or timing ran.

## luna-high

- Session start: 2026-10-16 06:47:54, launched by Claude Code with
  `admin/launch_codex_run.sh luna-high gpt-6-luna high`, in a detached tmux
  session (socket `turbomicm`, session `luna-high`); the goal was sent at
  06:47:58. The session record shows gpt-6-luna at high, base instructions from
  gpt-6-luna, Codex 0.163.0, and no model switch; the process's environment holds
  no Claude Code variables besides the SSE's own `CLAUDE_CODE_DISABLE_ALTERNATE_SCREEN`.
- Codex session `01a144c1-3ac3-7f31-b3de-4cf81d130c81`.
- Stopped at 08:11, after 1 h 23 min, by the SSE's decision: the Claude Code
  session had to end early, and luna-high was far from done. Its first turn
  ended at 07:51 with "The task remains incomplete: both compiler test suites
  report 615 failures, and the 10,000-cell seed 2026 benchmark still has 2,185
  out-of-tolerance values"; the goal had started a second turn. Claude Code
  interrupted the goal (Escape, then `/quit`) and closed the tmux session.
- Tokens at the stop: 1,677,749 total; input 1,247,782 (plus 136,146,944 cached);
  output 429,967 (212,664 reasoning).
- The stopped attempt, 23 commits ending `217bf6e` (08:10) plus 5 uncommitted
  paths, was moved intact to `~/NSF-NCAR/aborted-runs/luna-high-2026-10-16`,
  outside the run root, and `luna-high/` was recreated from its prompt with
  `admin/setup_runs.py` for a fresh run, overnight.

## luna-max

- Order: on 2026-10-19 the SSE chose to run Luna's three entries fresh, one at a
  time, at max, then xhigh, then high (luna-high's rerun last). Each is collected,
  verified, and timed with all entries before the next launches.
- Session start: 2026-10-19 19:14:52, launched by Claude Code with
  `admin/launch_codex_run.sh luna-max gpt-6-luna max`, in a detached tmux
  session (socket `turbomicm`, session `luna-max`); the goal was sent at
  19:14:56. The session record shows gpt-6-luna at max, base instructions from
  gpt-6-luna, Codex 0.163.0, and no model switch; the process's environment holds
  no Claude Code variables besides the SSE's own `CLAUDE_CODE_DISABLE_ALTERNATE_SCREEN`.
- Codex session `01a156e0-2dbb-7742-9165-1f8fc55480ba`.
- Progress, from Claude Code's checks every 20 minutes: 61 cases passing at 18
  minutes, 446 at 31, 528 at 44, 67 failing at 61; the TS1 gate's mismatch traced
  to MICM's clamp of negative concentrations after each interval, at 80. No sign
  of floundering at any check.
- Session end: the first turn ended at 20:43 with DONE: 671 PASS, 0 FAIL, 9 N/A
  with both compilers, the TS1 gate passing (largest relative difference
  2.78e-13), and its code slower than MICM in its own single runs (185,107 ms
  against 145,414 ms on one thread, 34,479 against 31,199 ms on eight). The goal
  started a second turn, which replaced shared smoke tests with case-specific
  ports of MICM's analytical tests and switched the LU updates to `fma`; it ended
  at 21:05:47 with DONE, and Codex showed "Goal achieved (1h 50m)". The record was
  idle 10 minutes by 21:16; Claude Code quit Codex (`/quit`) and closed the tmux
  session.
- Final commit `da2d722` (21:05), 5 commits; collected as `runs/luna/max`.
- Tokens: 2,918,597 total; input 2,262,991 (plus 122,233,984 cached); output
  655,606 (354,791 reasoning). Codex's weekly usage window went from 69 % to 75 %.
- Other load during the session: from 20:14 the SSE ran a separate Codex session
  (gpt-6-astra, xhigh) in another project, active through the end of luna-max's;
  Claude Code's read-only progress checks. No builds, tests, or timing by the
  experimenter.
- Transcript archived as `transcripts/luna-max.jsonl.gz`. The credential scan's
  only matches (28, of the `sk-` pattern) were inside encrypted reasoning and
  encrypted compaction history.
- Verified, 21:22, with `admin/verify_tier3.py --stages build,tests,ts1,audit`
  (`results/2026-10-19-luna-max-verification`): 671 PASS, 0 FAIL, 9 N/A with
  each compiler, every id once; the TS1 gate passes (2.78e-13 with clang, 2.2e-13
  with gcc-15); the audit finds nothing outside the rules. The library builds
  without warnings with clang and with three `-Wint-in-bool-context` warnings
  with gcc-15. The benchmark keeps MICM's default grouping, one group per thread,
  and does MICM's default solver work.
- Timed with all entries, 23:06 to 00:36, once the machine had been idle for 10
  minutes (`~/NSF-NCAR/tier1-work/time-when-idle.sh`; the wait, from 21:25, ran
  past the SSE's separate Codex session, which last wrote at 21:34, and an hour of
  Spotlight indexing that ended at 22:57). All 100 runs passed the gate. luna-max:
  30,914 ms at 10 threads and 183,244 ms on one, 0.90 and 0.79 times MICM's
  default (27,876 and 144,522 ms); the other six keep their order, their medians
  within 1.3 % of the 2026-10-16 session. Wiki updated: Tier 1 Results,
  Leaderboard, Home, Overview, Tier 3 Results, Tier 1 Implementations.

## luna-xhigh

- Wait: Claude Code's launcher, armed at 00:46 on 2026-10-20 to start luna-xhigh
  once no Codex session was running, waited on the SSE's separate Codex session,
  open and idle since 21:34; the laptop restarted at 09:28, ending both, before
  it launched.
- Session start: 2026-10-20 09:31:08, launched by Claude Code with
  `admin/launch_codex_run.sh luna-xhigh gpt-6-luna xhigh`, in a detached tmux
  session (socket `turbomicm`, session `luna-xhigh`); the goal was sent at
  09:31:13. The session record shows gpt-6-luna at xhigh, base instructions from
  gpt-6-luna, Codex 0.163.0, and no model switch; the process's environment holds
  no Claude Code variables besides the SSE's own `CLAUDE_CODE_DISABLE_ALTERNATE_SCREEN`.
- Codex session `01a159f0-2065-7583-b7c2-91934e443647`. Codex's weekly usage
  window stood at 9 % at the start (it had reset since luna-max; next reset
  2026-10-27 00:33).
- The SSE repeated, at 09:34, the authorization to stop a floundering Luna run.
- Progress, from Claude Code's checks every 20 minutes: the 10,000-cell TS1
  comparison passing at 38 minutes (10:09), with 43 cases passing at 41; 55 at
  61; 95 at 81, when the first turn ended (10:49) reporting the acceptance gate
  not met; 172 at 101; 245 at 121; 580 at 141, after dense and sparse matrix
  storage; 663 passing and none failing at 159. The count rose at every check;
  no sign of floundering.
- Session end: the second turn ended at 12:16:46 with DONE: 662 PASS, 0 FAIL, 18
  N/A with both compilers, TS1 serial and four-thread checks passing. The goal
  started a third, which found that its earlier `CC=gcc` runs had used Apple
  clang, reran the cases with `gcc-15`, refreshed its TS1 comparison (largest
  relative difference 1.77e-13; 21,397 ms against MICM's 27,851 ms in one
  10-thread run), and ended at 12:20:38 with DONE. The record was idle 10
  minutes by 12:31; Claude Code quit Codex (`/quit`) and closed the tmux
  session.
- Final commit `ada55d5` (12:20), 28 commits; collected as `runs/luna/xhigh`.
- Tokens: 4,278,572 total; input 3,223,268 (plus 227,209,088 cached); output
  1,055,304 (563,811 reasoning). Codex's weekly usage window read 9 % at the
  start and at the end.
- Other load during the session: Claude Code's read-only progress checks; the
  SSE noted at 10:28 that an agent in another session was working on the wiki's
  slides (no such process ran on this machine at that check). No builds, tests,
  or timing by the experimenter.
- Transcript archived as `transcripts/luna-xhigh.jsonl.gz`. The credential scan's
  only matches (53 of the `sk-` pattern, one of `ghp_`) were inside encrypted
  content.
- Verified, 12:32, with `admin/verify_tier3.py --stages build,tests,ts1,audit`
  (`results/2026-10-20-luna-xhigh-verification`): 662 PASS, 0 FAIL, 18 N/A with
  each compiler, every id once (11 of the N/A cases N/A in no other entry); the
  TS1 gate passes (1.77e-13 with clang, 2.33e-13 with gcc-15); the audit finds
  nothing outside the rules (the banned-construct scan's three matches are a
  test variable named `vector_size`); no warnings with either compiler. The
  benchmark keeps MICM's default grouping, one group per thread, and does
  MICM's default solver work.
- Timed with all entries, 12:45 to 14:24, once the machine had been idle for 10
  minutes (`~/NSF-NCAR/tier1-work/time-when-idle.sh`, waiting from 12:34). All
  110 runs passed the gate. luna-xhigh: 21,443 ms at 10 threads and 59,916 ms
  on one, 1.32 and 2.41 times MICM's default (28,198 and 144,435 ms), scaling
  2.79 times from 1 to 10 threads, the least of any program; the other seven
  keep their order, their medians within 2.2 % of the 2026-10-19 session.

## luna-high (rerun)

- The fresh run of luna-high, the last of Luna's three, launched once
  luna-xhigh was verified and timed. The first attempt is under "luna-high"
  above.
- Session start: 2026-10-20 14:27:39, launched by Claude Code with
  `admin/launch_codex_run.sh luna-high gpt-6-luna high`, in a detached tmux
  session (socket `turbomicm`, session `luna-high`); the goal was sent at
  14:28:11. The session record shows gpt-6-luna at high, base instructions from
  gpt-6-luna, Codex 0.163.0, and no model switch; the process's environment holds
  no Claude Code variables besides the SSE's own `CLAUDE_CODE_DISABLE_ALTERNATE_SCREEN`.
- Codex session `01a15aff-97c7-7251-a51e-5e14f2fa689d`. Codex's weekly usage
  window read 9 % at the start.
- Progress, from Claude Code's checks every 20 minutes: the first turn ended at
  14:35, seven minutes in, with 2 cases passing and "not complete"; a one-cell
  TS1 mismatch at 20 minutes; the 10,000-cell gate failing with 2,185 of
  2,110,000 values out of tolerance at 40 minutes (15:03), the first attempt's
  count exactly. Five full comparisons from 15:03 to 15:38, each after a
  different change (species order, pressure, falloff rates, shared step
  control), gave the same 2,185; 44 cases passed at 81 minutes, and three turns
  had ended "not complete".
- Decision, 15:49, not to stop it: the case count was rising (2 to 44 in 40
  minutes) and each attempt at the gate was a new hypothesis, not a repeat, so
  by the standard Claude Code had set (a count that stalls or falls while the
  agent circles the same ground) it was slow, not floundering; the unchanged
  2,185 was the warning sign, to be weighed again at the next checks. At 16:02
  the agent clamped the solver states, and the gate passed at 16:09 (largest
  relative difference 5.34e-14). Cases then climbed: 89 at 16:17, 171 at 16:37,
  340 at 16:57, 551 at 17:19, 621 at 17:54, and 662 passing, none failing, 18
  N/A at 18:09.
- Session end: the tenth turn ended at 18:09:56 with DONE; the goal started an
  eleventh, which audited the work against the prompt, sped up the sparse
  Jacobian construction (its own 10,000-cell time from about 298 s to 72.8 s,
  then 48.5 s, against MICM's 27.3 s), and ended at 18:20:21 with DONE. The
  record was idle 10 minutes by 18:30; Claude Code quit Codex (`/quit`) and
  closed the tmux session.
- Final commit `d3c2480` (18:20), 128 commits; collected as `runs/luna/high`.
- Tokens: 4,475,652 total; input 3,325,961 (plus 372,281,344 cached); output
  1,149,691 (533,305 reasoning). Codex's weekly usage window read 10 % at the
  end.
- Other load during the session: Claude Code's read-only progress checks, its
  writing of luna-xhigh's records and wiki pages, and, in the first minutes, a
  read-only subagent drafting luna-xhigh's section of Tier 1 Implementations.
  No builds, tests, or timing by the experimenter.
- Transcript archived as `transcripts/luna-high.jsonl.gz`. The credential scan's
  only matches (55 of the `sk-` pattern, one of `gho_`) were inside encrypted
  content.
- Verified, 18:31, with `admin/verify_tier3.py --stages build,tests,ts1,audit`
  (`results/2026-10-20-luna-high-verification`): 662 PASS, 0 FAIL, 18 N/A with
  each compiler, every id once (10 of the N/A cases N/A in no other entry); the
  TS1 gate passes (5.34e-14 with clang, 1.43e-13 with gcc-15); the audit finds
  nothing outside the rules (the banned-construct scan's 75 matches are all a
  parameter or variable named `vector_size`); no warnings with either compiler.
  One group per thread, with MICM's default solver work at ten threads.
- Not yet timed: at 18:33 the SSE deferred the timing session, all nine entries
  with MICM's three programs, to a follow-on session.

# Tier 2

Launch record for the Tier 2 runs (specialization for a frozen mechanism; see
`tier2/README.md`). Times are EDT, the clock of the machine Tier 2 runs on.

Each run works in `benchmarks/tier2/<run>/`, a copy of its own Tier 1 repository
at the Tier 1 run's final commit, made by `tier2/admin/setup_runs.py` from the
TurboMICM branch `runs/<model>/<effort>`, which holds the Tier 1 repository
unchanged: branch `main`, the Tier 1 history, no remote, no tags, and PROMPT.md
holding the Tier 2 prompt. Claude Code launches every run with

    admin/launch_codex_run.sh --tier2 <run> <model> <effort>

which types `codex -m <model> -c model_reasoning_effort=<effort> -c features.memories=false
-c approval_policy=never -c sandbox_mode=danger-full-access` (no web search, as in
Tier 1; the last two settings give the full access every contest run has had,
which this machine's `config.toml` does not) into a detached tmux session (socket
`turbomicm`, session `tier2-<run>`), and then the goal:

    /goal read PROMPT.md and work until done, commit

Order
- 2026-10-22, 09:00: the SSE: "orchestrate the Tier 2 experiments. all
  contestants, in reverse order of capability, that is start with Astra Max,
  then Astra Extra High, ..., followed by Sol, Luna, ... continuously update the
  wiki as results come in." So all nine Tier 1 entrants, one run each, one at a
  time: astra-max, astra-xhigh, astra-high, sol-max, sol-xhigh, sol-high,
  luna-max, luna-xhigh, luna-high. Each is collected, archived, verified, and
  timed before the next launches.

The machine
- Tier 2 runs on a second laptop: an Apple M5 Pro (LALAL0726071122; 6 Super and
  12 Performance cores, 24 GB; macOS 26.7.1; Command Line Tools 27, no Xcode;
  Homebrew in `~/homebrew`, installed without administrator rights; user
  `dfillmor`). Tier 1 and Tier 3 ran on the SSE's other Mac, the Apple M5
  laptop Gemini.local (4 Super and 6 Efficiency cores, 16 GB, macOS 27.0.1,
  Homebrew in `/opt/homebrew`, user `fillmore`). The SSE, 2026-10-22: "note the
  machine change, Tier 1 was run on my other Mac".
- Setup, 2026-10-22 from 09:00, the machine having restarted at 08:49:
  - tmux 3.8 (Homebrew bottle) and GCC 15.3.0 (`gcc@15`, built from source,
    since Homebrew's bottles need the `/opt/homebrew` prefix); libomp 23.1.3,
    cmake 4.4.4, and conda 26.5.3 were already installed. Codex CLI is 0.163.0,
    as for every run since astra-xhigh.
  - MICM e15898ce cloned into `benchmarks/micm` and built with Apple clang
    (`SDKROOT` set), Release, double precision, OpenMP through Homebrew's libomp
    (MICM's CMake finds it with `brew --prefix libomp`), and
    `FETCHCONTENT_TRY_FIND_PACKAGE_MODE=NEVER`: 67 of 67 ctest tests pass. The
    case list regenerated from this build is identical to
    `reference/micm_cases.txt` (680 cases).
  - MICM's TS1 driver, built here, reproduces
    `reference/micm_ts1_final_100_seed1.csv` byte for byte, and its 10,000-cell
    seed-2026 final state has the old machine's SHA-256 (`13329009…6ba7d`).
  - The five disclosed day-night scenarios, regenerated here, match Tier 3's
    checksums (`tier3/heldout/SHA256SUMS`, the 10,000-cell versions) byte for
    byte, and MICM's runs at the gate's settings take the 2026-10-16
    calibration's internal steps exactly.
  - Apple clang links here with or without `SDKROOT`.
- Every Tier 1 Makefile names `/opt/homebrew/opt/libomp`, which does not exist
  here, so each agent must repoint its Makefile; the prompt says so.

Prompt changes
- `tier2/prompt.v2` (2026-10-22), from `prompt.v1` (2026-10-16, never used):
  - MACHINE AND TOOLCHAIN describes this machine, and says that it is not
    Tier 1's, that Homebrew's paths differ, and that
    `reference/README.md`'s timings and `SDKROOT` note describe Tier 1's machine;
  - the honor rule names Tier 1: "do not look at any other agent's solution, to
    this tier or to Tier 1, on this machine or online" (v1: "other solutions to
    this task"), since every Tier 1 solution is on this machine and on GitHub,
    and two of them already generate mechanism-specific code;
  - `make test-specialized` runs the generator's tests, so the gate's check of
    them is mechanical, and it joins the targets that must work with both
    compilers;
  - the timing and RESULTS.md use this machine's 18 cores (v1: 10).

Undisclosed inputs, fixed before any Tier 2 run, in `~/NSF-NCAR/.tier2-heldout/`
(a git repository with no remote, outside the run root and this repository)
- 2026-10-22, 09:30:00 (commit `eb4ae5f`): `make_heldout.sh`, five held-out
  day-night scenarios. SHA-256
  `c8dfef98d54b09c52910a764cffd7791836278323f73c7b9f6e325543365750b`. Tier 3's
  held-out set is public, so Tier 2 has its own; the scenarios vary what the
  agents' scenarios README announces (seeds, latitudes, seasons, start times,
  and initial NOx; steps of 5 to 30 minutes; runs of up to 48 hours). Claude
  Code chose them, under the SSE's instruction to orchestrate Tier 2; the
  script is published when the last Tier 2 run ends.
- 2026-10-22, 09:30:31 (commit `e10fd84`): the TS1 ranking input, 10,000 cells
  from `reference/make_cells_csv.py` with a 63-bit seed from Python's
  `secrets`. The cells file's SHA-256 is
  `be7ccfdbc8c74a68c1112122e14bb2407009eba50d84485a784af8597e263772`, and MICM's
  final state for it `cb7d0cc7014b598515a70805364e3854d5f318ccc1466ebd90e01b1b5bac62ce`.

Launch changes
- 2026-10-22: `admin/launch_codex_run.sh` gained `--tier2`: the checks of
  `--tier3` (an untouched copy at the Tier 1 commit the prompt names, no remote,
  nothing changed but PROMPT.md) plus no tags; no `--search`; tmux session
  `tier2-<run>`. For this machine it finds tmux in `~/homebrew`, waits for the
  login shell's prompt before typing (a smart-card PIN prompt from the SSE's
  `.zshrc`, which ends with `ssh-add -s`, would get Ctrl-C, never a PIN; with
  the launcher's empty environment `ssh-add` finds no agent and exits at once),
  and waits up to 15 s for Codex's status line, which can be drawn after the
  composer.
- 2026-10-22, 09:21: the Tier 2 mode was tested on a throwaway copy in Claude
  Code's scratch space, with gpt-6.1-sol at low effort and a one-line task
  (sessions `01a163c5-b393-7fc3-b6e6-700501dafde1`, which the status-line race
  stopped before its goal, and `01a163c6-4b55-7673-bcef-dd66501794c6`). The
  second launched and verified, but ran with approval `on-request` and a
  `workspace-write` sandbox without network, and stopped to ask leave to `git
  commit`: this machine's `config.toml` does not grant the full access (approval
  `never`, sandbox `danger-full-access`) that every Tier 1 and Tier 3 run had,
  from the old machine's `config.toml`, as their transcripts' turn contexts
  show. The request was declined and the session closed. Claude Code's
  auto-mode safety check refused to let Claude Code add those two settings to
  the launch command, and the SSE decided.
- 2026-10-22, about 09:40, the SSE: "I have a Codex usage reset available, so
  now is a good time to run these. go with option 1." (the same conditions as
  Tiers 1 and 3), and "approved to pass -c approval_policy=never -c
  sandbox_mode=danger-full-access". The launcher now passes both, for every
  tier, and its check of the session record requires approval `never` and
  sandbox `danger-full-access` in every turn context, waiting for the first.
- 2026-10-22, about 09:41: the SSE had started a Codex session of its own at
  about 09:34 (`codex --yolo`, pid 15599, in `~/NSF-NCAR/CERES-MIP`), which the
  launcher's check of other Codex sessions stops at. Asked whether to launch
  alongside it or wait for it to end, the SSE chose to launch alongside it. The
  launcher's new `ALLOW_CODEX_PIDS` names the sessions allowed to run alongside,
  and prints each one's command and directory; each run's record lists them as
  other load, as luna-max's did.
- 2026-10-22, 09:41:59: the throwaway test, rerun with both settings (session
  `01a163d8-e198-7023-b5c8-1baa718aee9a`), verified: gpt-6.1-sol at low,
  approval `never`, sandbox `danger-full-access`, no model switch, the goal. It
  committed `hello.txt` and ended with DONE 18 s after the goal; Claude Code
  quit it.

Calibration on this machine
- 2026-10-22, 09:43–09:45: MICM's scenario driver built with g++-15 (GCC 15.3.0)
  and run on the five disclosed scenarios at the gate's settings takes exactly
  the clang build's internal steps on each and agrees with it to 14.37, 14.28,
  14.24, 14.07, and 14.35 digits overall (worst species isoprene, 13.47 at the
  least): the 2026-10-16 calibration's numbers, to the hundredth. The prompt's
  "14 digits overall and at least 13.4 for every species" stands.

| Run | Model | Effort | Starts from | Harness | Session | Final commit | Status |
|---|---|---|---|---|---|---|---|
| astra-max | gpt-6-astra | max | `0dc7f67` | Codex CLI 0.163.0, interactive | 2026-10-22 09:58–11:06 | `cf16213` | done, verified, timed |
| astra-xhigh | gpt-6-astra | xhigh | `fbd365c` | Codex CLI 0.163.0, interactive | 2026-10-22 12:06–13:20 | `1c9112c` | done, verified, timed |
| astra-high | gpt-6-astra | high | `4055d7c` | Codex CLI 0.163.0, interactive | 2026-10-22 13:56–15:03 | `f478080` | done, verified, timed |
| sol-max | gpt-6.1-sol | max | `001fb64` | Codex CLI 0.163.0, interactive | 2026-10-22 15:39–17:21 | `ca0119e` | done, verified, timed |
| sol-xhigh | gpt-6.1-sol | xhigh | `171f54a` | Codex CLI 0.163.0, interactive | 2026-10-22 17:58–18:33 | `7af0738` | done, verified, timed |
| sol-high | gpt-6.1-sol | high | `1d944f4` | Codex CLI 0.163.0, interactive | 2026-10-22 19:20–20:34 | `07848a8` | done, verified, timed |
| luna-max | gpt-6-luna | max | `da2d722` | Codex CLI 0.163.0, interactive | 2026-10-22 21:11– | | running |
| luna-xhigh | gpt-6-luna | xhigh | `ada55d5` | Codex CLI 0.163.0, interactive | | | set up |
| luna-high | gpt-6-luna | high | `d3c2480` | Codex CLI 0.163.0, interactive | | | set up |

## Tier 2: astra-max

- Held-out data first: MICM's results at the gate's settings and the
  references for the five held-out scenarios were made from 09:40 to 09:57,
  outside the run root (references converged to 10.10–10.44 digits; MICM at the
  gate's settings 4.73–5.10 digits from them; MICM's default fails on one and
  keeps 2.36–3.64 digits on the others), so that nothing heavy ran during the
  run. The SSE's own Codex session in CERES-MIP had ended by the launch.
- Session start: 2026-10-22 09:58:25, launched by Claude Code with
  `admin/launch_codex_run.sh --tier2 astra-max gpt-6-astra max`, in a detached
  tmux session (socket `turbomicm`, session `tier2-astra-max`); the goal was
  sent at 09:58:30. The session record shows gpt-6-astra at max, base
  instructions from gpt-6-astra, Codex 0.163.0, approval `never` and sandbox
  `danger-full-access`, and no model switch; the process's environment holds no
  Claude Code variables besides the SSE's own `CLAUDE_CODE_DISABLE_ALTERNATE_SCREEN`.
- Codex session `01a163e7-f388-7153-baf9-d8a8478cce81`; process 15395. Codex's
  weekly usage window read 17 % before the launch (resets 2026-10-27 02:33);
  the SSE has a usage reset available.
- `caffeinate -ims -w 15395` keeps the machine awake while the session lives.
- Progress, from Claude Code's read-only checks every 20 minutes: at 20
  minutes, a first generator and specialized build passing the 100-cell MICM
  check but slower than the generic solver, and a first day-night scenario
  agreeing with MICM to 14.29 digits; at 40, its first commit (10:26) and a
  parameter sweep; at 60, full-size timings and both compilers' gates. The
  compliance scan (`admin/scan_commands.py`) found no command reaching outside
  the run directory and its three inputs, none naming another run or a
  private folder, and no install or download.
- Session end: the one turn ended at 11:06:18 with DONE (1 h 7 min 48 s from
  the goal), no context compaction, no warning; Codex showed "Goal achieved
  (1h 7m)". Every turn context shows gpt-6-astra at max, approval `never`,
  sandbox `danger-full-access`. The record was idle 10 minutes by 11:16; Claude
  Code quit Codex (`/quit`, 11:16:56) and closed the tmux session.
- Final commit `cf16213` (11:05), 2 commits after `0dc7f67` (`30c5ebf` at 10:26,
  "Generate mechanism-specific chemistry, sparse solves, and scenario runner";
  `cf16213`, "Tune generated kernels and document complete compiler and
  scenario validation"): 22 files changed, 3,911 lines added and 127 removed,
  the generic library (`src/`, `include/`) untouched. The working tree is clean.
  Collected as `tier2/astra/max`.
- Tokens: 182,784 total; input 106,808 (plus 7,022,592 cached); output 75,976
  (33,300 reasoning). Codex's weekly usage window went from 17 % to 18 %.
- Other load during the session: Claude Code's read-only progress checks and
  its writing of scripts and records. No builds, tests, or timing by the
  experimenter.
- Self-reported, not yet verified: both compilers pass the 680 cases, the
  generator's tests (TS1, MICM's three tutorial mechanisms, every rate family,
  a mechanism with no reactions), and every disclosed scenario at more than 14
  digits; on 10,000 cells at 18 threads, 3.063 s specialized against 3.872 s
  generic and 11.906 s for MICM (1.26 and 3.89 times), and on one thread 37.4,
  45.1, and 136.2 s.
- Transcript archived as `transcripts/tier2/astra-max.jsonl.gz` (585 lines,
  1.8 MB). The credential scan's one match (of the `sk-` pattern) was inside
  encrypted reasoning.
- Verified, 11:18–11:22, with `admin/verify_tier2.py`
  (`results/2026-10-22-tier2-astra-max-verification`): both compilers build
  without warnings (clean builds 17.8 s with clang, 52.6 s with gcc-15); 680
  PASS, 0 FAIL, 0 N/A with each, every id once; the Tier 1 TS1 gate passes
  (2.96e-13 with clang, 2.75e-13 with gcc-15); bin/ts1_bench_specialized passes
  on the seed-2026 cells (the same) and on the undisclosed input (1.4e-13,
  1.5e-13), with MICM's 281,700 internal steps; `make test-specialized` passes
  (TS1 and five other mechanisms); the day-night gate passes on all ten
  scenarios with both compilers, 14.05 to 14.45 digits overall and at least
  13.43 per species, with MICM's internal steps exactly; the audit finds
  nothing outside the rules in 44 C files, generated included. Generated C for
  TS1: 114,209 lines, 2.86 MB, from a 337-line generator.
- Timed, 11:31–11:38, provisionally at 18 threads (`admin/time_tier2.sh`, 5
  rounds, the undisclosed input; idle for 10 minutes first; all 30 runs pass
  the gate): specialized 3,084.5 ms; its generic bin/ts1_bench_cpu 3,797.3 ms
  and its Tier 1 commit's 3,803.8 ms (built here with the libomp path on the
  make command line); MICM 12,003.8 ms by default, 10,603.9 ms with 32-cell
  groups, 4,426.7 ms in its fastest configuration (vector8, 32-cell groups).
  Specialization gains 1.23 times over its own generic solver; 3.89 times
  MICM's default and 1.44 times its fastest. `time_tier2.sh` then failed to
  write its exit status (`status` is read-only in zsh; fixed), which delayed
  the next launch by about 25 minutes.

## Tier 2: astra-xhigh

- Session start: 2026-10-22 12:06:52, launched by Claude Code with
  `admin/launch_codex_run.sh --tier2 astra-xhigh gpt-6-astra xhigh`, once
  astra-max was collected, archived, verified, and timed; the goal was sent at
  12:06:57. The session record shows gpt-6-astra at xhigh, base instructions
  from gpt-6-astra, Codex 0.163.0, approval `never` and sandbox
  `danger-full-access`, and no model switch; the process's environment holds no
  Claude Code variables besides the SSE's own `CLAUDE_CODE_DISABLE_ALTERNATE_SCREEN`.
- Codex session `01a1645d-8d6b-7ea1-adde-c65151f4d896`; process 40578;
  `caffeinate -ims -w 40578`. No other Codex session was running; the machine
  was idle (load averages 1.35, 1.41, 2.12) and on AC power.
- Progress, from Claude Code's read-only checks every 20 minutes: at 20
  minutes, a first commit (12:19) extending its Tier 1 generator to the LU
  factorization, triangular solves, and the stages, passing the TS1 gate and
  diurnal-30min (14.37 digits) but slower than Tier 1 (3.16 against 3.00 s);
  at 40, four-cell tiles, reciprocal pivots, and Doolittle LU at 1.99 s, and a
  second commit (12:45); at 60, both compilers' full checks. The compliance
  scan found nothing outside the run directory and its inputs.
- Codex's usage window: 18 % at the launch, 21 % at 13:04; at 13:11:55 it
  read 0 %, with a new window ending 2026-10-29 13:11: the SSE's usage reset.
- Session end: the one turn ended at 13:20:57 with DONE (1 h 14 min 0 s from
  the goal), no context compaction, no warning; Codex showed "Goal achieved
  (1h 13m)". Every turn context shows gpt-6-astra at xhigh, approval `never`,
  sandbox `danger-full-access`. The record was idle 10 minutes by 13:31;
  Claude Code quit Codex (`/quit`, 13:31:34) and closed the tmux session.
- Final commit `1c9112c` (13:20), 3 commits after `fbd365c` (`c64e3c1`, 12:19,
  "Generate fixed-mechanism Rosenbrock kernels and add scenario driver";
  `e47de9b`, 12:45, "Tile generated solvers and accumulate sparse LU entries
  locally"; `1c9112c`, "Verify both compilers, scenario gates and benchmark
  results"): 33 files changed, 3,842 lines added and 182 removed. In `src/`
  and `include/`: the benchmark's input and output moved to `src/bench_io.h`,
  a new day-night driver `src/scenario.c`, a specialized-solver header, and one
  line of `src/state.c`. The working tree is clean. Collected as
  `tier2/astra/xhigh`.
- Tokens: 369,502 total; input 266,838 (plus 13,663,872 cached); output
  102,664 (57,149 reasoning).
- Other load during the session: Claude Code's read-only progress checks and
  record keeping. No builds, tests, or timing by the experimenter.
- Self-reported, not yet verified: both compilers pass the 680 cases, 70
  specialized cases, and every disclosed scenario at more than 14 digits;
  2.012 s at 18 threads, 1.55 times its Tier 1 benchmark and 5.90 times MICM.
- Transcript archived as `transcripts/tier2/astra-xhigh.jsonl.gz` (925 lines,
  2.7 MB). The credential scan's three matches were inside encrypted content.
- Verified, 13:32–13:37, with `admin/verify_tier2.py`
  (`results/2026-10-22-tier2-astra-xhigh-verification`): both compilers build
  without warnings (clean builds 71.2 s with clang, 73.6 s with gcc-15); 680
  PASS, 0 FAIL, 0 N/A with each; the Tier 1 TS1 gate passes (2.03e-13,
  2.14e-13); bin/ts1_bench_specialized passes on the seed-2026 cells (2.96e-13,
  2.75e-13) and the undisclosed input (1.4e-13, 1.5e-13), with MICM's 281,700
  internal steps; `make test-specialized` passes (70 cases, seven mechanisms);
  the day-night gate passes on all ten scenarios with both compilers, 14.05 to
  14.45 digits overall and at least 13.31 per species, with MICM's internal
  steps exactly; the audit finds nothing outside the rules in 82 C files.
  Generated C for TS1: `generated/ts1_specialized.c`, 75,172 lines, 2.61 MB;
  the Tier 2 generator adds 492 lines of Python and a 112-line C template to
  Tier 1's generator.
- Timed, 13:46–13:54, provisionally at 18 threads (5 rounds, the undisclosed
  input, after 10 idle minutes; all 35 runs pass the gate): specialized
  2,098.4 ms; its generic 3,170.3 ms and its Tier 1 commit's 3,152.0 ms; MICM
  11,969.9 ms by default, 10,630.6 ms with 32-cell groups, 4,423.0 ms in its
  fastest configuration; astra-max's specialized 3,045.5 ms (its own session:
  3,084.5). Specialization gains 1.51 times over its own generic solver; 5.70
  times MICM's default and 2.11 times its fastest; 1.45 times astra-max's
  specialized solver.

## Tier 2: astra-high

- Session start: 2026-10-22 13:56:22, launched by Claude Code with
  `admin/launch_codex_run.sh --tier2 astra-high gpt-6-astra high`, once
  astra-xhigh was collected, archived, verified, and timed; the goal was sent
  at 13:56:25. The session record shows gpt-6-astra at high, base instructions
  from gpt-6-astra, Codex 0.163.0, approval `never` and sandbox
  `danger-full-access`, and no model switch; the process's environment holds no
  Claude Code variables besides the SSE's own `CLAUDE_CODE_DISABLE_ALTERNATE_SCREEN`.
- Codex session `01a164c1-c86d-74f1-a385-74d24b64d3a1`; process 58197;
  `caffeinate -ims -w 58197`. No other Codex session was running; AC power.
- To watch: astra-high's Tier 1 benchmark recomputes every rate constant before
  each timed step (`src/ts1_bench.c:154`), which MICM's driver and the other
  entries do once, before the warm-up.
- Progress, from Claude Code's read-only checks every 20 minutes: at 20
  minutes, a first commit (14:05) of generated chemistry, Jacobian, LU, and
  solve kernels, matching MICM's 100-cell result exactly and diurnal-10min to
  14.29 digits; at 40, a second commit (14:29) adding fixed 32-cell strides,
  a fused error-estimate loop, and an exact rate cache, which skips the rate
  constants only when every condition and custom parameter is byte for byte
  unchanged (a generator option, `GEN_OPTIONS=--no-cache`); at 66, its
  documents. The compliance scan found nothing outside the run directory and
  its inputs (its one apparent hit, markdown text naming `../micm` in
  backticks, led to a fix of the scanner).
- Session end: the one turn ended at 15:03:42 with DONE (1 h 7 min 16 s from
  the goal), no context compaction, no warning; Codex showed "Goal achieved
  (1h 7m)". Every turn context shows gpt-6-astra at high, approval `never`,
  sandbox `danger-full-access`. Claude Code quit Codex (`/quit`, 15:14:32) and
  closed the tmux session.
- Final commit `f478080` (15:03), 3 commits after `4055d7c` (`e78e3a1`, 14:05,
  "Generate mechanism-specific chemistry and sparse solver kernels with
  scenario driver"; `617b634`, 14:29, "Tune generated kernels and add exact
  rate caching and regeneration validation"; `f478080`, "Document
  specialization experiments and verify both compiler gates"): 20 files
  changed, 1,656 lines added and 203 removed, 454 and 193 of them in `src/`
  and `include/` (the Rosenbrock step moved into a template shared with the
  generated solver, the benchmark driver into `src/ts1_driver.h`, the cache's
  fields in the state, and a new day-night driver). The working tree is
  clean. Collected as `tier2/astra/high`.
- Tokens: 329,804 total; input 242,781 (plus 9,947,264 cached); output 87,023
  (45,333 reasoning). Codex's usage window went from 1 % to 2 % (of the
  window that began with the SSE's reset).
- Other load during the session: Claude Code's read-only progress checks and
  record keeping. No builds, tests, or timing by the experimenter.
- Self-reported, not yet verified: both compilers pass the 680 cases, the
  specialized tests, and every disclosed scenario; 1.25 times its Tier 1
  benchmark and 3.53 times MICM at 18 threads.
- Transcript archived as `transcripts/tier2/astra-high.jsonl.gz` (768 lines,
  2.3 MB). The credential scan's one match was inside encrypted content.
- Verified, 15:14–15:18, with `admin/verify_tier2.py`
  (`results/2026-10-22-tier2-astra-high-verification`): both compilers build
  without warnings (clean builds 25.0 s with clang, 37.0 s with gcc-15); 680
  PASS, 0 FAIL, 0 N/A with each; the Tier 1 TS1 gate passes (0 with clang,
  1.43e-13 with gcc-15); bin/ts1_bench_specialized passes on the seed-2026
  cells (2.31e-14, 1.43e-13) and the undisclosed input (1.95e-14, 1.52e-13),
  with MICM's 281,700 internal steps; `make test-specialized` passes (TS1 and
  five other mechanisms, one a mutated TS1); the day-night gate passes on all
  ten scenarios with both compilers, 14.08 to 14.46 digits overall and at
  least 13.42 per species, with MICM's internal steps exactly; the audit finds
  nothing outside the rules in 75 C files. Generated C for TS1:
  `build/specialized/ts1.c`, 44,974 lines, 2.05 MB, from a 406-line generator.
- Timed, 15:27–15:37, provisionally at 18 threads (5 rounds, the undisclosed
  input, after 10 idle minutes; all 45 runs pass the gate): specialized
  3,332.7 ms, and 3,407.1 ms without its rate cache (a variant built here with
  `GEN_OPTIONS=--no-cache`); its generic 4,172.8 ms and its Tier 1 commit's
  4,183.9 ms; MICM 11,807.4, 10,516.0, and 4,402.3 ms; astra-xhigh's
  specialized 2,053.4 ms and astra-max's 2,975.6 ms. Specialization gains 1.25
  times over its own generic solver (1.22 without the cache); 3.54 times MICM's
  default and 1.32 times its fastest. Every program ran 1.4 to 3.5 % faster
  than in the two earlier sessions.

## Tier 2: sol-max

- Session start: 2026-10-22 15:39:19, launched by Claude Code with
  `admin/launch_codex_run.sh --tier2 sol-max gpt-6.1-sol max`, once astra-high
  was collected, archived, verified, and timed; the goal was sent at 15:39:24.
  The session record shows gpt-6.1-sol at max, base instructions from
  gpt-6.1-sol, Codex 0.163.0, approval `never` and sandbox `danger-full-access`,
  and no model switch; the process's environment holds no Claude Code variables
  besides the SSE's own `CLAUDE_CODE_DISABLE_ALTERNATE_SCREEN`.
- Codex session `01a16520-0dd1-7032-bd32-99fef62218eb`; process 74257;
  `caffeinate -ims -w 74257`. No other Codex session was running; AC power.
- sol-max's Tier 1 entry already generates TS1-specific kernels
  (`chem_load_ts1(1)`, `src/bench.c:103`), as astra-xhigh's does.
- Progress, from Claude Code's read-only checks every 20 minutes: at 20
  minutes, a first commit (15:50) of complete generated solvers and a scenario
  driver, diurnal-10min at 14.29 digits, and tests on seven mechanisms (TS1
  with an added reaction and a changed coefficient among them); at 40, four-cell
  tiles inside each 32-cell group (16:18); at 66, a fix for regeneration when
  an included configuration file changes, and 42 generator checks; at 96,
  both compilers' final passes. The compliance scan found nothing outside the
  run directory and its inputs.
- Session end: the one turn ended at 17:21:30 with DONE (1 h 42 min 6 s from
  the goal), no context compaction, no warning. Every turn context shows
  gpt-6.1-sol at max, approval `never`, sandbox `danger-full-access`. The record
  was idle 10 minutes by 17:31; Claude Code quit Codex (`/quit`, 17:32:25) and
  closed the tmux session, with the new `admin/collect_tier2.sh`, which also
  collects the branch and archives the transcript.
- Final commit `ca0119e` (17:20), 4 commits after `001fb64` (`2b71739`, 15:50,
  "Generate complete mechanism-specialized Rosenbrock solvers and scenario
  driver"; `7ef74b6`, 16:18, "Tile independent cell stages while preserving
  grouped Rosenbrock error control"; `38792b6`, 16:50, "Track recursive
  configuration inputs and strengthen specialization validation"; `ca0119e`,
  "Record compiler gates, full-size timings, build costs, and day-night
  scores"): 22 files changed, 1,858 lines added and 322 removed; in `src/`, the
  benchmark's input and output moved to `src/bench_io.h` and a new day-night
  driver, `src/scenario.c`. The working tree is clean. Collected as
  `tier2/sol/max`.
- Tokens: 549,677 total; input 350,337 (plus 27,473,152 cached); output 199,340
  (133,422 reasoning). Codex's usage window went from 2 % to 3 %.
- Other load during the session: Claude Code's read-only progress checks and
  record keeping. No builds, tests, or timing by the experimenter.
- Self-reported, not yet verified: both compilers pass 675 cases (5 N/A, as in
  Tier 1), 42 generator checks, and every disclosed scenario; 2,330 ms at 18
  threads against 2,928 ms for its Tier 1 benchmark and 9,990 ms for MICM,
  1.26 and 4.29 times.
- Transcript archived as `transcripts/tier2/sol-max.jsonl.gz` (1,582 lines,
  4.4 MB). The credential scan's four matches were inside encrypted content.
- Verified, 17:32–17:37, with `admin/verify_tier2.py`
  (`results/2026-10-22-tier2-sol-max-verification`): both compilers build
  without warnings (clean builds 31.3 s with clang, 56.1 s with gcc-15); 675
  PASS, 0 FAIL, 5 N/A with each (the five C++-only cases N/A in its Tier 1
  entry too), every id once; the Tier 1 TS1 gate passes (1.26e-13, 1.16e-13);
  bin/ts1_bench_specialized passes on the seed-2026 cells (the same) and the
  undisclosed input (9.98e-14, 1.48e-13), with MICM's 281,700 internal steps;
  `make test-specialized` passes (42 checks); the day-night gate passes on all
  ten scenarios with both compilers, 14.06 to 14.47 digits overall and at least
  13.44 per species, with MICM's internal steps exactly; the audit finds
  nothing outside the rules in 37 C files. Generated C for TS1:
  `generated/ts1_specialized.c`, 66,075 lines, 2.22 MB.
- Timed, 17:47–17:56, provisionally at 18 threads (5 rounds, the undisclosed
  input, after 10 idle minutes; all 45 runs pass the gate): specialized
  2,511.8 ms; its generic 3,065.5 ms and its Tier 1 commit's 3,116.1 ms; MICM
  12,061.1, 10,625.4, and 4,475.8 ms; astra-xhigh's specialized 2,085.1 ms,
  astra-max's 3,005.4 ms, astra-high's 3,256.1 ms. Specialization gains 1.22
  times over its own generic solver; 4.80 times MICM's default and 1.78 times
  its fastest.

## Tier 2: sol-xhigh

- Session start: 2026-10-22 17:58:32, launched by Claude Code with
  `admin/launch_codex_run.sh --tier2 sol-xhigh gpt-6.1-sol xhigh`, once sol-max
  was collected, archived, verified, and timed; the goal was sent at 17:58:36.
  The session record shows gpt-6.1-sol at xhigh, base instructions from
  gpt-6.1-sol, Codex 0.163.0, approval `never` and sandbox `danger-full-access`,
  and no model switch; the process's environment holds no Claude Code variables
  besides the SSE's own `CLAUDE_CODE_DISABLE_ALTERNATE_SCREEN`.
- Codex session `01a1659f-81c3-7d00-b307-1f319db2589c`; process 90297;
  `caffeinate -ims -w 90297`. No other Codex session was running; AC power.
- Progress, from Claude Code's read-only checks every 20 minutes: at 20
  minutes, a first commit (18:09) of generated fixed-index solvers and a
  day-night driver, its 100-cell result bit-identical to MICM's, diurnal-10min
  at 14.29 digits, generated-against-generic tests on five mechanisms; tiling
  tried and left off (slower here). The turn ended before the 40-minute check.
  The compliance scan found nothing outside the run directory and its inputs
  (its apparent hits, relative paths after a `cd` inside the command, led to
  the scanner's following such `cd`s).
- Session end: the one turn ended at 18:33:12 with DONE (34 min 35 s from the
  goal, the shortest Tier 2 run so far), no context compaction, no warning.
  Every turn context shows gpt-6.1-sol at xhigh, approval `never`, sandbox
  `danger-full-access`. The record was idle 10 minutes by 18:43; Claude Code
  quit Codex (18:43:36) and collected and archived it with
  `admin/collect_tier2.sh`.
- Final commit `7af0738` (18:32), 3 commits after `171f54a` (`e1b1338`, 18:09,
  "Generate fixed-index mechanism solvers and add day-night driver";
  `148be07`, 18:24, "Bound generated kernels, select measured defaults and
  verify optional transformations"; `7af0738`, "Record cross-compiler gates,
  benchmark results and completed specialization audit"): 20 files changed,
  2,103 lines added and 227 removed; in `src/`, the benchmark's common code
  moved to `src/bench_common.h`, and new `src/scenario.c`,
  `src/specialized.c`, and `src/specialized_solver.h`. The working tree is
  clean. Collected as `tier2/sol/xhigh`.
- Tokens: 207,782 total; input 122,093 (plus 10,847,872 cached); output 85,689
  (36,129 reasoning). Codex's usage window stayed at 3 %.
- Other load during the session: Claude Code's read-only progress checks,
  record keeping, and copying two mechanisms from MUSICA for the review's
  exercises. No builds, tests, or timing by the experimenter.
- Self-reported, not yet verified: both compilers pass the Tier 1 suite, 70
  specialized checks, and every disclosed scenario; 2.61 s at 18 threads and
  36.17 s on one, 1.32 and 1.21 times its generic solver.
- Transcript archived as `transcripts/tier2/sol-xhigh.jsonl.gz` (728 lines,
  2.2 MB). The credential scan's two matches were inside encrypted content.
- Verified, 18:43–18:47, with `admin/verify_tier2.py`
  (`results/2026-10-22-tier2-sol-xhigh-verification`): both compilers build
  without warnings (clean builds 18.3 s with clang, 46.1 s with gcc-15); 675
  PASS, 0 FAIL, 5 N/A with each (as in Tier 1); the Tier 1 TS1 gate passes (0
  with clang, 1.43e-13 with gcc-15); bin/ts1_bench_specialized passes on the
  seed-2026 cells (2.31e-14, 1.43e-13) and the undisclosed input (1.95e-14,
  1.52e-13), with MICM's 281,700 internal steps; `make test-specialized` passes
  (70 checks); the day-night gate passes on all ten scenarios with both
  compilers, 14.08 to 14.46 digits overall and at least 13.42 per species, with
  MICM's internal steps exactly; the audit finds nothing outside the rules in
  65 C files. Generated C for TS1: `generated/sp_ts1.c`, 112,387 lines, 3.43 MB.
- Its outputs are byte-identical to astra-high's in all six TS1 runs and 18 of
  the 20 day-night runs (the other two, clean-10min with each compiler, agree
  to 15.4 digits). Checked as a possible copy, it is convergence: the two
  generators are independent work (374 and 406 lines sharing 7 substantive
  lines, all boilerplate), no command of sol-xhigh's session lists the run root
  or reaches into another run, and both entries reproduce MICM's operation
  order and emit literal-coefficient rate formulas in it.
- Timed, 19:08–19:18, provisionally at 18 threads (5 rounds, the undisclosed
  input; the idle wait restarted after more than one core of other load at
  18:57; all 50 runs pass the gate): specialized 2,741.5 ms; its generic
  3,672.0 ms and its Tier 1 commit's 3,664.9 ms; MICM 11,977.9, 10,993.3, and
  4,452.4 ms; sol-max's specialized 2,474.6, astra-xhigh's 2,068.9,
  astra-max's 3,007.4, astra-high's 3,233.4 ms. Specialization gains 1.34 times
  over its own generic solver; 4.37 times MICM's default and 1.62 times its
  fastest.

## Tier 2: sol-high

- Session start: 2026-10-22 19:20:06, launched by Claude Code with
  `admin/launch_codex_run.sh --tier2 sol-high gpt-6.1-sol high`, once sol-xhigh
  was collected, archived, verified, and timed; the goal was sent at 19:20:11.
  The session record shows gpt-6.1-sol at high, base instructions from
  gpt-6.1-sol, Codex 0.163.0, approval `never` and sandbox `danger-full-access`,
  and no model switch; the process's environment holds no Claude Code variables
  besides the SSE's own `CLAUDE_CODE_DISABLE_ALTERNATE_SCREEN`.
- Codex session `01a165ea-3043-75a0-876c-f74b06d91fbf`; process 2771;
  `caffeinate -ims -w 2771`. No other Codex session was running; AC power.
- Progress, from Claude Code's read-only checks every 20 minutes: at 20
  minutes, a first commit (19:36), the TS1 gate (4.25e-14) and diurnal-10min
  (14.29 digits, MICM's step count) passing, and generated-against-generic
  tests on seven configurations; at 40, tile and layout experiments; at 60, a
  second commit (20:10) passing both compilers' gates. Its final review found
  a generator bug, negative fixed coefficients emitted as an invalid C token,
  which it fixed and tested before its last commit. The compliance scan found
  nothing outside the run directory and its inputs (its one lead, a path
  inside a Python string that writes its own test scenario, reviewed).
- Session end: the one turn ended at 20:34:51 with DONE (1 h 14 min 40 s from
  the goal), no context compaction, no warning. Every turn context shows
  gpt-6.1-sol at high, approval `never`, sandbox `danger-full-access`. Claude
  Code quit Codex (20:45:16) and collected and archived it with
  `admin/collect_tier2.sh`.
- Final commit `07848a8` (20:34), 3 commits after `1d944f4` (`2587d2e`, 19:36,
  "Add mechanism generator and specialized benchmark and scenario drivers";
  `0590848`, 20:10, "Optimize generated sparse plans and stage arithmetic; pass
  both compiler gates"; `07848a8`, "Cover negative rate coefficients and record
  complete compiler gates and timings"): 23 files changed, 3,536 lines added
  and 328 removed; in `src/`, the Rosenbrock algorithm moved into a template
  shared by the generic and generated solvers (`src/rosenbrock_impl.h`), the
  benchmark's input and output into `src/bench_io.h`, and a new day-night
  driver. The working tree is clean. Collected as `tier2/sol/high`.
- Tokens: 451,395 total; input 294,639 (plus 18,595,456 cached); output 156,756
  (103,978 reasoning). Codex's usage window stayed at 3 %.
- Other load during the session: Claude Code's read-only progress checks and
  writing of the blind review's instructions. No builds, tests, or timing by
  the experimenter.
- Self-reported, not yet verified: both compilers pass 675 cases, 36 generator
  comparisons, and every disclosed scenario; 3.01 s at 18 threads against 4.33 s
  for its Tier 1 benchmark (1.44 times), and 35.05 against 46.48 s on one.
- Transcript archived as `transcripts/tier2/sol-high.jsonl.gz` (1,140 lines,
  3.6 MB). The credential scan's four matches were inside encrypted content.
- Verified, 20:45–20:48, with `admin/verify_tier2.py`
  (`results/2026-10-22-tier2-sol-high-verification`): both compilers build
  without warnings (clean builds 17.2 s with clang, 23.8 s with gcc-15); 675
  PASS, 0 FAIL, 5 N/A with each (as in Tier 1); the Tier 1 TS1 gate passes
  (2.03e-13, 2.14e-13); bin/ts1_bench_specialized passes on the seed-2026 cells
  (2.96e-13, 2.75e-13) and the undisclosed input (1.4e-13, 1.5e-13), with
  MICM's 281,700 internal steps; `make test-specialized` passes (36 checks);
  the day-night gate passes on all ten scenarios with both compilers, 14.05 to
  14.45 digits overall and at least 13.31 per species, with MICM's internal
  steps exactly; the audit finds nothing outside the rules in 43 C files.
  Generated C for TS1: `build/ts1_specialized.c`, 26,300 lines, 0.99 MB, the
  least so far (its LU is a table run by shared kernels).
- Timed, 20:58–21:10, provisionally at 18 threads (5 rounds, the undisclosed
  input, after 10 idle minutes; all 55 runs pass the gate): specialized
  2,832.7 ms; its generic 3,865.4 ms and its Tier 1 commit's 3,891.8 ms; MICM
  11,863.6, 10,689.0, and 4,339.1 ms; sol-xhigh's specialized 2,792.9, sol-max's
  2,426.7, astra-xhigh's 2,042.8, astra-max's 3,001.8, astra-high's 3,193.4 ms.
  Specialization gains 1.37 times over its own generic solver; 4.19 times
  MICM's default and 1.53 times its fastest.

## Tier 2: luna-max

- Session start: 2026-10-22 21:11:19, launched by Claude Code with
  `admin/launch_codex_run.sh --tier2 luna-max gpt-6-luna max`, once sol-high
  was collected, archived, verified, and timed; the goal was sent at 21:11:24.
  The session record shows gpt-6-luna at max, base instructions from
  gpt-6-luna, Codex 0.163.0, approval `never` and sandbox `danger-full-access`,
  and no model switch; the process's environment holds no Claude Code variables
  besides the SSE's own `CLAUDE_CODE_DISABLE_ALTERNATE_SCREEN`.
- Codex session `01a16650-05af-7cf2-b421-d5800f93e165`; process 22711;
  `caffeinate -ims -w 22711`. No other Codex session was running; AC power.
- luna-max's Tier 1 entry keeps MICM's default grouping (one group per thread)
  and was slower than MICM's default (0.90 times at ten threads on the M5).

# Tier 3

Launch record for the Tier 3 runs (algorithm design; see `tier3/README.md`).

Each run works in `benchmarks/tier3/<run>/`, a clone of its own Tier 1 repository
at the Tier 1 run's final commit, made by `tier3/admin/setup_runs.py`: it shares
no files with the Tier 1 repository, has no remote, and its PROMPT.md holds the
Tier 3 prompt. Claude Code launches every run with

    admin/launch_codex_run.sh --tier3 <run> <model> <effort>

which types `codex --search -m <model> -c model_reasoning_effort=<effort> -c features.memories=false`
into a detached tmux session (socket `turbomicm`, session `tier3-<run>`), and
then the goal:

    /goal read PROMPT.md and work until done, commit

Order
- 2026-10-18: the runs go one at a time, in the order of their Tier 1 runs:
  sol-xhigh, sol-max, astra-xhigh, astra-max (the SSE's decision). Luna has no
  Tier 3 run: its arm waits on its Tier 1 runs, which have not been rerun since
  luna-high was stopped on 2026-10-16.

Undisclosed scenarios
- 2026-10-18, 01:50: the SSE approved five held-out scenarios, and their
  parameters were fixed, before any Tier 3 result, in a script kept outside the
  run root and this repository. Its SHA-256 is
  `08016cac9dfb41c1824648e843b2962686de94a1fcc2ddefb58a9445268e195c`. They vary
  what the agents' scenarios README announces: seeds, latitudes, seasons, and
  48-hour runs, with steps of 5 to 30 minutes and initial NOx from x0.1 to x10, as
  in the disclosed set. The script will be published here once the runs end; its
  scenario files and reference solutions are made between runs, never during one.
- 2026-10-18, 03:45–04:12, between sol-xhigh and sol-max: the held-out files
  (1,000 cells, and 10,000 for timing) and their references were made outside the
  run root, as for the disclosed set (`tier3/admin/make_references.sh`, MICM's
  scenario driver built afresh, 10 threads). The references converge to 10.16 to
  10.69 digits overall, with no species below 8.98. MICM's default three-stage
  Rosenbrock fails on two of the five (a NaN in step 1; step size too small) and
  reaches 1.43 to 3.14 digits on the other three. The 10,000-cell versions of the
  disclosed scenarios were made at the same time. Results on the held-out
  scenarios stay out of this repository and the wiki, by name, until the last
  run ends.

Launch changes
- 2026-10-18: `admin/launch_codex_run.sh` gained `--tier3`. In that mode it checks
  that the run directory is an untouched clone at the Tier 1 commit its prompt
  names (HEAD at that commit, no remote, nothing changed or added but PROMPT.md),
  adds `--search`, and names the tmux session `tier3-<run>`. Its check of the
  status line now looks for the model's name: the old check, for "GPT-6-", would
  have stopped every Sol launch, since Sol is GPT-6.1-Sol.
- 2026-10-18, 01:43: the Tier 3 mode was tested on a throwaway repository in
  Claude Code's scratch space, with gpt-6.1-sol at low effort and a one-line task
  that needed web search (session `01a14df7-2727-70b3-a00e-4769b27dc0a0`). The
  session searched the web (`web.run`, recorded as `web.search` items), opened two
  pages, committed, and ended with DONE in 21 s.

| Run | Model | Effort | Starts from | Harness | Session | Final commit | Status |
|---|---|---|---|---|---|---|---|
| sol-xhigh | gpt-6.1-sol | xhigh | `171f54a` | Codex CLI 0.163.0, interactive, `--search` | 2026-10-18 01:45–03:34 | `63d246f` | done, verified |
| sol-max | gpt-6.1-sol | max | `001fb64` | Codex CLI 0.163.0, interactive, `--search` | 2026-10-18 04:23–06:03 | `c122d12` | done, verified |
| astra-xhigh | gpt-6-astra | xhigh | `fbd365c` | Codex CLI 0.163.0, interactive, `--search` | 2026-10-18 06:20–08:18 | `db44838` | done, verified |
| astra-max | gpt-6-astra | max | `0dc7f67` | Codex CLI 0.163.0, interactive, `--search` | 2026-10-18 08:39–11:14 | `3cab602` | done, verified |

## Tier 3: sol-xhigh

- Session start: 2026-10-18 01:45:12, launched by Claude Code with
  `admin/launch_codex_run.sh --tier3 sol-xhigh gpt-6.1-sol xhigh`; the goal was
  sent at 01:45:16. The session record shows gpt-6.1-sol at xhigh, base
  instructions from gpt-6.1-sol, Codex 0.163.0, and no model switch; the `codex`
  command line carries `--search`.
- Codex session `01a14df8-cf84-7851-b9c3-3c76c4eef53e`.
- The machine had restarted at 01:31. At launch, macOS's photo analysis
  (`mediaanalysisd`) was still using part of one core.
- The whole run was one turn, from the goal at 01:45:16 to its end at 03:34:07
  (1 h 48 min 51 s), with no context compaction. Its one turn context shows
  gpt-6.1-sol at xhigh. The session marked its goal complete and ended with
  DONE. Claude Code closed it with `/quit` at 03:45, after the record had been
  idle for 10 minutes.
- Tokens: 661,749 total; input 421,143 (plus 24,919,424 cached); output 240,606
  (170,859 reasoning). Codex's weekly limit went from 19 % used to 23 %.
- Web: 13 web items: 3 searches (9 queries), 4 pages opened, 2 finds in a page,
  and 4 page reads. It searched for Sandu et al. (1997), KPP's Rosenbrock
  methods, Söderlind and Gustafsson's step-size control, SEULEX, and CVODE's
  Jacobian reuse. It read Sandu et al. (1997, part II), Söderlind (2002), KPP's
  documentation of its Rosenbrock methods (the source of its RODAS3.1
  coefficients), Dreger et al. (2025, GMD; the H211b controller), and Hairer and
  Lubich's overview. REFERENCES.md cites every page it opened.
- 14 commits after `171f54a`, from 01:49 to 03:33, ending `63d246f`, "Complete
  integration research, work-precision reports and final validation audit":
  68 files changed, 4,263 lines added and 281 removed. No remote; the working
  tree is clean.
- Compliance: no command reads outside its directory and the read-only inputs,
  apart from Codex's own check for AGENTS.md in the parent folders. It built its
  own copies of MICM's drivers in its `bin/` and made its own 48-hour test
  conditions with references. It installed nothing outside its conda
  environment.
- Self-reported, not yet verified: both compilers pass the Tier 1 gate; five
  operating points, four RODAS3.1 (rtol 2e-3 to 3e-4) and one six-stage RODAS4;
  the recommended RODAS3.1 point (rtol 3e-4, atol 1e-18) reaches 3.70 digits
  overall and 2.73 for every species on the disclosed scenarios, 3.3 to 8.6
  times faster than MICM's three-stage Rosenbrock tuned to 3 digits, on 1,000
  cells.
- Other load during the session: Claude Code checked progress read-only, wrote
  scripts and the run log, and at 01:54 ran one dry check of the verification
  script's audit stage (`make -n` and a grep, a few seconds at low priority). No
  builds, tests, or timing ran.
- Collected as `tier3/sol/xhigh`; transcript `transcripts/tier3/sol-xhigh.jsonl.gz`
  on the `transcripts` branch.
- Verified, 04:17–04:22, with `admin/verify_tier3.py`
  (`results/2026-10-18-tier3-sol-xhigh-verification`): both compilers build
  without warnings; 675 PASS, 0 FAIL, 5 N/A with each; the Tier 1 TS1 gate
  passes (largest relative difference 0 with clang, 1.43e-13 with gcc-15); the
  audit finds nothing outside the rules. All five points pass the gate on the 6
  disclosed and 5 held-out scenarios. L2 and L2_margin reach level 2 on all 11,
  and L3, recommended, and L4 reach level 3 on all 11. The disclosed scores match
  the agent's RESULTS.md to the third decimal. A first pass at 04:12 found two
  bugs in the verification script, since fixed: it counted the case list's 4
  comment lines as missing ids, and its audit scanned the conda package cache.
- The session record's world state lists the SSE's approved Codex command
  prefixes, which name other projects; the Tier 1 transcripts carry the same
  list. It holds no credential values.

## Tier 3: sol-max

- Session start: 2026-10-18 04:23:36, launched by Claude Code with
  `admin/launch_codex_run.sh --tier3 sol-max gpt-6.1-sol max`; the goal was sent
  at 04:23:40. The session record shows gpt-6.1-sol at max, base instructions
  from gpt-6.1-sol, Codex 0.163.0, and no model switch; the `codex` command line
  carries `--search`.
- Codex session `01a14e89-d573-77d2-8286-e1d020bdb20d`.
- Launched after sol-xhigh was collected, archived, and verified, and after the
  held-out scenarios were made. The machine was idle, on AC power, with load
  averages falling from that work (2.05 over 1 minute).
- The whole run was one turn, from the goal at 04:23:40 to its end at 06:03:04
  (1 h 39 min 23 s), with no context compaction. Its one turn context shows
  gpt-6.1-sol at max. The session marked its goal complete and ended with DONE.
  Claude Code closed it with `/quit` at 06:14, after the record had been idle
  for 10 minutes.
- Tokens: 775,162 total; input 598,037 (plus 26,008,576 cached); output 177,125
  (106,937 reasoning). Codex's weekly limit went from 24 % used to 27 %.
- Web: 10 web items: 3 searches (6 queries), 2 pages opened, 3 finds in a page,
  and 2 page reads. It searched for Sandu et al. (1997), Rosenbrock W-methods
  (ROS3P, ROS3PRL2), and accuracy targets for chemistry solvers, and read Verwer
  et al. (1999) on ROS2 and Dreger et al. (2025, GMD). REFERENCES.md cites both
  pages it opened.
- 4 commits after `001fb64`, from 04:33 to 06:02, ending `c122d12`, "Finalize
  compiler gates, operating-point benchmarks, and research reports": 88 files
  changed, 26,469 lines added and 241 removed. No remote; the working tree is
  clean.
- Compliance: no command reads outside its directory and the read-only inputs,
  apart from Codex's own check for AGENTS.md in the parent folders; no install or
  download command.
- Self-reported, not yet verified: both compilers pass the Tier 1 gate (675 PASS,
  5 N/A); five RODAS3 points at rtol 1e-2 to 1e-4 (atol 1e-18), all with 16-cell
  groups and dynamic scheduling; the recommended point (rtol 3e-4) reaches 3.57
  digits overall and 2.55 for every species on the disclosed scenarios, 3.5 to
  11.6 times faster than MICM's three-stage Rosenbrock at rtol 3e-4 and atol
  1e-18, on 1,000 cells.
- Other load during the session: Claude Code made two read-only progress checks
  and edited the run log and wiki. No builds, tests, or timing ran.
- Collected as `tier3/sol/max`; transcript `transcripts/tier3/sol-max.jsonl.gz`.
- Verified, 06:15–06:19, with `admin/verify_tier3.py`
  (`results/2026-10-18-tier3-sol-max-verification`): both compilers build
  without warnings; 675 PASS, 0 FAIL, 5 N/A with each; the Tier 1 TS1 gate
  passes (largest relative difference 1.26e-13 with clang, 1.16e-13 with
  gcc-15); the audit finds nothing outside the rules. All five points pass the
  gate on the 6 disclosed and 5 held-out scenarios. fast and balanced reach
  level 2 on all 11, and accurate, recommended, and precise reach level 3 on all
  11. The disclosed scores match the agent's OPERATING_POINTS.md to the third
  decimal.

## Tier 3: astra-xhigh

- Session start: 2026-10-18 06:20:30, launched by Claude Code with
  `admin/launch_codex_run.sh --tier3 astra-xhigh gpt-6-astra xhigh`; the goal was
  sent at 06:20:34. The session record shows gpt-6-astra at xhigh, base
  instructions from gpt-6-astra, Codex 0.163.0, and no model switch; the `codex`
  command line carries `--search`.
- Codex session `01a14ef4-db46-7743-8cf5-57c66c1aa74c`.
- Launched after sol-max was collected, archived, and verified. The machine was
  on AC power, with load averages falling from the verification (3.44 over 1
  minute).
- The whole run was one turn, from the goal at 06:20:34 to its end at 08:18:28
  (1 h 57 min 54 s), with one context compaction, at 07:29. Both turn contexts
  show gpt-6-astra at xhigh. The session marked its goal complete and ended with
  DONE. Claude Code closed it with `/quit` at 08:29, after the record had been
  idle for 10 minutes. Codex's one warning was its startup notice that `-c` and
  `--search` run it without the shared background server, as for every launch.
- Tokens: 684,296 total; input 467,979 (plus 18,585,984 cached); output 216,317
  (133,384 reasoning). Codex's weekly limit went from 28 % used to 33 %.
- Web: 7 web items: 3 searches (6 queries), 3 page openings (6 pages, one of
  them refused by a CAPTCHA), and 1 page read. It searched for Sandu et al.
  (1997), Rosenbrock W-methods with Jacobian reuse, Verwer and Blom's chemistry
  integrators, Alexander's (1977) diagonally implicit Runge–Kutta methods, and
  linearly implicit Euler extrapolation (Deuflhard; Hairer). It read Sandu et
  al. (1997, part II), Verwer et al. (1999) on ROS2, Butcher's tutorial on
  implicit Runge–Kutta methods, Hairer and Lubich's overview, and the abstract
  of Hairer and Lubich (1988). REFERENCES.md cites every page it read, and part
  I of Sandu et al. from a search result.
- 6 commits after `fbd365c`, from 06:26 to 08:17, ending `db44838`, "Complete
  operating-point evaluation, work-precision figures and final research
  reports": 106 files changed, 20,069 lines added and 277 removed. No remote;
  the working tree is clean.
- Compliance: no command reads outside its directory and the read-only inputs,
  apart from Codex's own check for AGENTS.md in the parent folders; no install or
  download command, and no path into another run or the held-out folder.
- Timer: `src/scenario.c` times each block of up to four steps (cut at
  snapshots) whole: the photolysis update at the start of every step, the rate
  constants (all of them in the first step, then only those that depend on
  photolysis, since all other conditions stay constant), the integration, and a
  per-step check for non-finite or negative values. Setup, reading the input,
  and writing snapshots and output are outside it.
- Self-reported, not yet verified: both compilers pass the Tier 1 gate; five
  points, all with atol 1e-18, 32-cell groups, cached rate constants, four-step
  blocks, and dynamic scheduling: two RODAS3 points (rtol 1e-2 and 3e-3) for
  level 2, and three fourth-order linearly implicit Euler extrapolation points
  (rtol 4e-3 to 1e-3) for level 3. The recommended extrapolation point (rtol
  2e-3) reaches 3.45 digits overall and 2.55 for every species on the disclosed
  scenarios, 2.38 times faster (summed medians on 1,000 cells) than its own C
  three-stage Rosenbrock at rtol 1e-3 and atol 1e-18, the loosest setting it
  found to reach level 3.
- Other load during the session: Claude Code made two read-only progress checks.
  No builds, tests, or timing ran.
- Collected as `tier3/astra/xhigh`; transcript
  `transcripts/tier3/astra-xhigh.jsonl.gz`.
- Verified, 08:32–08:38, with `admin/verify_tier3.py`
  (`results/2026-10-18-tier3-astra-xhigh-verification`): both compilers build
  without warnings; 680 PASS, 0 FAIL, 0 N/A with each, as in Tier 1; the Tier 1
  TS1 gate passes (largest relative difference 2.03e-13 with clang, 2.14e-13
  with gcc-15); the audit finds nothing outside the rules. All five points pass
  the gate on the 6 disclosed and 5 held-out scenarios. fast and balanced reach
  level 2 on all 11, and accurate, recommended, and strict reach level 3 on all
  11. The disclosed scores match the agent's OPERATING_POINTS.md to the third
  decimal. fast's worst species reaches 1.035 digits, against a threshold of 1.

## Tier 3: astra-max

- Session start: 2026-10-18 08:39:53, launched by Claude Code with
  `admin/launch_codex_run.sh --tier3 astra-max gpt-6-astra max`; the goal was
  sent at 08:39:57. The session record shows gpt-6-astra at max, base
  instructions from gpt-6-astra, Codex 0.163.0, and no model switch; the `codex`
  command line carries `--search`.
- Codex session `01a14f74-7739-7a31-b3f8-1042c4a35f62`.
- Launched after astra-xhigh was collected, archived, and verified. The machine
  was on AC power, with load averages near 3 from the verification that had
  just ended (2.96 over 1 minute).
- The whole run was one turn, from the goal at 08:39:57 to its end at 11:14:18
  (2 h 34 min 21 s), with two context compactions, at 09:36 and 11:08. Every
  turn context shows gpt-6-astra at max. The session marked its goal complete
  and ended with DONE. Claude Code closed it with `/quit` at 11:26, after the
  record had been idle for 10 minutes; Codex's one warning was its usual
  startup notice about `-c` and `--search`.
- Tokens: 795,268 total; input 535,558 (plus 19,467,136 cached); output 259,710
  (176,408 reasoning). Codex's weekly limit went from 33 % used to 40 %.
- Web: 6 web items, all in its first seven minutes: 2 searches (5 queries), 3
  page openings, and 1 page read. It searched for Sandu et al. (1997), KPP's
  Rosenbrock methods, tolerances in atmospheric chemistry (Verwer, 1999),
  Alexander's (1977) DIRK methods, and Gustafsson's (1991) step-size control.
  It read a CWI report (`ir.cwi.nl/pub/125/0125D.pdf`), KPP's documentation of
  its Rosenbrock methods and its references page, TMAP8's page on an L-stable
  DIRK method, and the 2025 GMD paper on step-size control in KPP's Rosenbrock
  solvers; a PubMed Central page was refused by a CAPTCHA. REFERENCES.md cites
  every page it read.
- 6 commits after `0dc7f67`, from 08:50 to 11:14, ending `3cab602`, "Publish
  final work-precision results and reproducibility evidence": 87 files
  changed, 13,267 lines added and 232 removed. No remote; the working tree is
  clean.
- Compliance: no command reads outside its directory and the read-only inputs,
  apart from its check for AGENTS.md in the parent folders; no install or
  download command, no path into another run or the held-out folder, and no
  credential in the branch. Late in the run it searched its own files for model
  names and the home path, as a check for the blind review.
- Timer: `tools/ts1_scenario.c` times each step whole: the photolysis update,
  every rate constant (recomputed each step), and the integration. Setup,
  reading the input, a per-step scan for non-finite values, snapshot copies,
  and writing the output are outside it.
- Self-reported, not yet verified: both compilers pass the 680 cases and the
  Tier 1 gates; 960 research runs, 120 of them its final 10,000-cell timings;
  four RODAS3 points at rtol 3e-3 and 32-cell groups, differing in absolute
  tolerance and a cap on the internal step: economy for level 2, and fast3,
  balanced (recommended), and conservative for level 3. Balanced reaches at
  least 3.433 digits overall and 2.490 for every species on the six disclosed
  scenarios and its own three 48-hour cases.
- Other load during the session: Claude Code made read-only progress checks,
  and at 10:02 deleted untracked output files from other finished runs'
  directories (see "Disk space"). No builds, tests, or timing ran.
- Collected as `tier3/astra/max`; transcript
  `transcripts/tier3/astra-max.jsonl.gz`.
- Verified, 11:27–11:30 and then on the held-out scenarios, with
  `admin/verify_tier3.py` (`results/2026-10-18-tier3-astra-max-verification`):
  both compilers build (gcc-15 without warnings; clang with one, that it could
  not vectorize a loop marked `#pragma omp simd` in the error norm); 680 PASS,
  0 FAIL, 0 N/A with each; the Tier 1 TS1 gate passes (2.96e-13 with clang,
  2.75e-13 with gcc-15); the audit finds nothing outside the rules. All four
  points pass the gate on all 11 scenarios, and the disclosed scores match the
  agent's RESULTS.md to the third decimal. economy reaches level 2, and
  balanced and conservative level 3, on all 11. **fast3, declared for level 3,
  reaches it on 9 of 11**: it misses clean-48h-5min (2.46 overall, phenol 1.57)
  and winter60-10min (2.94, isoprene 1.87), the first declared point to miss
  its level on a held-out scenario.

## Tier 3: after the runs

- 2026-10-18, 11:14: astra-max, the last Tier 3 run, ended, so the held-out
  scenarios were published: `tier3/heldout/` holds `make_heldout.sh` (SHA-256
  `08016cac…195c`, as recorded above), the driver that made the data, the
  reference log, and checksums of the files, and each entry's verification
  folder has its held-out report.
- 2026-10-18, 11:30: `admin/time_entries.py --scenario` was tested with
  `admin/time_tier3.sh` on the six disclosed scenarios: astra-max's economy and
  MICM at rtol 1e-2 and atol 1e-15, 2 rounds at 1,000 cells. Names, gates, and
  scores were right, and the second round reused the first round's scores.
- The SSE moved the rest of Tier 3 to another machine: the MICM baselines, the
  timing session, the ranking, the blind review, and publishing. The steps,
  and the decisions already made, are in `tier3/HANDOFF.md`.

## Tier 3: the rest, on the same machine (2026-10-18 evening)

The SSE decided on 2026-10-18 at 19:32 to finish Tier 3 on the same laptop
after all, following `tier3/HANDOFF.md` step by step; no other machine was
set up. All times MDT.

- 19:40–19:55: **re-verification** of the four entries with
  `admin/verify_tier3.py`, both compilers, disclosed and held-out scenarios,
  into `~/NSF-NCAR/tier3-work/verify/`. Every table row of the regenerated
  reports (builds, the 680 cases, the TS1 gate, every point's digits, worst
  species, gate, and levels on all 11 scenarios) is identical to the
  committed reports apart from times; the clang binaries the timing uses
  are those builds.
- 19:43–20:03: **blind review**. The anonymized copies were made with
  `admin/anonymize_tier3.py` (letters drawn with `secrets`; the one
  leftover it found was bytes in a PNG), and two fresh Claude Code sessions
  were started from the review folder with `claude -p` at effort max, one
  with model opus reading A to D and one with model fable reading D to A,
  with the tools limited to Read, Grep, Glob, and Bash restricted to `wc`,
  `grep`, `head`, `tail`, and `cat`, and no web. Opus finished at 20:03:21
  after 159 turns, Fable at 20:01:06 after 180. The reports, the session
  records, the launch scripts, the key, and a README with the decoded
  scores and ranks are in `results/2026-10-18-tier3-blind-review/`.
  Decoded: A sol-xhigh, B astra-xhigh, C astra-max, D sol-max. Both
  reviewers rank astra-max first for code and sol-max last for research;
  Opus ranks astra-xhigh first for research and Fable astra-max.
- 20:13–: **MICM's baseline scan**, `admin/time_tier3.sh` with
  `admin/tier3_micm_scan.txt` (18 tolerance settings of the three-stage
  Rosenbrock, 32-cell groups) on the 11 scenarios at 1,000 cells, 10
  threads, one round, in a detached tmux session under `caffeinate`, into
  `~/NSF-NCAR/tier3-work/micm-scan/`. Nothing else ran on the machine
  during it. The baselines it picks go into the timing session's programs
  file (`~/NSF-NCAR/tier3-work/programs.txt`, built by `make_programs.sh`).
- 20:13–20:41: the **MICM scan** finished (198 runs; `ranking.md` in
  `results/2026-10-19-tier3-ranking/micm-scan/`). Every setting reaches
  level 2 on all 11 scenarios; the fastest is rtol 1e-2, atol 1e-15 (52.5 s
  summed, worst species isoprene 1.24 digits). Level 3 on all 11 is reached
  from rtol 3e-4 at atol 1e-15 (61.3 s; 3.29 digits, isoprene 2.05) and
  from 1e-3 at 1e-18 (70.3 s); the next setting, rtol 2e-4 at 1e-15
  (63.5 s, isoprene 2.21), is within one round's noise of it. MICM's time
  barely depends on the tolerance (52 to 97 s across the 18 settings)
  because its loop restarts every interval at 1e-6 of the step.
- 21:19: the **timing session** started (`~/NSF-NCAR/tier3-work/
  run-timing.sh` in tmux under `caffeinate`, with a load audit every
  minute): the 19 declared points and three MICM programs, `micm/L2`
  (rtol 1e-2, atol 1e-15), `micm/L3` (3e-4, 1e-15), and `micm/L3b` (2e-4,
  1e-15, so that a level-3 baseline survives if L3's 0.05-digit margin on
  isoprene does not hold at 10,000 cells; the ranking takes the fastest
  MICM point that reaches each level), on the 10,000-cell versions of the
  11 scenarios, 5 rounds, 10 threads: 1,210 runs. Between 20:41 and 21:19
  the machine built and tested the Fable entry (next item), which is why
  the session did not start at 20:41.
- 19:38 on: **Fable's Tier 3 reference solution**, commissioned by the SSE
  ("while the codex entrants work, you may also begin the Fable Tier 3
  reference solution, begin your own web search (use Sonnet researchers and
  review yourself, follow up yourself), then your structure scans";
  "structure scans are Fable at max effort"; goal "work until you judge you
  have arrived at the best algorithm, commit, push"). Set up like the
  contest runs: the Tier 1 slot `benchmarks/fable-max` received the
  history of TurboMICM's `fable/` (`git subtree split`, 47 commits, head
  6d2fcde), `tier3/prompts/make_prompts.py fable-max` wrote its prompt
  ("You are Fable."), and `tier3/admin/setup_runs.py fable-max` cloned it
  into `benchmarks/tier3/fable-max` with PROMPT.md (ODE 94ac532; the prompt
  mirrored in `tier3/prompts/`). Six Sonnet literature reviews ran in
  parallel (19:48–20:16) and are filed in the entry's `research/`. The
  entry is not a contestant; its collection branch will be
  `tier3/fable/max`.
- 05:46 on 2026-10-19: the **timing session ended**: 1,210 runs, no
  errors, every run passing the gate, every program's output byte-identical
  across its five rounds. `admin/rank_tier3.py` flagged nothing. The
  ranking (`results/2026-10-19-tier3-ranking/`), as sums over the 11
  scenarios of the median elapsed time: level 2, astra-xhigh/fast 34.72 s
  (13.4 times MICM at relative 1e-2, 465.40 s), sol-max/fast 39.80,
  sol-xhigh/L2 54.52, astra-max/economy 62.13; level 3, astra-xhigh/accurate
  74.47 s (7.3 times MICM at 3e-4, 542.66 s), sol-xhigh/L3 93.67,
  sol-max/accurate 106.57, astra-max/balanced 109.29. The held-out scenarios
  changed no rank; level 3's third place is a 2.6 % gap, inside the noise.
  Scores at 10,000 cells matched the 1,000-cell verification within 0.01
  digit; astra-max's fast3 misses level 3 on the same two held-out
  scenarios; `micm/L3` kept its isoprene margin, so `micm/L3b` was not
  needed. Load audit: AC power throughout; macOS's media analysis and
  Spotlight daemons ran at about half a core in a quarter of the minute
  samples, spread over all programs by the rotation; round-to-round spreads
  were 1.3 to 3.9 % on average and 8.4 % at most. Published on the wiki
  (Tier-3-Results, Home); `tier3/HANDOFF.md` retired, its rules moved to
  `tier3/README.md`.
- 05:52: the Fable entry's S4 scan (controllers) started, the first work on
  the machine after the session.

## Fable's Tier 3 reference solution: scans and operating points (2026-10-19 morning)

Continuation of the 2026-10-18 evening item above, in
`~/NSF-NCAR/ODE/benchmarks/tier3/fable-max/` (collected afterwards as
`tier3/fable/max`). Every scan ran with `tools/scan.py` on 1,000-cell
scenarios at 10 threads with nothing else timing on the machine; the entry's
SCANS.md has the tables and DECISIONS.md the choices.

- 20:41–21:19 on 2026-10-18 (before the timing session): S1 (the first
  step: the carried, quantized step is 35 % cheaper than MICM's restart),
  S2 (ten methods against tolerance on the six scenarios: RODAS3 leads to
  three digits, extrapolation beyond), S6a (factorization reuse: not a
  lever, the photolysis changes the matrix every daytime interval).
- 05:52–06:02: S4, the step-size controller (106 points on diurnal-10min
  and polluted-10min): H211b with b = 4 is the one controller above the
  elementary controller's work-precision curve.
- 06:02–06:08: S5 (the error norm: a tolerance rescaling, MICM's group RMS
  kept), S3 (the absolute tolerance: a lever at level 3, since isoprene
  and phenol sit near the 1e-14 floor the score counts from), S7
  (clamping: no effect, kept off), S9 (grouping at 10 threads: 32 kept).
- 06:07–06:12: S4b, the H211b grid over safety and tolerance: 30 to 35 %
  less time than the elementary controller above 1.5 digits per species;
  06:12–06:19: S8 (extrapolation sequences: the shifted sequence 2, 3, 4,
  ... beats the harmonic one; orders above 4 useless; the even sequence
  fails at orders 5 and 6) and S9 at one thread (32 confirmed).
- 06:19–06:27: S12, the joint grid rtol x atol x controller under the
  final conventions (quantized step), and S4c, the loose end of the H211b
  curve: the S4b floor at rtol 1e-2 was the absolute tolerance's.
- 06:27–06:39: S13, 21 candidates on all six disclosed scenarios, three
  rounds; 06:39–06:43: S13b, the loose H211b points for level 2 (all meet
  a 2.04-digit floor on diurnal-5min, BrO 1.15, the contest's floor:
  one step per 5-minute interval, so the controller cannot act).
- 06:47–06:53: S10 (generated against generic kernels at the final points:
  5.3 to 5.9 times, identical steps and digits, roundoff-different outputs)
  and S14 (the final points at one thread: 10 threads buy 4.6 to 5.1
  times; the kernel cost model explains about three quarters of the time).
- 06:53–07:01: MICM's three-stage Rosenbrock at the ranking's baseline
  settings on the six scenarios, three rounds (`tools/micm_baseline.py`;
  two reruns after a path and a snapshot-file fix in the runner): 24.6 s at
  level 2 and 29.0 s at level 3 against the five points' 1.16 to 3.45 s.
- **Operating points** (OPERATING_POINTS.md, with ACCURACY.md's rule of
  0.2 digit of margin for the worst species on every disclosed scenario):
  L2, RODAS3 with MICM's elementary controller at rtol 5e-3 and atol 1e-20
  (2.42 digits, isoprene 1.29, 1.16 s over the six scenarios, 21 times
  MICM at level 2); balanced, RODAS3 with H211b at 2e-3 and 1e-20
  (3.17 / 2.17, 2.21 s, declared for level 2, three digits with 0.17 of
  margin); L3, extrapolation of order 3 with the shifted sequence at 1e-3
  and 1e-20 (3.61 / 2.49, 2.65 s, 11 times MICM at level 3); recommended,
  order 4 at 5e-4 and 1e-22 (3.62 / 2.69, 3.04 s, 9.5 times); accurate,
  order 4 at 2e-4 and 1e-20 (3.87 / 3.00, 3.45 s, 8.4 times). The loose
  RODAS3 points at rtol 1e-2 reach the contest's 2.04-digit floor on
  diurnal-5min and are reported, not declared.
- 07:04–07:33: the gcc-15 and clang gates with build timings; two gcc
  warnings fixed (a const-qualified array parameter, a truncating snprintf)
  and a Makefile race (an object compiled in the same second as the flags
  file's rewrite passed as current after a compiler switch; the rule now
  removes every object when the flags change). Final commit `ec34c99`
  (26 commits after the Tier 1 head `6d2fcde`), collected as
  `tier3/fable/max` and pushed at 07:34.
- 07:34–07:44: **verification** with `admin/verify_tier3.py`, both
  compilers, disclosed then held-out
  (`results/2026-10-19-tier3-fable-max-verification/`): no warnings; 673
  PASS, 0 FAIL, 7 N/A on the 680 ids with each compiler; TS1 gate 1.04e-13
  and 9.97e-14; every point passes the gate on all 11 scenarios and reaches
  its declared level on all 11 (L2's thinnest held-out margin: isoprene 1.16
  on winter60-10min); balanced reaches three digits on 10 of 11 (not on
  clean-48h-5min: 2.71, phenol 1.79), as its level-2 declaration
  anticipated. Published: the wiki page Fable-Tier-3-Reference-Solution, a
  section on Tier-3-Results, Home; CLAUDE.md's Tier 3 notes and open items.
- 07:35: the SSE's other agent renamed the task force on the wiki (Agentic
  Coding Task Force); pulled before the wiki commit, and CLAUDE.md renamed
  to match at the SSE's request.
- 07:49–08:01: the goal's Stop hook asked for the judgment that the best
  algorithm had been reached; a **closing scan** (S15) tried what the
  search had left open, ROS3 and ROS4 under H211b and the tolerances between
  the Rosenbrock points, on two scenarios, then ROS3 on all six in three
  rounds (S15b). ROS3 under H211b at rtol 8e-3 and atol 1e-20 is the
  cheapest configuration meeting the level-3 margin rule (3.27 digits,
  isoprene 2.37, 2.34 s over the six, 12 % under the extrapolation point);
  it replaced the balanced point as **L3r**, declared for level 3. S16
  measured it with generic kernels and at one thread. The judgment is
  written in the entry's DECISIONS.md. Final commit `cb0f8a3` (28 commits
  after the Tier 1 head), re-collected as `tier3/fable/max` at 08:09 and
  re-verified (next line).
- 08:04–08:13: **re-verification** of `cb0f8a3` with both compilers on the
  disclosed and held-out scenarios: no warnings; 673 PASS, 7 N/A; TS1 gate
  1.04e-13 and 9.97e-14; all five points at their declared levels on all 11
  scenarios, L3r's held-out lowest 3.18 digits (HYDRALD 2.42 on
  tropics-20min). The results folder now holds this verification; the first
  one, with the balanced point, is recorded above.
- 08:15–08:19: a **comparison session** was started
  (`~/NSF-NCAR/tier3-work/run-timing-fable.sh`: Fable's five points, the
  contest's two fastest level-2 and two fastest level-3 points, and MICM's
  two baselines on the 10,000-cell scenarios, 3 rounds, 10 threads, 363
  runs, about 1 h 45 min) and stopped after five runs at the SSE's request,
  deferred to a follow-up session. To run it: that launcher, with
  `~/NSF-NCAR/tier3-work/programs-fable-vs-contest.txt` (Fable's binary from
  the verification worktree `tier3-work/verify/fable-max-2/clang/fable-max/`),
  then `admin/rank_tier3.py <out> --threads 10`; idle machine, AC power.
- 08:25: a **provisional placement** from the verification runs alone
  (single runs, 1,000 cells, 10 threads; the contest's points from their
  verifications of 2026-10-18, Fable's from 2026-10-19, so different
  sessions): at level 2, L2 1,429 ms on the five held-out scenarios against
  astra-xhigh/fast 2,110 and sol-max/fast 2,230 (disclosed: 1,196 against
  1,681 and 1,914); at level 3, L3r 2,942 ms against astra-xhigh/accurate
  4,519 and sol-xhigh/L3 5,551 (disclosed: 2,358 against 3,608 and 4,730),
  about 1.4 to 1.6 times faster at each level. The deferred session is the
  measurement; this is the estimate, on the wiki page with its caveat.
- 11:42–14:50: the **comparison session**, resumed at the SSE's request
  ("resume our deferred Tier 3 work") in a fresh Claude Code session. The
  laptop was on battery after a reboot, so the launcher was armed at 11:42
  behind a gate that waited for the adapter
  (`~/NSF-NCAR/tier3-work/run-timing-fable-on-ac.sh`); the adapter was
  attached at 12:33 and the session ran from 12:34:06 to 14:50:42: the 11
  programs of `programs-fable-vs-contest.txt` (Fable's five points from the
  `cb0f8a3` build in `verify/fable-max-2/`, astra-xhigh/fast, sol-max/fast,
  astra-xhigh/accurate, sol-xhigh/L3, and MICM's two baselines) on the
  10,000-cell versions of the 11 scenarios, 3 rounds, 10 threads, 363 runs,
  every one through the gate, every output deterministic
  (`results/2026-10-19-tier3-fable-vs-contest/`). Sums over the 11
  scenarios of the median elapsed time: at level 2, L2 24.55 s against
  astra-xhigh/fast 34.10 and sol-max/fast 39.45 (MICM 464.74), 1.39 and 1.61
  times faster and 18.9 times MICM; at level 3, L3r 49.08 s against
  astra-xhigh/accurate 73.46 and sol-xhigh/L3 95.69 (MICM 549.07), 1.50 and
  1.95 times faster and 11.2 times MICM; L3 54.60, recommended 62.14,
  accurate 71.91. L2 clears level 2 by 0.27 digit overall and 0.16 per
  species (the winners by 0.04 to 0.06); L3r clears level 3 by 0.18 and
  0.33 (astra-xhigh/accurate by 0.20 and 0.39). The contest's points and
  MICM's baselines reproduced their ranking-session sums within 2.2 % with
  identical scores, so the contest's ranking stands as published. Spread
  2.7 % of the median on average, 14.3 % at most; the gaps between adjacent
  places 16 % or more. Load audit: 137 samples, AC attached throughout, the
  usual Spotlight and media-analysis daemons, this session's status checks.
  Published: the measured table replaces the provisional one on
  Fable-Tier-3-Reference-Solution; the sentence on Tier-3-Results; Home.
  `tier3/FOLLOW-UP.md` deleted; CLAUDE.md's open item closed.
- Transcripts: this session and the six research sessions are on the
  `transcripts` branch under `transcripts/tier3/fable-max*`.

# Disk space

- 2026-10-18, 10:02, while astra-max ran, on the SSE's go-ahead to clean up the
  output of completed and recorded runs: Claude Code deleted 1,742 untracked CSV
  files of 1 MiB or more, 19.7 GiB, from the run directories of five finished
  Tier 1 runs (sol-high, sol-xhigh, sol-max, astra-high, astra-xhigh) and three
  finished Tier 3 runs (sol-xhigh, sol-max, astra-xhigh). 14.8 GiB of it was
  Tier 3 astra-xhigh's scan output, in `build/scan-runs/`. The files were solver
  outputs (`ts1_final.csv`, `ts1_snapshots.csv`, and the like) and generated cell
  inputs (`cells.csv`), none tracked by git. Every tracked file and log was kept,
  and so were the agents' own reference and check files (`reference_*.csv`,
  `check_*.csv`) from their stress and scaling tests. With Claude Code's scratch
  verification outputs and gcc-15 builds also deleted, free space went from 19 to
  45 GiB. Tier 1 astra-max's directory waited until Tier 3 astra-max ended.
- 2026-10-18, 11:58, after Tier 3 astra-max was collected, archived, and
  verified: the same rule deleted 158 files, 2.70 GiB, from Tier 3 astra-max's
  run directory (1.64 GiB) and Tier 1 astra-max's (1.06 GiB).
