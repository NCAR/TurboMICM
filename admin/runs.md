# Run log

Launch record for the MICM-to-C runs. Times are local (MDT). Tier 1 comes first;
the Tier 3 runs follow, under [Tier 3](#tier-3).

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
  the session's context was compacted twice along the way. Every turn context
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
| astra-xhigh | gpt-6-astra | xhigh | `fbd365c` | | | | queued |
| astra-max | gpt-6-astra | max | `0dc7f67` | | | | queued |

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
