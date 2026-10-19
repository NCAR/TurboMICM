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
