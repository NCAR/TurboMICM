# Run log

Launch record for the MICM-to-C runs. Times are local (MDT).

Each run is an interactive session started in its run directory,
`benchmarks/<run>/`, and launched with one instruction:

    /goal read PROMPT.md and work until done, commit, push

The run repositories have no remote, so "push" has nowhere to go.

A run is `<model>-<effort>`; repeat instances of the same model and effort, run as
time permits, add a suffix: `sol-high-2`, `sol-high-3`. Each finished run is kept
on the TurboMICM branch `runs/<model>/<effort>`, or `runs/<model>/<effort>-<n>`
for a repeat.

Contestant changes
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

Reference changes
- 2026-10-15: `reference/micm_ts1_ref.cpp` gained `--group-size G`, so MICM can be
  timed with the same cell grouping as an entry. Its default, one group per
  thread, is unchanged and gives bit-identical results.

Verification
- 2026-10-15: the three Sol entries were rebuilt from their branches with clang
  and gcc-15, their tests rerun and counted against the case list, their builds
  audited against the portability rules, and their TS1 timings replicated with
  `admin/time_entries.py` alongside MICM. The results are in
  `results/2026-10-15-sol-verification/` and on the wiki's Tier-1-Results page.

Transcripts
- Each finished run's session transcript is archived on the branch `transcripts`,
  in `transcripts/<run>.jsonl.gz`, kept apart from `main` and the run branches
  because a transcript names the model and shows every step. The files are made
  with `admin/archive_transcript.py`, which redacts the account identifiers in
  the session header and gzips the rest unchanged.
- 2026-10-15: sol-high, sol-xhigh, and sol-max archived. Their reasoning is
  encrypted by OpenAI and cannot be read; messages, tool calls with their full
  output, and token counts are complete.

| Run | Model | Effort | Harness | Session | Final commit | Status |
|---|---|---|---|---|---|---|
| sol-high | gpt-6.1-sol | high | Codex CLI 0.162.1, interactive | 2026-10-14 21:55-23:32 | `1d944f4`, 23:32 (6 commits) | finished |
| sol-xhigh | gpt-6.1-sol | xhigh | Codex CLI 0.162.1, interactive | 2026-10-15 03:54-05:29 | `171f54a`, 05:29 (6 commits) | finished |
| sol-max | gpt-6.1-sol | max | Codex CLI 0.162.1, interactive | 2026-10-15 07:26-09:13 | `001fb64`, 09:13 (7 commits) | finished |
| astra-high | gpt-6-astra | high | Codex CLI 0.162.1, interactive | 2026-10-15 12:10- | | running |

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
  `benchmarks/astra-high`, with `~/.codex/config.toml` set to gpt-6-astra at high
  effort; the session record confirms both.
- Codex session `01a140c1-b679-7082-aef7-62cd2c85941c`.
- Other load during the session: Tier 2 preparation ran from 13:24 to 13:34, using
  at most two threads at nice 10–15. It included two compiles, a few seconds of
  test runs, and the start of the reference computation. Astra ran no timing in
  that window; it was porting the sparse-matrix tests. At 13:34 the reference
  computation was suspended (SIGSTOP) until the run finished, to keep the machine
  idle. No other work from the experiment ran during the Sol sessions.
