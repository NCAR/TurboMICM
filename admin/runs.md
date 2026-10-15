# Run log

Launch record for the MICM-to-C runs. Times are local (MDT).

Each run is an interactive session started in its run directory,
`benchmarks/<run>/`, and launched with one instruction:

    /goal read PROMPT.md and work until done, commit, push

The run repositories have no remote, so "push" has nowhere to go.

| Run | Model | Effort | Harness | Session | Final commit | Status |
|---|---|---|---|---|---|---|
| sol-high | gpt-6.1-sol | high | Codex CLI 0.162.1, interactive | 2026-10-14 21:55-23:32 | `1d944f4`, 23:32 (6 commits) | finished |
| sol-xhigh | gpt-6.1-sol | xhigh | Codex CLI 0.162.1, interactive | 2026-10-15 03:54-05:29 | `171f54a`, 05:29 (6 commits) | finished |
| sol-max | gpt-6.1-sol | max | Codex CLI 0.162.1, interactive | 2026-10-15 07:26-09:13 | `001fb64`, 09:13 (7 commits) | finished |

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
