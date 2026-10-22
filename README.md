# Session transcripts

Each finished contest run's full session transcript, as the agent's harness saved it,
compressed with gzip. They are kept on this branch, apart from `main` and the
run branches, for two reasons: they name the model, so they must stay away from
the blind review, and they show every step of a solution, so they must stay away
from agents that are still running.

| File | Run | Harness | Session id | Original size |
|---|---|---|---|---:|
| `transcripts/sol-high.jsonl.gz` | sol-high | Codex CLI 0.162.1 | `01a13db3-1129-76d1-8797-c3b73194043b` | 10.0 MB |
| `transcripts/sol-xhigh.jsonl.gz` | sol-xhigh | Codex CLI 0.162.1 | `01a13efc-1abd-7381-a20e-142540e76a25` | 7.7 MB |
| `transcripts/sol-max.jsonl.gz` | sol-max | Codex CLI 0.162.1 | `01a13fbc-d221-7b81-88a8-3d73454dd2b7` | 7.4 MB |
| `transcripts/astra-high.jsonl.gz` | astra-high | Codex CLI 0.162.1 | `01a140c1-b679-7082-aef7-62cd2c85941c` | 10.5 MB |
| `transcripts/astra-xhigh.jsonl.gz` | astra-xhigh | Codex CLI 0.163.0 | `01a142b6-b159-7d52-801a-7ac70fa3fed4` | 9.9 MB |
| `transcripts/astra-max.jsonl.gz` | astra-max | Codex CLI 0.163.0 | `01a14386-84c0-70a3-867a-6cf69b22b824` | 7.1 MB |

## Tier 2 runs

Each Tier 2 run starts from a copy of the same model's Tier 1 repository, without
web search, on a second machine, an Apple M5 Pro (Tiers 1 and 3 ran on an Apple
M5). The file format and the redaction are those described under Contest runs,
below.

| File | Run | Harness | Session id | Original size |
|---|---|---|---|---:|
| `transcripts/tier2/astra-max.jsonl.gz` | Tier 2 astra-max | Codex CLI 0.163.0 | `01a163e7-f388-7153-baf9-d8a8478cce81` | 1.8 MB |
| `transcripts/tier2/astra-xhigh.jsonl.gz` | Tier 2 astra-xhigh | Codex CLI 0.163.0 | `01a1645d-8d6b-7ea1-adde-c65151f4d896` | 2.7 MB |
| `transcripts/tier2/astra-high.jsonl.gz` | Tier 2 astra-high | Codex CLI 0.163.0 | `01a164c1-c86d-74f1-a385-74d24b64d3a1` | 2.3 MB |
| `transcripts/tier2/sol-max.jsonl.gz` | Tier 2 sol-max | Codex CLI 0.163.0 | `01a16520-0dd1-7032-bd32-99fef62218eb` | 4.4 MB |

## Tier 3 runs

Each Tier 3 run starts from a clone of the same model's Tier 1 repository and
runs with web search on (`codex --search`). The file format and the redaction
are those described under Contest runs, below.

| File | Run | Harness | Session id | Original size |
|---|---|---|---|---:|
| `transcripts/tier3/sol-xhigh.jsonl.gz` | Tier 3 sol-xhigh | Codex CLI 0.163.0 | `01a14df8-cf84-7851-b9c3-3c76c4eef53e` | 8.2 MB |
| `transcripts/tier3/sol-max.jsonl.gz` | Tier 3 sol-max | Codex CLI 0.163.0 | `01a14e89-d573-77d2-8286-e1d020bdb20d` | 5.9 MB |
| `transcripts/tier3/astra-xhigh.jsonl.gz` | Tier 3 astra-xhigh | Codex CLI 0.163.0 | `01a14ef4-db46-7743-8cf5-57c66c1aa74c` | 7.5 MB |
| `transcripts/tier3/astra-max.jsonl.gz` | Tier 3 astra-max | Codex CLI 0.163.0 | `01a14f74-7739-7a31-b3f8-1042c4a35f62` | 7.8 MB |

## Fable's reference entry (not a contest run)

The Claude Code session in which Claude Fable 5.1 built its own solution to
the Tier 1 prompt outside the contest (wiki: Fable's reference entry), and the
nine subagents it dispatched to port MICM's tests. Nothing is redacted: the
session is not part of the contest, and the SSE chose to keep it whole. The
session ran from 2026-10-16 23:14 to the morning of 2026-10-17 (local time);
its id is `d57cf5ee-4422-4a52-8b7f-2a60d5ce7ae2`.

| File | Role | Harness | Original size |
|---|---|---|---:|
| `transcripts/fable-reference-entry.jsonl.gz` | the main session, to 08:45 on 2026-10-17 | Claude Code 2.1.301 | 44.8 MB |
| `transcripts/fable-reference-entry-agents/port-dense.jsonl.gz` | matrix, vector_matrix, grouped_views, view_category_sfinae tests | Claude Code subagent | 7.2 MB |
| `transcripts/fable-reference-entry-agents/port-solver.jsonl.gz` | rosenbrock, backward_euler, state, linear_solver, lu_decomposition, jacobian_verification, chapman_ode_solver tests | Claude Code subagent | 10.3 MB |
| `transcripts/fable-reference-entry-agents/port-constraints.jsonl.gz` | constraint and DAE tests | Claude Code subagent | 9.9 MB |
| `transcripts/fable-reference-entry-agents/port-external.jsonl.gz` | external model and aerosol model tests | Claude Code subagent | 7.3 MB |
| `transcripts/fable-reference-entry-agents/port-rates.jsonl.gz` | rate-constant, process, system, and solver-builder tests | Claude Code subagent | 9.9 MB |
| `transcripts/fable-reference-entry-agents/port-sparse.jsonl.gz` | sparse matrix tests; stalled after an hour and was replaced | Claude Code subagent | 1.2 MB |
| `transcripts/fable-reference-entry-agents/port-integration.jsonl.gz` | integration, regression, and tutorial tests; stalled and was replaced | Claude Code subagent | 4.7 MB |
| `transcripts/fable-reference-entry-agents/port-sparse-2.jsonl.gz` | sparse matrix tests (replacement) | Claude Code subagent | 3.1 MB |
| `transcripts/fable-reference-entry-agents/port-integration-2.jsonl.gz` | integration, regression, and tutorial tests (replacement) | Claude Code subagent | 3.9 MB |

Claude Code's transcripts hold every message, tool call, and tool output as
JSON lines, with the model's reasoning as the harness records it; the
subagents' briefs are their first user messages. Each subagent's
`.meta.json` sits beside its transcript.

## Tier 3: Fable's reference solution (not a contest run)

The Claude Code session in which Claude Fable 5.1 built its own solution to
the Tier 3 prompt outside the contest (branch `tier3/fable/max`; wiki:
Tier-3-Results), and the six literature-research sessions it started on
2026-10-18 at 19:48 local time with Claude Sonnet 5.5, whose reports are
`research/R1` to `R6` in the entry. Nothing is redacted, as for the Tier 1
reference entry. The main session also carried the contest's own Tier 3
work that evening and night (the re-verification, the blind review, the
timing session, the ranking, and their publication), since the SSE asked
for the reference solution while the timing session ran; its id is
`3b43ec7e-eef6-4039-b6da-5097df839d0b`.

| File | Role | Harness | Original size |
|---|---|---|---:|
| `transcripts/tier3/fable-max.jsonl.gz` | the main session, from 19:31 on 2026-10-18 to its archiving at 08:22 on 2026-10-19 | Claude Code 2.1.303 | 31.1 MB |
| `transcripts/tier3/fable-max-agents/R1-rosenbrock-w-methods.jsonl.gz` | brief R1: Rosenbrock and Rosenbrock-W methods | Claude Code session, Claude Sonnet 5.5 | 6.1 MB |
| `transcripts/tier3/fable-max-agents/R2-accuracy-requirements.jsonl.gz` | brief R2: accuracy requirements of atmospheric chemistry solvers | Claude Code session, Claude Sonnet 5.5 | 18.9 MB |
| `transcripts/tier3/fable-max-agents/R3-extrapolation-exponential-sdirk.jsonl.gz` | brief R3: extrapolation, exponential, and SDIRK methods | Claude Code session, Claude Sonnet 5.5 | 7.0 MB |
| `transcripts/tier3/fable-max-agents/R4-controllers-norms-reuse-positivity.jsonl.gz` | brief R4: step-size controllers, error norms, factorization reuse, positivity | Claude Code session, Claude Sonnet 5.5 | 10.4 MB |
| `transcripts/tier3/fable-max-agents/R5-alternative-schemes-precision.jsonl.gz` | brief R5: alternative schemes and reduced precision | Claude Code session, Claude Sonnet 5.5 | 8.7 MB |
| `transcripts/tier3/fable-max-agents/R6-ts1-camchem-musica-sda.jsonl.gz` | brief R6: TS1 in CAM-chem and MUSICA, and the solver's use there | Claude Code session, Claude Sonnet 5.5 | 15.1 MB |

Each research session's `.meta.json` beside it has the session id, the
model, and the start and end times (UTC); the sessions ran in parallel from
01:48 to 02:17 UTC on 2026-10-19 (19:48 to 20:17 local on 2026-10-18).

## Contest runs

Each line is one JSON event:
- the session settings: model, effort, and the harness's own instructions;
- every message;
- every tool call with its full output;
- the token counts after each turn.

The models' reasoning is stored only as encrypted content, which cannot be read.

`admin/archive_transcript.py` on `main` made these files. It replaced the two
OpenAI account identifiers in each session header (`creator_account_id`,
`creator_user_id`) with `<redacted>` and left every other line as Codex wrote it.
A scan for credentials found none in readable text; its only pattern matches were
inside the encrypted reasoning.

To read one:

    gunzip -c transcripts/sol-max.jsonl.gz | head
