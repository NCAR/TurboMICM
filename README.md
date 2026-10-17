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

## Fable's reference entry (not a contest run)

The Claude Code session in which Claude Fable 5.1 built its own solution to
the Tier 1 prompt outside the contest (wiki: Fable's reference entry), and the
nine subagents it dispatched to port MICM's tests. Nothing is redacted: the
session is not part of the contest, and the SSE chose to keep it whole. The
session ran from 2026-10-16 23:14 to the morning of 2026-10-17 (local time);
its id is `d57cf5ee-4422-4a52-8b7f-2a60d5ce7ae2`.

| File | Role | Harness | Original size |
|---|---|---|---:|
| `transcripts/fable-reference-entry.jsonl.gz` | the main session | Claude Code 2.1.301 | 38.8 MB |
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
