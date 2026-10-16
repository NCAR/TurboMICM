# Session transcripts

Each finished run's full session transcript, as the agent's harness saved it,
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
