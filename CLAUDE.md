# TurboMICM

Controlled experiment, for the NSF NCAR Agentic AI Task Force: LLM coding agents each
port MICM's CPU solver to portable C (C11, libm, OpenMP) from the same prompt. Each
entry is judged by:
- MICM's 680 test cases;
- TS1 agreement with MICM to 1e-10 relative;
- TS1 wall time;
- a blind code review.

The public record is the wiki: NCAR/TurboMICM.wiki, checked out at
`~/NSF-NCAR/TurboMICM.wiki`, branch `master`. Update it as runs come in.

## Where things live

- **Run root, the prompts' ROOT:** `~/NSF-NCAR/ODE/benchmarks`, in the ODE repo on
  branch `micm-c-port-experiment`. It holds:
  - `micm/`: MICM e15898ce, with a Release, double-precision, OpenMP build;
  - `reference/`, `prompts/`, `admin/`;
  - one git repo per run, `<model>-<effort>/`, with no remote.

  `micm/` and the run repos are gitignored.
- **This repo:**
  - `main`: `prompt.v3` (v1 and v2 for history), `prompts/`, `reference/`, `admin/`
    (scripts and `runs.md`, the run log), and `results/<date>-<topic>/`;
  - `runs/<model>/<effort>[-<n>]`: each finished run's repo, fetched unchanged;
  - `transcripts`: an orphan branch of gzipped, redacted session transcripts.
- **Keep the two copies identical:** `prompt.v3`, `prompts/`, `reference/`, and the
  prompt and setup scripts exist in both the run root and this repo. Check with
  `diff -r`.
- **Contestants:**
  - Codex CLI 0.162.1: gpt-6.1-sol, gpt-6-astra, and gpt-6-luna. Luna replaced
    gpt-5.6-terra on 2026-10-15.
  - Claude Code: Sonnet, Opus, Fable.
  - Effort levels: high, xhigh, max.
- **Run names:** `<model>-<effort>`; a repeat is `<model>-<effort>-<n>`.
- **Codex files:**
  - sessions: `~/.codex/sessions/YYYY/MM/DD/rollout-*-<id>.jsonl`;
  - config: `~/.codex/config.toml`;
  - model list: `~/.codex/models_cache.json`.

## Procedures

- **Prompts:** in the run root, run `python3 prompts/make_prompts.py [--instance N RUN ...]`,
  then `python3 admin/setup_runs.py [RUN ...]`. Copy any changes to this repo.
  - `setup_runs.py` refuses unfilled placeholders.
  - It never touches a started run or rewrites an unchanged PROMPT.md, so it is safe
    while a run is live.
- **Launch:** the SSE launches every run interactively, on subscriptions; never the
  API or `codex exec`.
  - Start Codex with
    `cd <run> && codex -m <model> -c model_reasoning_effort=<effort> -c features.memories=false`.
  - Then give the session `/goal read PROMPT.md and work until done, commit`.
  - Always pass the effort: Codex's model picker saves the last choice in config.toml.
- **Record:** take each run's start, end, model, effort, and tokens from its session
  JSONL, then add them to `admin/runs.md`.
  - Tokens are the last `token_count` event's `info.total_token_usage`.
  - total = input − cached + output.
- **Collect:**
  `git fetch ~/NSF-NCAR/ODE/benchmarks/<run> main:runs/<model>/<effort> && git push origin runs/<model>/<effort>`.
- **Transcript:**
  - Run `admin/archive_transcript.py <session.jsonl> transcripts/<run>.jsonl.gz`;
    it redacts the account IDs and keeps every other line.
  - Commit the file to the `transcripts` branch from a temporary worktree, and push.
  - In a fresh orphan checkout, clear the index with `git rm -rfq .`.
- **Verify an entry:**
  - Make a detached worktree of `runs/*` in a directory whose `micm/` and
    `reference/` symlink to the run root's. The entries' Makefiles use `..`.
  - For each of `CC=clang` and `CC=gcc-15`: run `make env`, then
    `make cpu bench OPTFLAGS="-O3 -DNDEBUG"`, then `make test-cpu`.
  - Count the `<id> PASS|FAIL|N/A` lines against `reference/micm_cases.txt`: 680 ids,
    each exactly once.
  - Audit the flags (`make -n -B`) and the linked libraries (`otool -L`), and grep the
    sources for banned constructs.
- **MICM reference:**
  `clang++ -std=c++20 -O3 -DNDEBUG -I micm/include -I micm/test/benchmark reference/micm_ts1_ref.cpp -o micm_ts1_ref`.
  - The default is one cell group per thread.
  - `--group-size 32` matches the entries' grouping and solver work.
- **Timing:**
  `admin/time_entries.py --input cells.csv --reference micm_final.csv --out DIR micm=… "micm-g32=… --group-size 32" <run>=<dir>/bin/ts1_bench_cpu …`.
  - It runs 5 rotating rounds at 10 and 1 threads, and checks every run against MICM.
  - Run it on an idle machine on AC power, with no agent running.
- **Inputs:** `reference/make_cells_csv.py 10000 <seed>` writes the cells file;
  `reference/compare_final.py YOURS MICM` applies the gate.

## Conventions

- **Results:** the wiki shows verified numbers only, with tables updated in place.
  Agents' self-reports appear only as a one-line finding.
- **Naming:** call the user "the SSE" (scientific software engineer). Name them once,
  in the preamble of Role-of-Claude-in-the-Design, and use no pronouns for them.
- **Conflict of interest:** Claude models are contestants, so Claude must not grade or
  review entries alone. Keep the disclosure on the Role page.
- **Git:**
  - Commits end with the Co-Authored-By trailer.
  - Push this repo and the wiki; the SSE wants them pushed for backup.
  - Push the run-root branch too, after its files change.

## Pitfalls

- **Placement of agent instructions:** never put CLAUDE.md or AGENTS.md in the run
  root or any folder above it (`~/NSF-NCAR/ODE`, `~/NSF-NCAR`, `~`), because the
  harnesses load them into contestants' sessions. The `.claude/settings.local.json`
  files in `~/NSF-NCAR/ODE` and `~/NSF-NCAR`, and the user plugins in `~/.claude`,
  already apply to Claude runs.
- **Isolation is honor-system:** agents run with full access and the SSE's GitHub and
  SSH credentials. Keep solutions and transcripts off `main`'s working tree, and
  never edit a running session's directory.
- **Clang linking:** Apple clang links only with
  `SDKROOT=$(xcrun --sdk macosx --show-sdk-path)`. The default Command Line Tools 27
  SDK lists an `arm64e.x1` architecture that Xcode's linker rejects.
- **OpenMP with Apple clang:**
  `-Xpreprocessor -fopenmp -isystem /opt/homebrew/opt/libomp/include … -L/opt/homebrew/opt/libomp/lib -lomp`.
  With `-I` instead of `-isystem`, `omp.h` fails `-pedantic`.
- **MICM build:** cmake needs `-DFETCHCONTENT_TRY_FIND_PACKAGE_MODE=NEVER`, because
  conda's GTest breaks the tests, plus explicit `OpenMP_*` hints for libomp.
- **MICM tests:** ctest records absolute paths, so run the test executables from a
  copy of `build/`.
- **Shell quirks:**
  - `ls` and `grep` are aliased to `eza` and `ugrep`; use `command ls` and
    `command grep` in pipelines, where ugrep can hang.
  - There is no `timeout` command.
  - zsh reads `$var:c` as a history modifier; write `${var}:` instead.
- **Timing drift:** MICM on one thread took 137.4 s on 2026-10-14 and 144.6 s on
  2026-10-15. Compare times only within one interleaved session.
  - 10-thread runs vary 5–10 % on the M5's mixed cores; 1-thread runs vary 0.3 % or
    less.
- **Transcripts:** Codex stores reasoning only as encrypted content. The transcripts
  have every message, tool call, and output, but no readable reasoning.

## Open items

- Keep the SSE's Claude Code setup (plugins, hooks, memory, `.claude` settings) out
  of the Claude runs.
- Rank with an undisclosed seed, and run the blind review.
- Consider Claude Haiku as the fast-tier match for Luna.
- Time the entries on an x86-64 Linux machine.
- Current status is in `admin/runs.md` and on the wiki's Results page.
