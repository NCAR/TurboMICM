# TurboMICM

Controlled experiment, for the NSF NCAR Agentic Coding Task Force: OpenAI's coding agents
(Codex CLI) each port MICM's CPU solver to portable C (C11, libm, OpenMP) from the
same prompt. Anthropic's Claude Opus and Claude Fable, at max reasoning effort,
design the experiment, orchestrate the runs, and validate, evaluate, and judge the
entries; afterwards they build on the results a framework usable from MUSICA,
CheMPAS-A, CAM-chem, and other host models. Each entry is judged by:
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
  - `derecho/` on `main`: the Fable entry tuned for NCAR's Derecho (gcc 14, AMD
    EPYC 7763), with its job scripts, results, and README; merged 2026-10-18 from
    `derecho/fable` and `derecho/fable-2`. Wiki: Fable-on-Derecho. New Derecho
    work goes on a branch off `main` under `derecho/`, never in `fable/`;
  - `runs/<model>/<effort>[-<n>]`: each finished run's repo, fetched unchanged;
  - `transcripts`: an orphan branch of gzipped, redacted session transcripts.
- **Keep the two copies identical:** `prompt.v3`, `prompts/`, `reference/`, and the
  prompt and setup scripts exist in both the run root and this repo. Check with
  `diff -r`. `tier3/` likewise mirrors the committed files of
  `ODE/benchmarks/tier3` (`git ls-files benchmarks/tier3` in ODE); its README.md,
  HANDOFF.md, `heldout/`, and `review/` are this repo's only.
- **Contestants:**
  - Codex CLI: gpt-6.1-sol, gpt-6-astra, and gpt-6-luna. Luna replaced
    gpt-5.6-terra on 2026-10-15. Codex was 0.162.1 through astra-high and is
    0.163.0 from astra-xhigh on; the SSE upgraded it and chose not to revert.
  - No Claude model competes (the SSE's decision, 2026-10-16, given usage plans):
    Opus and Fable design, orchestrate, and judge. The run root's sonnet-*, opus-*,
    and fable-* directories and prompts are unused.
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
- **Launch:** every run is an interactive session, on subscriptions; never the API
  or `codex exec`. The SSE launched Sol and Astra by hand. From luna-high on (the
  SSE's decision, 2026-10-16), Claude launches with
  `admin/launch_codex_run.sh <run> <model> <effort>`: a detached tmux session,
  watched with `tmux -L turbomicm attach -t <run>`. The steps it automates:
  - Start Codex with
    `cd <run> && codex -m <model> -c model_reasoning_effort=<effort> -c features.memories=false`.
  - Then give the session `/goal read PROMPT.md and work until done, commit`.
  - Always pass the effort: Codex's model picker saves the last choice in config.toml.
  - Tier 3: `admin/launch_codex_run.sh --tier3 <run> <model> <effort>` adds
    `--search`, names the tmux session `tier3-<run>`, and requires an untouched
    clone at the Tier 1 commit the prompt names.
- **Monitor:** follow the session record read-only, in the background:
  `until admin/run_state.py <rollout.jsonl> <codex pid>; do sleep 60; done`.
  - A run is done when its last event is `task_complete` and the record has been
    idle for 10 minutes; the goal starts each next turn within a second.
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
  - Run it on an idle machine on AC power, with no agent running. Check
    `pmset -g batt` first: after a reboot the laptop may be on battery; a
    launcher can wait for the adapter, as `tier3-work/run-timing-fable-on-ac.sh`
    did on 2026-10-19.
- **Inputs:** `reference/make_cells_csv.py 10000 <seed>` writes the cells file;
  `reference/compare_final.py YOURS MICM` applies the gate.

## Conventions

- **Results:** the wiki shows verified numbers only, with tables updated in place.
- **Wiki layout:** each tier has Design, Prompt, and Results pages (`Tier-1-*`,
  `Tier-2-*`, `Tier-3-*`), plus Home, Running-the-Experiment, and
  Role-of-Claude-in-the-Design.
  - Inside a table, use `[text](Page)` links: the `|` in `[[text|Page]]` splits
    the cell.
  Agents' self-reports appear only as a one-line finding.
- **Naming:** call the user "the SSE" (scientific software engineer). Name them once,
  in the preamble of Role-of-Claude-in-the-Design, and use no pronouns for them.
- **Judging:** since 2026-10-16 no Claude model competes, and Opus and Fable judge
  the entries. Keep the blind review blind: fresh sessions, anonymized copies.
  The Role page describes the roles; its old conflict-of-interest disclosure is
  gone.
- **Git:**
  - Commits end with the Co-Authored-By trailer.
  - Push this repo and the wiki; the SSE wants them pushed for backup.
  - Push the run-root branch too, after its files change.

## Pitfalls

- **Placement of agent instructions:** never put CLAUDE.md or AGENTS.md in the run
  root or any folder above it (`~/NSF-NCAR/ODE`, `~/NSF-NCAR`, `~`), because the
  harnesses load them into contestants' sessions.
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

## Tier 2: specialization for a frozen mechanism (adopted 2026-10-16; prompt not yet written)

How fast can a solver for TS1 alone get? Reference mechanisms like TS1 stay in
production unchanged for years; prototyping still needs run-time configurability.
The wiki pages are Tier-2-Design, Tier-2-Prompt (not yet written), and
Tier-2-Results.

- Start from each agent's own Tier 1 solution; the generic path and the Tier 1
  gate stay.
- For TS1 anything goes: literal indices and constants, a fully unrolled LU,
  fused kernels; same Rosenbrock method and tolerances; portability rules
  unchanged.
- The specialized C must come from a generator the agent keeps, so a revised
  mechanism costs a regeneration, not a new run. The generator must handle any
  mechanism (the SSE's decision, 2026-10-16), and the tests run it on others,
  such as MICM's tutorial mechanisms, so the result carries over to TS1's
  successors, such as TS2 and T4S.
- Verify by differential testing against the entry's generic path and MICM at
  the same settings: held-out seeds and day–night runs (with a tolerance that
  runs there).
- Measure speed against the entry's generic and generated Tier 1 paths and
  MICM, plus generation and compile time and tokens. Optionally a variant with
  rate parameters as data, for ensembles.
- Arms to be chosen, for example sol-max and astra-xhigh, which already
  generate kernels, and one entry that does not.
- Prompt: `tier2/prompt.v1`, its choices settled by the SSE on 2026-10-16: two
  benchmarks (the specialized one ranked, the Tier 1 one timed alongside); run-time
  data versus code as drafted; the linear algebra free (any exact factorization or
  direct solve); a generator for any mechanism; no build-time cap, unless builds
  start to grow without bound (build times are reported); the day-night gate; no
  web search; a blind review of the generator. Its
  day-night gate, agreement with MICM to 10 digits overall and 9 per species at
  rtol 1e-6, atol 1e-15, 32-cell groups, is calibrated in
  `results/2026-10-16-tier2-day-night-calibration` (MICM agrees with itself
  across compilers to 14 digits).
- Still to do: settle the draft's choices, choose the arms, give Tier 2 its own
  `scenarios/` (copies of the disclosed day-night scenarios, plus MICM's results
  at the gate's settings; not Tier 3's folder, which leaves tolerances to its
  agents), port `make_prompts.py` and `setup_runs.py` from `tier3/`, and give
  `admin/launch_codex_run.sh` a Tier 2 mode like its `--tier3`.

## Tier 3: algorithm design (approved 2026-10-15 as Tier 2; renumbered 2026-10-16; launched 2026-10-18)

- **Files:** `tier3/prompt.v1` is the template and `tier3/README.md` is the plan
  and build list; the wiki pages are Tier-3-Design, Tier-3-Prompt, and
  Tier-3-Results.
  - Renamed from `tier2/` on 2026-10-16, here and in the run root (ODE `a219270`):
    the template says Tier 3, and the unstarted sol-max and sol-xhigh runs in
    `ODE/benchmarks/tier3/` have regenerated prompts. Tier 3 runs will go on
    branches `tier3/<model>/<effort>`.
- **Order:** the SSE started Tier 3 on 2026-10-18, ahead of Luna's Tier 1 reruns
  and Tier 2: one run at a time, in the order of the Tier 1 runs (sol-xhigh,
  sol-max, astra-xhigh, astra-max), each collected, archived, and gate-checked
  before the next launches. All four ended and were verified that day (details
  in `admin/runs.md`).
- **Finished on the same laptop** (the SSE's decision of 2026-10-18 19:32,
  reversing a planned handoff): MICM's baselines, the blind review
  (`results/2026-10-18-tier3-blind-review/`), the timing session (21:19 to
  05:46), and the ranking (`results/2026-10-19-tier3-ranking/`; wiki:
  Tier-3-Results). The rules are in `tier3/README.md`, "How the ranking was
  made".
- **Arms:**
  - Sol and Astra at xhigh and max, a 2 x 2 matrix, each started from a clone of
    the run's own Tier 1 repository (`<tier1_commit>` placeholder);
  - Luna: one run at max, if its Tier 1 passes comfortably;
- **Rules:**
  - any mechanism-generic method; no mechanism reduction or learned surrogates;
  - the Tier 1 gate still passes;
  - web search for the literature (`codex --search`), no copied code, every source
    in REFERENCES.md;
  - Python prototypes measure accuracy and operation counts, while speed comes only
    from C, with a test checking each C method against its prototype;
  - scans logged in scans/ and SCANS.md;
  - the agent justifies its own accuracy target in ACCURACY.md.
- **Yardstick:** significant digits (Sandu et al., 1997) against converged,
  tight-tolerance MICM references, over species at or above 1e-14 mol m-3.
  - The gate: no NaN or infinity, and no value below -1e-20 mol m-3.
  - Rank by the fastest declared operating point that reaches L digits overall
    and L - 1 digits per species on every scenario, at L = 2 and L = 3.
  - Scenarios: f(t) photolysis applied per step, 24-hour day and night, 5–30
    minute steps, clean and polluted cells, some held out.
- **Before the ranking:** the tooling, the six disclosed scenarios, and their
  references (converged to 10.2–10.8 digits) are done; see `tier3/README.md`.
  - MICM's default fails on 3 of the 5 day–night scenarios (absolute tolerance
    1e-3 mol m-3); 1e-15 fixes it.
  - Undisclosed scenarios: five, approved by the SSE on 2026-10-18 and fixed in
    `make_heldout.sh` (SHA-256 `08016cac…195c` in `admin/runs.md`), kept
    outside the run root until the last run ended, then published in
    `tier3/heldout/` with the driver that makes their files (1,000 and 10,000
    cells), references, and the 10,000-cell disclosed set, in
    `~/NSF-NCAR/.tier3-heldout/`.
  - MICM's baselines: 32-cell groups, relative 1e-2 with absolute 1e-15 for
    level 2 and 3e-4 with 1e-15 for level 3, the loosest settings of an
    18-setting scan that reach the level on all 11 scenarios.
- **Fable's reference solution (2026-10-18/19, not a contestant):** built in
  `ODE/benchmarks/tier3/fable-max` from the Tier 1 `fable/` entry, collected
  as `tier3/fable/max`. Its levers: the step carried across intervals and
  quantized to them, the H211b controller (b = 4), an absolute tolerance of
  1e-20 (isoprene and phenol sit near the score's 1e-14 floor), RODAS3 for
  two digits, ROS3 under H211b for the cheapest three, extrapolated linearly
  implicit Euler with the shifted substep sequence for three with margin and
  more. Five declared points with 0.2 digit of margin per species on the
  disclosed scenarios (OPERATING_POINTS.md in the entry); 21 times MICM at
  level 2 and 8.4 to 12.4 times at level 3 on the six disclosed scenarios at
  1,000 cells in its own session. Beside the contest's winners under the
  ranking's protocol (2026-10-19, `results/2026-10-19-tier3-fable-vs-contest/`):
  L2 24.6 s over the 11 scenarios against astra-xhigh/fast 34.1 s, and L3r
  49.1 s against astra-xhigh/accurate 73.5 s (1.4 and 1.5 times faster, 18.9
  and 11.2 times MICM); the contest's ranking stands.
- **Ranking (2026-10-19):** astra-xhigh first at both levels: fast (RODAS3,
  relative 1e-2) at 34.7 s summed over the 11 scenarios, 13.4 times MICM;
  accurate (extrapolated linearly implicit Euler, order 4, 4e-3) at 74.5 s,
  7.3 times MICM. Then sol-max, sol-xhigh, astra-max at level 2 and
  sol-xhigh, sol-max, astra-max at level 3 (the last two 2.6 % apart). The
  held-out scenarios changed no rank; astra-max's fast3 misses level 3 on two
  of them. Every entry's 10,000-cell scores matched its 1,000-cell
  verification within 0.01 digit.

## Open items

- Rank with an undisclosed seed, and run the blind review.
- Time the entries on an x86-64 Linux machine.
- Tier 3: Luna's arm, if its Tier 1 passes (it needs new held-out scenarios,
  since `tier3/heldout/` is public). Fable's reference solution (not a
  contestant) is done and collected as `tier3/fable/max`, verified in
  `results/2026-10-19-tier3-fable-max-verification/` and timed beside the
  contest's winners in `results/2026-10-19-tier3-fable-vs-contest/`.
- Verification worktrees of `runs/*`, MICM's timed programs, and the seed-2026
  gate files are in `~/NSF-NCAR/tier1-work` (outside the run root), with its
  scripts: `prepare.sh` rebuilds them, `time-when-idle.sh LABEL RUN ...` waits
  for an idle machine (no Codex record written for 10 minutes, CPU below one
  core, AC) and then runs `run-timing.sh`, and `progress.py` summarizes a live
  run. If it is gone, `git worktree prune` and rebuild per "Verify an entry".
- Time each new entry with all the others in one session (`admin/time_entries.py`,
  5 rounds at 10 and 1 threads), alongside `micm-standard-in-place-mozart`, its
  `-g32` variant, and `micm-vector8-in-place-mozart-g32`, as in
  `results/2026-10-19-luna-max-verification`.
  - `admin/leaderboard.py` regenerates the wiki's Leaderboard: at 10 threads from
    every session's `timings.csv` except the Sol verification's, whose MICM names
    it does not parse; at 1 thread from the sessions that ran it.
  - Add each entry's grouping to `ENTRY_GROUPING` in `admin/leaderboard.py`.
  - When building in parallel, make the `micm/` and `reference/` links first, or
    with `ln -sn`: a racing `ln -s` writes a link inside the run root.
- Luna: run all three fresh, one after another, at max, xhigh, then high (the
  SSE's order, 2026-10-19), each verified and timed with all entries before the
  next starts. luna-max is done: finished 21:05, verified, and timed overnight
  (`results/2026-10-19-luna-max-verification`), correct but 0.90 times MICM's
  default. luna-xhigh is next. The first luna-high attempt was stopped at 08:11
  on 2026-10-16, after 1 h 23 min, far from done (615 failing cases at its 07:51
  checkpoint); it is kept, outside the run root, at
  `~/NSF-NCAR/aborted-runs/luna-high-2026-10-16`, and `luna-high/` was recreated.
  For each run: `admin/launch_codex_run.sh <run> gpt-6-luna <effort>`, then arm
  `admin/run_state.py` (a monitor does not outlive its Claude Code session).
- Next: an analysis pass on every entry. The wiki's Tier-1-Implementations page
  describes each entry; this pass measures, and, now that Claude judges, it may
  evaluate as well.
  - Instruction counts and cycles from `/usr/bin/time -l`, at 1 thread, differenced
    between two step counts to cancel setup.
  - Time by solver phase (`xcrun xctrace`, Time Profiler or CPU Counters), and code
    size.
  - Generated against generic kernels within one entry: astra-xhigh's benchmark
    takes `--generic`; sol-max's loader takes `chem_load_ts1(0)` or `(1)`
    (`src/bench.c:103`), so time a scratch copy built with 0.
  - astra-high's rate constants: their share of its timed loop, which recomputes
    them every step (`src/ts1_bench.c:154`).
  - Compile time of the generated kernels (16,000 to 21,000 lines of C), with
    clang and gcc-15.
- Framework, after the contest: Opus and Fable build on its results a framework
  usable from MUSICA, CheMPAS-A, CAM-chem, and other host models (wiki:
  From-Contest-to-Framework). First, a MUSICA integration study, after the
  analysis pass: how mechanism-specific kernels could reach MUSICA's API.
  - Ahead of time: a CMake function that generates and links a mechanism's
    kernels when the host model is built, as KPP-based models do.
  - At run time: generate C, compile it once with the system compiler, cache the
    shared library by a hash of mechanism, compiler, and flags, and `dlopen` it.
  - Embedded JIT (LLVM ORC or MLIR): MICM removed its LLVM JIT in PR #886,
    2025-12-10, as unused, not shipped with musica, costly to keep in step with
    interface changes, and pinned to an old LLVM.
  - In any design: keep the generic kernels as the fallback, check a fingerprint
    of the index plan when loading generated kernels, register kernels through
    function pointers, and keep one implementation of the species ordering.
  - Write it up on the wiki once measured.
- Current status is in `admin/runs.md` and on the wiki's Tier-1-Results and
  Tier-2-Results pages.
