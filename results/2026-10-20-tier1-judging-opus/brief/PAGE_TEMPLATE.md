# Entrant page: template and rules

You are writing one page of the TurboMICM wiki (GitHub wiki, Markdown with
Mermaid): the page for one Tier 1 entrant. Nine agents write the nine pages at
the same time from this template, so the pages must match in format and in
content: the same sections, in the same order, with the same table rows and the
same six diagrams. Where something does not apply to your entrant, keep the
section or row and say so in a few words ("none", "not used").

## The contest, in short

TurboMICM Tier 1: each entrant (one OpenAI model at one reasoning effort, in
Codex CLI) received the same prompt: port MICM's CPU chemistry solver (C++) to
portable ISO C11 with libm and OpenMP, keep MICM's numerical methods, port
MICM's 680 test cases, read MICM's YAML/JSON configurations with a Python
reader, translate the TS1 mechanism (210 species, 547 reactions) from MICM's
`ts1_mechanism.hpp`, and make the TS1 benchmark fast (10,000 cells, 100 steps
of 30 s). Read `~/NSF-NCAR/TurboMICM.wiki/Tier-1-Design.md` first
for the rules, and `~/NSF-NCAR/TurboMICM.wiki/Tier-1-Prompt.md` if
you need the prompt's exact wording.

## Sources

1. The entry at its final commit, read-only:
   `~/NSF-NCAR/tier1-work/verify/<run>/clang/<run>/` (a git
   worktree; `git log`, `git show` work there; `bin/` and generated files are
   from a verification build). Read README.md, DESIGN.md, DECISIONS.md,
   RESULTS.md, the Makefile, and the C and Python sources.
2. Its section of `~/NSF-NCAR/TurboMICM.wiki/Tier-1-Implementations.md`
   (heading `### <run>`): written earlier to a common template, with about 20
   of its citations checked against the code. Reuse its content freely;
   re-check every citation you carry over.
3. Its section of `~/NSF-NCAR/TurboMICM/admin/runs.md` (the run
   log), and the notes about it in `~/NSF-NCAR/TurboMICM.wiki/Tier-1-Results.md`.
4. Its session, condensed: `<scratch>/judging/digests/<run>.txt`
   has one line per agent message (SAY), shell command (RUN), other tool call
   (TOOL), and output line that reports a time or a test count (OUT), with
   times in MDT. Use it for "How it got there": the order in which the design
   and its optimizations came, what the agent tried and kept or dropped, and
   why it said it did so.
5. MICM itself, for comparisons: `~/NSF-NCAR/ODE/benchmarks/micm/`
   (for example `include/micm/solver/rosenbrock.inl`,
   `include/micm/solver/lu_decomposition_mozart_in_place.hpp`).

Work read-only: do not edit, build, or run anything in the entry's directory or
anywhere outside your scratch output folder. Light commands (`grep`, `wc`,
`git log`, `git show`, `sed -n`, `python3` reading files) are fine. Do not run
the benchmark or the tests.

## Rules for the content

- **Describe; do not judge.** The judges' scores are on another page. Do not
  write evaluative words (clean, elegant, poor, sloppy, impressive,
  unfortunately, well-designed). State what the code does, with evidence.
- **Cite the code.** Every statement about the code carries `path:line` or
  `path:start-end` at the final commit, in backticks, as
  Tier-1-Implementations does. Check each one.
- **Verified numbers only**, from the facts you are given. The agent's own
  measurements may appear in one line only, in "How it got there", introduced
  as its own ("By its own measurements, ..." or "Its RESULTS.md reports ...").
- **Style of this wiki:** plain declarative sentences; American spelling; the
  person who runs the experiment is "the SSE", with no pronouns; times in MDT;
  thousands separators (7,097 ms); `×` for speedups (3.97×); en dashes for
  ranges (5–10 %). Name entries by run name (`sol-max`), and MICM's
  configurations as Tier-1-Results does ("MICM's default": standard ordering,
  in-place Mozart LU, one group per thread; "MICM's fastest configuration":
  vector8 ordering, in-place Mozart LU, 32-cell groups).
- **Links:** wiki links as `[[Text|Page-Name]]` outside tables and
  `[Text](Page-Name)` inside tables. Code links to GitHub as
  `https://github.com/NCAR/TurboMICM/blob/runs/<model>/<effort>/<path>#L<n>`
  only in the Records section; elsewhere plain `path:line` in backticks.
- Do not name or describe Claude's own reference entry or the other entrants,
  except for a comparison that is a plain fact and needed (for example "the
  same layout as the six Sol and Astra entries").
- Keep the page between 220 and 380 lines of Markdown.

## Diagrams: the same six on every page

Mermaid only, and only `flowchart` (TB or LR), `sequenceDiagram` and `pie`.
GitHub renders Mermaid in wiki pages; keep to plain syntax.

Colors (the wiki's key: blue for the agents' work, grey for MICM and reference
material, green for checks, gold for results, red for problems). In every
flowchart, end with exactly these lines, keeping only the classes you use:

    classDef agent fill:#ddf4ff,stroke:#54aeff,color:#0a3069
    classDef ref fill:#eaeef2,stroke:#8c959f,color:#1f2328
    classDef check fill:#dafbe1,stroke:#4ac26b,color:#0f5323
    classDef result fill:#fff8c5,stroke:#d4a72c,color:#4d2d00
    classDef warn fill:#ffebe9,stroke:#ff8182,color:#82071e

and style every subgraph with
`style <id> fill:none,stroke:#8c959f,stroke-dasharray:4 3`.

Shapes: code (a file or module) is a rectangle `id["label"]`; data and
generated files are cylinders `id[("label")]`. Labels always in double
quotes; use `<br/>` for line breaks; no `#`, `;`, `&`, `{`, `}`, or bare `"`
inside labels; never use `end`, `graph`, or `style` as a node id. Keep labels
short (at most about 40 characters per line).

1. **Architecture** (`flowchart TB`): three subgraphs, `inputs["Inputs from
   MICM"]`, `build["Build time: Python"]`, `run["Run time: C"]`, plus the
   benchmark, its output (class `result`), and the tests (class `check`).
   Nodes are the real files and functions on the TS1 path, with arrows for the
   data flow, as in this shape:

   ```mermaid
   flowchart TB
     subgraph inputs["Inputs from MICM"]
       direction LR
       hpp[("test/benchmark/ts1_mechanism.hpp")]:::ref
       cfg[("tutorial configurations")]:::ref
     end
     subgraph build["Build time: Python"]
       direction LR
       tr["tools/translate_ts1.py"]:::agent --> js[("configs/ts1.json")]:::agent
       js --> gen["tools/configure.py"]:::agent --> gc[("generated/ts1.c")]:::agent
     end
     subgraph run["Run time: C"]
       direction LR
       model["src/model.c<br/>model, LU plan"]:::agent --> solve["src/rosenbrock.c<br/>chem_solve"]:::agent
     end
     hpp --> tr
     cfg --> gen
     gc --> model
     solve --> bench["src/bench.c<br/>bin/ts1_bench_cpu"]:::agent --> out[("ts1_final.csv, elapsed_ms")]:::result
     solve --> tests["tests: make test-cpu"]:::check
     style inputs fill:none,stroke:#8c959f,stroke-dasharray:4 3
     style build fill:none,stroke:#8c959f,stroke-dasharray:4 3
     style run fill:none,stroke:#8c959f,stroke-dasharray:4 3
     classDef agent fill:#ddf4ff,stroke:#54aeff,color:#0a3069
     classDef ref fill:#eaeef2,stroke:#8c959f,color:#1f2328
     classDef check fill:#dafbe1,stroke:#4ac26b,color:#0f5323
     classDef result fill:#fff8c5,stroke:#d4a72c,color:#4d2d00
   ```

2. **Data layout** (`flowchart LR`): the arrays one solver call works on for
   one cell group, each with its index form and size for TS1 (for example
   `y[s·G + c]`, 210 × G), any copy between layouts, and a `check` node naming
   the loop that is innermost and whether it carries `omp simd`.
3. **The timed loop** (`sequenceDiagram`, with `autonumber`): participants for
   the benchmark driver, the OpenMP threads, the solver call, and the
   kernels/LU. Show what runs before the timer (setup, rate constants,
   warm-up), the timed region in `rect rgb(255, 248, 197)` with a `Note over`
   naming it, the loops over 30 s steps and internal steps with the work in
   each, accept/reject, the clamp, and what runs after the timer.
4. **Sparse linear algebra** (`flowchart LR`): two subgraphs, `once["Once,
   when ... is built"]` (ordering, fill, plan or lists) and `each["Each
   attempted step"]` (the diagonal shift, the factorization, the stage solves),
   with the counts you can establish (nonzeros after fill, update operations)
   and the alternatives the code offers (for example Doolittle).
5. **How it got there** (`flowchart LR`): five to nine milestones from the git
   history and the digest, left to right: `"HH:MM<br/>label"` or
   `"HH:MM · abc1234<br/>label"`, class `ref` for the start, `agent` for
   design and optimization steps, `check` for the TS1 gate first passing and the
   test suite complete, `warn` for a setback (if any), `result` for the final
   commit.
6. **Code composition**: leave the marker `<!-- CODE COMPOSITION -->`; the
   judge adds this pie chart from line counts made the same way for all nine.

Check every diagram renders before you finish:

    python3 <scratch>/judging/checkmmd.py <your page>

## The page

Write exactly these sections, in this order. Text in angle brackets is for
you to fill; everything else is literal.

```
# <run>

*This page names the model behind the entry. A blind reviewer should not read it.*

**<Model display name> at <effort word> effort.** One of the nine entrants in
Tier 1 of TurboMICM, which asks each agent to port MICM's CPU chemistry solver
to portable C and make it fast while keeping MICM's numerical methods.
[[Tier 1 Design|Tier-1-Design]] has the rules, [[Tier 1 Results|Tier-1-Results]]
the verified results, and [[Tier 1 Judging|Tier-1-Judging]] the judges' review
of the secondary criteria.

<Summary: one paragraph of four to six sentences. What the entry is (library,
integrators, how a mechanism enters), its defining design choices for speed,
where it stands on the verified numbers, and the one or two facts a reader
most needs (for example something unusual in what it times or tests).>

## At a glance

| | |
|---|---|
| **Model and effort** | <from the facts> |
| **Harness** | <from the facts> |
| **Session** | <from the facts> |
| **Tokens** | <from the facts> |
| **Final commit** | <from the facts> |
| **Tests, clang and gcc-15** | <from the facts> |
| **TS1 gate** | <from the facts> |
| **TS1, 10 threads** | <from the facts> |
| **TS1, 1 thread** | <from the facts> |
| **Speedup, 1 to 10 threads** | <from the facts> |
| **C library and benchmark** | <from the facts> |
| **Tests** | <from the facts> |
| **Python tooling** | <from the facts> |
| **Generated for TS1** | <from the facts> |
| **Cell groups and threads** | <short phrase, cited> |
| **Data layout** | <short phrase, cited> |
| **Forcing and Jacobian** | <short phrase, cited> |
| **Sparse LU** | <short phrase, cited> |

## Architecture

<Diagram 1>

<One sentence on how a mechanism reaches the solver.>

| Part | Files | What it does |
|---|---|---|
| Configuration | | |
| TS1 translation | | |
| Configuration reader and generator | | |
| Generated code | | |
| Mechanism or model | | |
| Solver state | | |
| Rate constants | | |
| Forcing and Jacobian | | |
| Sparse LU and solves | | |
| Rosenbrock step and step control | | |
| Other integrators and APIs | | |
| Benchmark driver | | |
| Tests and case runner | | |

## Data layout

<Diagram 2>

<Three to six bullets, cited: the index form of each array, padding and
alignment, where the group size comes from, copies between layouts.>

## The timed loop

<Diagram 3>

<Three to six bullets, cited: what is inside and outside the timer, how
threads take groups, what each solver call does, and how it compares with
what MICM's driver times (setup and rate constants before the warm-up, one
warm-up step, then the timed steps).>

## Sparse linear algebra

<Diagram 4>

<Three to six bullets, cited: ordering, fill, the plan or lists, the
factorization and solves, the storage, the alternatives.>

## How it optimized the code

### What it does for speed

| Technique | In this entry | Where |
|---|---|---|
| Layout with the cells contiguous | | |
| Cell groups | | |
| Threading | | |
| SIMD: `omp simd` and `restrict` | | |
| Precomputed LU plan | | |
| Species reordering | | |
| Mechanism-specific generated code | | |
| No allocation on the solve path | | |
| Work kept out of the step loop | | |
| Reuse within a step | | |
| Other | | |

### How it got there

<Diagram 5>

<Two to five paragraphs, from the git history (commit hash and time) and the
digest: when the solver took shape, when the gate first passed, when and in
what order each optimization came, what the agent measured and tried and
dropped, and how the rest of the session went (test ports, documentation).
The agent's own numbers in one line at most, introduced as its own.>

## Correctness and tests

<Four to seven bullets, cited: in what language and how MICM's cases were
ported (number of files, programs, or modules), how `make test-cpu` produces
the `<id> PASS|FAIL|N/A` lines (does each id come from a check that runs for
that id?), which cases are N/A and why, what checks TS1 or the benchmark, and
the verified result.>

## Relation to MICM

### From C++ templates to C

- **Matrix layout and vector width.** <cited>
- **LU storage.** <cited>
- **LU algorithm.** <cited>
- **Solver parameters.** <cited>

### Where it departs from MICM

<Bullets, cited: any difference from MICM's methods or from what MICM's
driver times, each with MICM's own line for comparison. If you find none,
say "None found in the TS1 path:" and list what you compared.>

## Extension points

<!-- EXTENSION POINTS -->

## Code composition

<!-- CODE COMPOSITION -->

## Beyond Tier 1

<!-- BEYOND TIER 1 -->

## Records

- **Repository:** branch [`runs/<model>/<effort>`](https://github.com/NCAR/TurboMICM/tree/runs/<model>/<effort>), final commit `<hash>`.
- **Verification and timing:** <results folders, as links to https://github.com/NCAR/TurboMICM/tree/main/results/...>
- **Run log:** [`admin/runs.md`](https://github.com/NCAR/TurboMICM/blob/main/admin/runs.md), section "<run>".
- **Transcript:** `transcripts/<run>.jsonl.gz` on the [`transcripts`](https://github.com/NCAR/TurboMICM/tree/transcripts) branch.
- **Comparison with the other entries:** [[Tier 1 Implementations|Tier-1-Implementations]] and the [[Leaderboard]].
```

## What to hand back

1. The page: `<scratch>/judging/pages/Entrant-<run>.md`,
   with diagrams 1 to 5 passing `checkmmd.py`.
2. Notes for the judge, not for the wiki:
   `<scratch>/judging/notes/<run>.md`:
   facts you met that bear on readability and maintainability, each cited:
   how easy the TS1 path is to trace and how many functions it passes
   through; naming; large or duplicated functions; logic kept in two places
   (for example Python and C both computing the species order); dead or
   unreachable code; error handling (return codes, allocation failures,
   invalid input); places where DESIGN.md or README.md disagree with the
   code; build hygiene; anything that looks like a defect (say how sure you
   are, and how you checked). Facts only, no scores, at most 80 lines.
3. Your final message: at most 150 words: the two paths, the six diagrams'
   check status, any citation in Tier-1-Implementations that you found wrong,
   and any defect you found.
