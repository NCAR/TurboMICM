# Blind review of nine Tier 2 entries

You are one of two reviewers in the blind review of a controlled experiment that
compares LLM coding agents. Nine agents received the same task, PROMPT.md in this
folder: starting from their own earlier port of MICM's CPU chemistry solver to
portable C, write a code generator, in Python, that reads a mechanism
configuration and emits C specialized to that mechanism, keeping MICM's
numerical method, and use it to make the TS1 mechanism as fast as possible.
Their final repositories are in `A/` to `I/`, in an order and with labels chosen
at random. Each has had its git history, local paths, and run names removed;
`ROOT/entry` in a file stands for the entry's own directory. `scenarios/` and
`reference/` hold the README files the agents were given, and `mechanisms/`
holds two mechanism configurations none of them used, Chapman and Carbon Bond 5,
in the same configuration format.

Correctness and speed are judged separately, by measurement: every entry passed
the same gate (MICM's test cases, agreement with MICM on TS1 and on ten day-night
scenarios, and the portability rules). This review is the secondary ranking, of
the generator, as PROMPT.md defines it:

1. *Readability:* can a scientific software engineer follow the generator: how
   it reads a mechanism, what it computes symbolically (ordering, sparsity, the
   LU steps), what it emits, and how the generated C maps back to the
   chemistry? Structure, naming, comments, and documentation (README.md,
   DESIGN.md, SPECIALIZATION.md).
2. *Maintainability:* tests of the generator (what they compare, on which
   mechanisms, and how tightly), the build (does `make` regenerate whatever a
   change requires?), separation of the generator from the generic library and
   the solver, error handling for bad input, and how well the documentation
   matches the code.
3. *Extensibility:* how hard the prompt's three exercises are, done as the
   entry's own documentation says they should be: (a) add a reaction to TS1 and
   regenerate; (b) change a rate parameter of an existing TS1 reaction; (c) run
   the generator on another mechanism, for example `mechanisms/chapman` or
   `mechanisms/carbon_bond_5`, and check its result against the generic
   library. Also: how hard it would be to add a new kind of rate constant, or
   a new specialization.

## How to work

- Read only inside this folder. Do not look anywhere else on this machine or
  online for these entries, their authors, or other solutions; do not try to
  work out which agent wrote which entry.
- Do not build or run the code: the machine may be timing programs. Read the
  files; `grep`, `wc`, and similar light commands are fine. Do the exercises
  on paper: for each one, name the files to edit and the commands to run, from
  the entry's own documentation and code, and say what you expect to break or
  to need care. Where the documentation is silent or wrong, say so.
- Judge each entry on what is in its folder. Verify claims where you can from
  its code (for example, that the generator emits what SPECIALIZATION.md says,
  or that a test really compares the generated solver with the generic one),
  and say when you could not.
- Review the entries in the order you are given, then compare them.

## Report

Write your report as Markdown, in this order:

1. For each entry, a section with: a summary of its generator, what it
   specializes, and how it is tested (three to five sentences); a score from 1
   to 10 for each of the three criteria, each with two to four sentences of
   justification citing files; the three exercises on paper, each as a short
   list of edits and commands with its difficulty (easy, moderate, hard, or not
   supported); and the entry's three main strengths and three main weaknesses.
2. A comparison: the nine entries ranked by the sum of their three scores, with
   the reasons for the order, and any ties.
3. At the end, one JSON block, exactly in this form:

```json
{"scores": {"A": [1, 2, 3], "B": [...], "C": [...], "D": [...], "E": [...],
            "F": [...], "G": [...], "H": [...], "I": [...]},
 "exercises": {"A": ["easy", "easy", "moderate"], "B": [...], ...},
 "rank": ["?", "?", "?", "?", "?", "?", "?", "?", "?"]}
```

with the three scores in the order of the criteria above, the exercises'
difficulties in the order (a), (b), (c), and the rank from best to worst.
