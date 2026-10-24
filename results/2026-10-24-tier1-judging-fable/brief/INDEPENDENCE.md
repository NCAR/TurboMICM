# Independence of the second judgment

The second judge, Claude Fable 5.1 at max effort, scores after the first judge
(Claude Opus 5.5) has reported. To keep the second judgment independent, every
helper of the second judge works from the code and MICM only, and does not read:

- `~/NSF-NCAR/TurboMICM/results/2026-10-20-tier1-judging-opus/` (the first
  judge's notes, exercise reports, and test traces);
- the wiki page `Tier-1-Judging-Opus.md`, and the "Results" section of
  `Tier-1-Judging.md`;
- the "Extension points" section of each wiki page `Entrant-<run>.md` (it
  summarizes the first judge's exercises). The rest of an entrant page, its
  Architecture table in particular, may serve as a map; check everything
  against the code.

Running anything in an entry's directory is not allowed. For E5 only, a helper
may copy an entry's Python tools and a configuration into its own scratch
folder and run the copy there with `python3 -I`; nothing in C is built or run.
