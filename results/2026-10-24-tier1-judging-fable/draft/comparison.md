## Compared with Judge 1

Judge 1's report page ([Tier-1-Judging-Opus](Tier-1-Judging-Opus)) and its
`scores.json` were read after the scores above were committed (TurboMICM
`1b1d404`, wiki `763bd8d`). Nothing above was changed afterwards. This section
records what the two judgments share, where their standards part, and what
Judge 1 found that this judgment had missed, checked in the code after the
fact.

### The two sets of scores

| Entrant | Judge 1 total (rank) | Judge 2 total (rank) | R 1 / 2 | M 1 / 2 | E1–E5 Judge 1 | E1–E5 Judge 2 |
|---|---:|---:|:-:|:-:|---|---|
| astra-high | 22.4 (3) | 24.8 (1) | 6 / 8 | 8 / 8 | 5, 4, 5, 3, 4 | 5, 4, 5, **4**, 4 |
| sol-xhigh | 22.6 (2) | 23.4 (2) | 8 / 7 | 7 / 8 | 5, 4, 4, 3, 3 | 5, 4, **5**, 3, **4** |
| astra-xhigh | 23.8 (1) | 23.2 (3) | 7 / 7 | 8 / 7 | 5, 5, 5, 4, 3 | 5, 5, 5, 4, **4** |
| sol-high | 20.6 (5) | 23.0 (4) | 6 / 7 | 7 / 8 | 5, 4, 3, 3, 4 | 5, 4, **4**, 3, 4 |
| astra-max | 22.4 (4) | 21.4 (5) | 7 / 6 | 7 / 7 | 5, 4, 5, 3, 4 | 5, 4, 5, 3, 4 |
| sol-max | 19.2 (6) | 20.0 (6) | 6 / 6 | 6 / 6 | 5, 4, 3, 3, 3 | 5, 4, **4**, 3, **4** |
| luna-max | 13.8 (8) | 16.2 (7) | 4 / 5 | 3 / 4 | 4, 4, 4, 3, 2 | 4, **3**, **5**, 3, **3** |
| luna-xhigh | 15.0 (7) | 15.4 (8) | 5 / 4 | 4 / 5 | 5, 3, 2, 3, 2 | 5, **2**, **3**, 3, **3** |
| luna-high | 9.4 (9) | 9.4 (9) | 3 / 3 | 2 / 2 | 2, 2, 2, 3, 2 | 2, 2, 2, 3, 2 |

Bold marks an exercise rating that differs from Judge 1's.

- **The same three groups.** The six Sol and Astra entries span 19.2 to 24.8
  in the two judgments together, the three Luna entries 9.4 to 16.2, and
  luna-high is last in both with identical component scores.
- **The same three entrants lead,** astra-xhigh, sol-xhigh and astra-high,
  in a different order; sol-xhigh is second and sol-max sixth in both;
  sol-high and astra-max exchange ranks 4 and 5, luna-max and luna-xhigh
  ranks 7 and 8. The Spearman rank correlation is 0.90.
- **This judgment's totals are higher by 0.8 on average** (19.6 against 18.8
  of 30); the largest differences are 2.4 (astra-high, sol-high, luna-max).
  Readability is within one point for eight entrants and identical for three;
  maintainability is within one point for all nine and identical for four; 32
  of the 45 exercise ratings are identical, the 13 others are one point
  apart, 11 of them with this judgment higher.

### Where the standards differ

- **E3, five entrants.** Judge 1 graded the context a replacement solver gets:
  per state or per solver 5, one context shared by every group 4, no context
  3, no hook 2. This judgment graded presence: a pointer pair with a context
  5, without one 4, an enum over built-in methods 3, nothing 2. Hence
  sol-xhigh and luna-max 5 against 4, sol-high and sol-max 4 against 3,
  luna-xhigh 3 against 2. The code cited is the same.
- **E5, five entrants.** Judge 1 charges the generated kernels' dependence on
  Python copies of the ordering to E5 (sol-max and astra-xhigh 3); this
  judgment charges it to maintainability (astra-xhigh 7 against 8; sol-max 6
  in both) and not again under E5. Judge 1 gives the three Luna readers 2 for
  their fixed unit exponents; this judgment gives 3 for wrong coefficients
  and 2 for wrong reactions (luna-high). Judge 1 also docks sol-xhigh for
  scaling user-defined rates, taken up below.
- **E2, two entrants.** luna-max 3 against 4: its one `callback_context`
  serves the controller, the factor and the solve for every group the solver
  runs in parallel, with no group identity, so a per-group PI state needs an
  API change. luna-xhigh 2 against 3: no hook, and the rules written in two
  integrators, which this judgment rated as it rated luna-high.
- **E4, one entrant.** astra-high 4 against 3; the callback scan below.
- **Readability.** The two-point difference on astra-high (8 against 6)
  decides rank 1. Judge 1 weighed the timed path's spread over about twenty
  functions in eight files, hooks that loop over empty lists, 43 lines over
  120 characters and a name reused for two types; this judgment weighed the
  rate union with MICM's field names, the 80-line integrator with MICM's
  names and structure, the highest comment share of the nine and the design
  document's file map. The one-point differences: sol-high 7 against 6 (the
  line-for-line integrator and its comments at the correspondences with MICM,
  against the positional `p[24]` both cite); sol-xhigh 7 against 8 (both cite
  the 297-line constructor; this judgment also counts what `process.c` holds
  unannounced and the 1.0 % comment share); astra-max 6 against 7 (the
  positional `p[32]` slots and the order table by enum position, against the
  extension guide Judge 1 credits); luna-max 5 against 4 and luna-xhigh 4
  against 5, adjacent in both.
- **Maintainability.** astra-xhigh 7 against 8 (the unchecked kernel
  indices, charged here); sol-xhigh 8 against 7 and sol-high 8 against 7 (one
  more faithful case each in this judgment's sample, sol-xhigh's reader tests
  and generated tutorials, sol-high's supplemental TS1 runs at four group
  sizes and its inspection of every candidate value); luna-max 4 against 3
  and luna-xhigh 5 against 4 (the per-id suites and the TS1 gate at one and
  four threads).

### The test sample

The same twelve cases. 81 of the 108 grades agree. Of the 27 that differ, 15
are F against P: this judgment's helpers graded faithful where Judge 1's
graded partial in the rejected-step case (sol-high, sol-xhigh, sol-max,
astra-max, luna-xhigh), the Oregonator (sol-max, astra-max, luna-max,
luna-xhigh), Robertson (sol-high, astra-max) and the regression Jacobian
(astra-high), and partial where Judge 1's graded faithful for the surface
store case in sol-high and the invalid-arguments state case in sol-max and
astra-max. The other 12 are S against P or F, all in the Luna entries, most
of them the rule stated above: an id whose verdict is the exit status of a
binary shared by many ids is S even where a matching check exists, so
luna-xhigh's and luna-high's Troe, surface-store and tutorial cases and
luna-high's state and Oregonator cases are S here and P or F there; in the
other direction, luna-xhigh's random LU and DAE cases and luna-high's
rejected-step case are P here and S there. Faithful of twelve, for sol-high,
sol-xhigh, sol-max, astra-high, astra-xhigh and astra-max: Judge 1 10, 8, 9,
11, 12, 7; this judgment 11, 9, 10, 12, 12, 9. For luna-max, luna-xhigh and
luna-high: Judge 1 1, 2, 4; this judgment 2, 3, 1. Both find the same single
M, the Chapman `factored_alpha_minus_jac` case in luna-high, and both read
every luna-high line as printed from a list.

### What Judge 1 found that this judgment had missed

Checked after its page was read, by reading the code and, for the second
item, compiling two declarations in a scratch folder; nothing in the entries
was built.

1. **luna-xhigh's three-stage `e[1]`** is −2.907955871680469821793… where
   MICM's is −2.9079558716805469821718… (`src/chem.c:197`; MICM
   `include/micm/solver/rosenbrock_solver_parameters.hpp:227`): 2.7e-14
   relative, confirmed. The gate's 1e-10 cannot see it.
2. **luna-high with a mechanism that has no custom parameters.** The
   generator then emits `custom_labels[] = {}` and
   `default_custom_parameters[] = {}` (`tools/generate_c.py:215-219`); with
   the entry's own flags, `-std=c11 -Wall -Wextra -Wpedantic`, gcc-15 rejects
   the file ("zero or negative size array") and clang compiles it with two
   warnings; and the benchmark divides by the custom-parameter count
   (`src/main.c:242`). Confirmed. Its warm-up loop reads `status` by cell and
   writes it by worker (`src/main.c:312-315`), a race; confirmed by reading.
3. **sol-xhigh scales user-defined rates by the unit conversion.** For a
   configuration in molecule cm-3 units, the default (`tools/configure.py:79-82`),
   the reader multiplies a user-defined reaction's scaling factor by
   6.022e17 for each reactant beyond the first (`:149-150`); MICM's
   user-defined rate constant is the per-cell value times the scaling factor
   as given (`include/micm/process/rate_constant/user_defined_rate_constant.hpp:17,26-28`).
   DESIGN.md states the rule (`DESIGN.md:121-126`), and the tutorials'
   user-defined reactions are first order, so the tests cannot see it; a
   MICM-format mechanism with a second-order user-defined reaction, run with
   values in MICM's units, would differ from MICM by that factor. This
   judgment's probe had no such reaction. Under the E5 standard above it
   would have cost a point: E5 3, extensibility 8.0, a total of 23.0, a tie
   with sol-high behind astra-xhigh. The scores stand as committed; the point
   is recorded for Judge 3.
4. **astra-high's `chem_set_rate_callback`** scans from the first callback to
   the end of the type-sorted process array (`src/mechanism.c:146`), so a
   rate type added after `CHEM_CALLBACK` in the enum would be scanned too,
   its label compared by `strcmp`, whenever a name matches no callback; the
   recipe (`DESIGN.md:155-158`) does not say to insert the new value before
   `CHEM_CALLBACK`. Confirmed. The E4 rating of 4 stands under the standard
   as applied, a coupling that bites an appended type's own evaluation, of
   which there is none here; the hazard is recorded.
5. **astra-max's `RATE_NAMES`** has ten names against eleven enum values with
   `CH_LAMBDA` last (`tools/configure.py:17-19`; `include/chem.h:28-31`), so a
   name appended to the list takes `CH_LAMBDA`'s value in C. Confirmed; it is
   the coupling the E4 table records for this entry, rated 3 by both judges.

Findings of this judgment that Judge 1's page does not record: luna-max's
four Troe ids evaluating one fixed case without MICM's parameters, and its
`DAESolveEnforcesConservation` id falling to an unrelated constraint step;
luna-xhigh's CAMP reader emitting every initial concentration as zero; the
E5 probe's result that luna-high emits all ten reactions wrong on the
species-reversed tutorial while the six Sol and Astra generators survive the
reversal; sol-high's inspection of every candidate value, which makes its
floored norm harmless; and the three hook hazards in the E3 table
(astra-xhigh's state copy, astra-max's `ch_state_set_lu`, astra-high's
half-NULL pair).

### For Judge 3

The two judgments agree on the groups, on luna-high, on sol-xhigh second and
sol-max sixth. They differ on which of astra-high and astra-xhigh leads,
which turns on astra-high's readability and on whether astra-xhigh's
unchecked kernel indices count once or twice; on sol-high against astra-max
at ranks 4 and 5; and on luna-max against luna-xhigh at 7 and 8, where
luna-max has the higher extensibility in both and the order turns on whether
its callbacks or luna-xhigh's per-id tests weigh more. Items 3 and 4 above
were not in this judgment's evidence; their effect on the scores is given
there.
