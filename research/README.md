# Research notes

The literature was searched before the methods were written, in six focused
reviews run in parallel by research assistants on 2026-10-18, each given a
written brief (the questions at the head of each report), web search, and
the instruction to read primary sources, to say for every source whether it
was read in full, in part, or seen only as an abstract or a citation, and
never to invent a citation. The six reports are kept here as written, with
their own open questions; the author then checked the sources that the design
rests on (the method tableaus, the controller formulas and constants, the
accuracy conventions and tolerances) and recorded in `../REFERENCES.md` what
was taken from each, with the reading status after the check. Where a report
and the author's check disagree, REFERENCES.md is authoritative.

| Report | Question |
|---|---|
| `R1-rosenbrock-w-methods.md` | Rosenbrock and W-methods for stiff kinetics: the catalog, Jacobian reuse, KPP's integrator, 2015 to 2026 |
| `R2-accuracy-requirements.md` | how accurate a chemistry solver must be, and how the solver comparisons measured and concluded |
| `R3-extrapolation-exponential-sdirk.md` | extrapolated linearly implicit Euler, exponential integrators, and SDIRK/ESDIRK, against Rosenbrock at low accuracy |
| `R4-controllers-norms-reuse-positivity.md` | step-size controllers, error norms and tolerances, Jacobian and LU reuse, restarts, positivity |
| `R5-alternative-schemes-precision.md` | TWOSTEP, EBI, QSSA, BDF/NDF, IMEX and partitioned schemes, reduced precision |
| `R6-ts1-camchem-musica-sda.md` | the TS1 mechanism, the CAM-chem/MOZART/MUSICA solvers and tolerances, the origin of the significant-digits metric, the size of other errors |

What the reviews settled for the design is summarized at the head of
`../REFERENCES.md`; the decisions they led to are in `../DECISIONS.md`.
