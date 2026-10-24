## Findings for the SSE

These bear on the primary gate or on the record, not on the scores above.
Each was confirmed by the judge in the code.

- **luna-high's 662 PASS lines are printed from a list.** `tools/test_cases.py`
  runs no check; it prints `PASS` with a prefix-chosen reason for each of 510
  literal ids once `bin/test_core`, a single program that returns at its first
  failing check, has exited 0 (`tools/test_cases.py:7-639,640-778`;
  `Makefile`, target `test-cpu`). Of the twelve sampled ids, one
  (`chapman_ode_solver:…factored_alpha_minus_jac`) has no corresponding check in
  `tests/test_core.c`; the others have checks of varying fidelity that gate the
  program as a whole.
- **luna-max reports some ids from checks that are not theirs.** All 52 LU and
  linear-solver ids run one dense check (`tests/test_cpu.c:4313-4316`); the four
  Troe ids evaluate one fixed Troe case with none of MICM's parameters
  (`:21-48,113-115`); `external_model_constraints:…DAESolveEnforcesConservation`
  has no branch and falls to an unrelated constraint step (`:3332-3426`).
- **luna-xhigh takes one binary's exit status for 30 rate-constant ids and the
  reaction-rate store ids** (`tools/test_cases.py:74,145-150`).
- **All three Luna entries accept a step whose error norm is NaN.** The norm
  floors with `fmax(…, 1e-10)`, which returns 1e-10 for a NaN argument, so the
  finiteness test after it never fires and the candidate is accepted
  (luna-max `src/solver.c:80,767-776`; luna-xhigh `src/chem.c:1045,1206`;
  luna-high `src/solver.c:907,916`). sol-high floors the same way but inspects
  every candidate value (`src/rosenbrock.c:98-101`); the other five propagate
  the NaN through a conditional or a custom maximum. TS1 does not reach this.
- **luna-high's generator names the wrong species for CAMP inputs.**
  `tools/read_config.py:83-86` records each term's input position and
  `tools/generate_c.py:102,126` uses it as the index into a species table it has
  permuted. The E5 helper's run on MICM's user-defined tutorial with the species
  listed in reverse order emitted all ten reactions with wrong species
  (`exercises/e5-check/run.log`). TS1 is translated by `tools/import_ts1.py`
  without positions and is unaffected.
- **The three Luna readers convert units with fixed exponents.** Troe and
  ternary coefficients are scaled by fixed powers of the molecule-per-cm3
  factor, and Arrhenius, tunneling and branched coefficients are not scaled at
  all (luna-max `tools/generate_c.py:88-94`; luna-xhigh
  `tools/generate.py:131-136`; luna-high `tools/read_config.py:75-79`). A
  bimolecular Arrhenius rate given in molecule cm-3 units stays as written: the
  E5 helper's probe shows luna-high leaving A = 1e-12 for an A + B reaction.
  MICM's tutorials happen to need only the Troe and ternary conversions.
- **luna-max's backward Euler is not MICM's.** One Newton loop over the whole
  interval per cell, with a dense LU; on non-convergence it restores the state
  and returns -1; no step reduction sequence (`src/backward_euler.c:35-117`).
- **Generated kernels and the run-time plan.** astra-xhigh installs generated
  kernels whose flat indices come from a Python copy of the ordering and fill,
  with no check against the model's own (`src/model.c:201-202`;
  `tools/generate.py:13-34,61-75`); sol-max checks only that the nonzero
  counts agree (`tools/configure.py:230`). Both agree with MICM on TS1 and on
  the tutorials; the hazard is for a new mechanism.
- **astra-high recomputes every rate constant inside the timed loop**
  (`src/ts1_bench.c:154`), work that MICM's driver and the other entries do
  before the warm-up; its TS1 time therefore includes it.
- **luna-xhigh builds its symbolic plan and allocates all work arrays inside
  the integrator call**, once per group inside the benchmark's timer
  (`src/chem.c:1080-1098`, `src/ts1_bench.c:364-365`).
