# Decisions

- The implementation uses double precision and C11. The numerical formulas are
  ported from the specified reference revision; no reference code is linked into
  the C library or benchmark.
- Configuration reading and code generation run only at build time, in the local
  prefix environment. JSON and YAML share one normalized representation.
- The TS1 translation preserves reaction order and floating-point constants.
  Extensions to CAMP configuration: species `third body`, `initial concentration`,
  mechanism `initial conditions`, and `initial rate parameters`. These represent
  the reference's parameterized air species and InitState defaults explicitly.
- The user corrected the delivery instruction: commit locally and do not push.
  The supplied repository has no remote.
- PROMPT.md is an experiment input and is ignored rather than committed because
  it contains identifying names forbidden in deliverable files.
