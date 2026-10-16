# Decisions

- The public C interface uses owned model and state handles, explicit error codes,
  and callback tables in place of C++ builders, exceptions, and policy templates.
  All numerical storage is double precision. Mechanisms built in C and mechanisms
  read from configuration use the same model finalization and solver paths.
- A group stores each species contiguously across cells. OpenMP distributes groups;
  `omp simd` expresses independent cell arithmetic. Each group has independent
  adaptive step control and private work arrays.
- Configuration parsing and optional C generation happen only at build time, in
  the local conda environment. The native library has no Python dependency.
- The original prompt stays an untracked input because it contains experiment
  identities. Generated and authored deliverables contain no participant identity.
