# Decisions

- The C library uses a mechanism descriptor with contiguous reaction arrays and
  a separate solver context per cell group. This keeps mechanism data immutable
  and lets the benchmark use the same construction path as library clients.
- YAML and JSON input are normalized in Python at build time. The generated C
  contains only standard C data and the library has no Python or parser runtime
  dependency.
- TS1's source builder is the translation input because the task explicitly
  provides it as the source of truth. The kept translation script emits the
  same YAML/JSON configuration shape that the public reader accepts.
- The config extension needed by TS1 is a parameterized third-body species
  (`air_density`), which participates in pressure-dependent rates but is
  omitted from the integrated state. User-defined rates also require named
  rate parameters, and surface rates require diffusion and aerosol parameters.
