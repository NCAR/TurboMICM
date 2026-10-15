# Architecture

The C library owns mechanism descriptions, rate calculations, process forcing
and analytic Jacobians, sparse symbolic/numeric LU, solver parameters, state,
and adaptive integration. Configuration tooling reads CAMP YAML/JSON and emits
C descriptions and optional kernels. Both descriptions built in C and generated
from configuration use the same public library API.

The benchmark owns CSV handling and partitions cells among independent solver
states. OpenMP parallelizes groups; species-major storage within each group
allows portable SIMD across cells. Group width and thread count are parameters.

Detailed extension instructions and measured defaults will be recorded alongside
implementation and validation.
