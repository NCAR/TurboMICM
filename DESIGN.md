# Architecture

The mechanism owns species, reaction, rate-constant, and constraint data. Finalization establishes the stable reaction order, custom-parameter labels, sparse Jacobian structure, fill pattern, and sparse factorization/substitution plans. Mutable state and solver workspaces belong to each cell group; the finalized mechanism is shared read-only across threads.

Species and sparse matrix elements are stored in separate contiguous arrays over cells. OpenMP SIMD loops traverse cells; outer OpenMP loops distribute cell groups. The group size is currently 32 by default, configurable through `CHEM_GROUP_SIZE` at build time or `--group N` at run time. The thread count defaults to `omp_get_max_threads()` and is configurable with `--threads N` or the OpenMP environment.

`tools/configure.py` reads the example configuration format in either YAML or JSON, resolving `camp-files` relative to each source. It emits mechanism construction through the public C API and specialized forcing/Jacobian kernels. The same mechanism can select the interpreted kernels, which provides a direct check of specialization. The benchmark selects generated kernels.

The Rosenbrock solver receives a parameter table, a step-size controller function, and factor/solve function pointers. A parameter set can be added to `src/parameters.c` without changing the integration loop. The default controller is isolated in `chem_default_controller`; callers may replace it and pass their own context. Linear algebra defaults to in-place Mozart elimination; Doolittle is also available through the factor callback. A new rate type adds an enum entry, parameter fields, and a branch in `chem_rate_value`, then parser support if it is used in configuration files.

The C API is declared in `include/chem.h`. All model storage is owned by the model, all state storage by the state, and states must be destroyed before their model. Reaction and constraint insertion copies caller-owned arrays. The model is immutable after finalization. Setters return error codes for invalid names and sizes.
