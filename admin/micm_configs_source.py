#!/usr/bin/env python3
"""Write a copy of reference/micm_ts1_ref.cpp whose MICM configuration is chosen at compile time.

    admin/micm_configs_source.py REFERENCE_SOURCE OUTPUT

The copy replaces the reference's solver type aliases (micm_bench's default:
standard ordering, in-place Mozart LU) with macros for micm_bench's other CPU
configurations, and changes nothing else:
    -DTS1_VECTOR=L      vector ordering, L cells interleaved (micm_bench: 1, 2, 4, 8, 128)
    -DTS1_LU_SEPARATE   separate L and U factors (CpuSolverBuilder), not in-place
    -DTS1_LU_DOOLITTLE  Doolittle LU, not Mozart
Built without them, the copy is the reference program. The reference itself is
left alone, because it is one of the agents' read-only inputs.
"""

import sys
from pathlib import Path

DEFAULT = """  using DenseMatrix = micm::Matrix<micm::Real>;
  using SparseMatrix = micm::SparseMatrix<micm::Real, micm::SparseMatrixStandardOrdering>;
  using SolverBuilder = micm::CpuSolverBuilderInPlace<
      micm::RosenbrockSolverParameters,
      DenseMatrix,
      SparseMatrix,
      micm::LuDecompositionMozartInPlace<SparseMatrix>>;
"""


def builder(kind, lu):
    return f"""  using SolverBuilder = micm::{kind}<
      micm::RosenbrockSolverParameters,
      DenseMatrix,
      SparseMatrix,
      micm::{lu}<SparseMatrix>>;
"""


CONFIGURABLE = (
    """  // Chosen at compile time (admin/micm_configs_source.py); the default is micm_bench's.
#ifdef TS1_VECTOR
  using DenseMatrix = micm::VectorMatrix<micm::Real, TS1_VECTOR>;
  using SparseMatrix = micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<TS1_VECTOR>>;
#else
  using DenseMatrix = micm::Matrix<micm::Real>;
  using SparseMatrix = micm::SparseMatrix<micm::Real, micm::SparseMatrixStandardOrdering>;
#endif
#if defined(TS1_LU_SEPARATE) && defined(TS1_LU_DOOLITTLE)
"""
    + builder("CpuSolverBuilder", "LuDecompositionDoolittle")
    + "#elif defined(TS1_LU_SEPARATE)\n"
    + builder("CpuSolverBuilder", "LuDecompositionMozart")
    + "#elif defined(TS1_LU_DOOLITTLE)\n"
    + builder("CpuSolverBuilderInPlace", "LuDecompositionDoolittleInPlace")
    + "#else\n"
    + builder("CpuSolverBuilderInPlace", "LuDecompositionMozartInPlace")
    + "#endif\n"
)

source, output = Path(sys.argv[1]), Path(sys.argv[2])
text = source.read_text()
if text.count(DEFAULT) != 1:
    sys.exit(f"{source}: the solver type aliases are not in the expected form; update {Path(__file__).name}")
output.write_text(text.replace(DEFAULT, CONFIGURABLE))
print(f"wrote {output}")
