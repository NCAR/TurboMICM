"""One-time, retained numerical-table translation, not used by normal builds."""
import pathlib
import re
import sys

source = pathlib.Path(sys.argv[1]).read_text()
methods = [(2, "TwoStageRosenbrockParameters"), (3, "ThreeStageRosenbrockParameters"),
           (4, "FourStageRosenbrockParameters"), (43, "FourStageDifferentialAlgebraicRosenbrockParameters"),
           (6, "SixStageDifferentialAlgebraicRosenbrockParameters")]
result = ['/* Numerical parameter tables translated from MICM (Apache-2.0). */',
          '#include "chem.h"', '#include <float.h>', '#include <math.h>', '#include <string.h>',
          'ck_status ck_rosenbrock_parameters_init(int method, ck_rosenbrock_parameters *p)', '{',
          '  if (!p) return CK_INVALID;', '  memset(p, 0, sizeof *p);',
          '  p->max_steps = 1000; p->round_off = DBL_EPSILON;',
          '  p->factor_min = 0.2; p->factor_max = 6; p->rejection_decrease = 0.1; p->safety = 0.9;',
          '  switch (method) {']
for number, name in methods:
    body = source.split('inline RosenbrockSolverParameters RosenbrockSolverParameters::' + name + '()')[1]
    body = body.split('return parameters;')[0]
    result.append(f'  case {number}: {{')
    if number == 2:
        result.append('    double g = 1.0 + 1.0 / sqrt(2.0);')
    for line in body.splitlines():
        line = line.strip()
        if not line.startswith('parameters.'):
            continue
        line = line.replace('parameters.', 'p->')
        line = re.sub(r'([a-z_]+)_', lambda m: m[0], line)  # replace suffixes below, not internal underscores
        line = re.sub(r'p->(\w+)_', lambda m: 'p->' + m[1], line)
        line = line.replace('new_function_evaluation', 'new_function').replace('estimator_of_local_order', 'order')
        match = re.fullmatch(r'p->(\w+)\.fill\(([^)]+)\);', line)
        if match:
            field, value = match.groups()
            length = 15 if field in ('a', 'c') else 6
            if value not in ('0', '0.0', 'false'):
                result.append(f'    for (size_t j = 0; j < {length}; ++j) p->{field}[j] = {value.replace("true", "1")};')
        else:
            result.append('    ' + line.replace('true', '1').replace('false', '0'))
    result.append('    return CK_OK; }')
result += ['  default: return CK_INVALID;', '  }', '}']
pathlib.Path(sys.argv[2]).write_text('\n'.join(result) + '\n')
