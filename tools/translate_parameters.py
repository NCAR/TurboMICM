"""One-time, auditable translation of the supplied numerical parameter tables."""
import argparse
import re
from pathlib import Path


def translate(source):
    text = Path(source).read_text()
    names = ["TwoStage", "ThreeStage", "FourStage", "FourStageDifferentialAlgebraic", "SixStageDifferentialAlgebraic"]
    result = ['/* Numerical coefficients: Copyright (C) 2023-2026 University Corporation for Atmospheric Research.',
              ' * SPDX-License-Identifier: Apache-2.0 */', '#include "internal.h"',
              'ChemRosenbrock chem_rosenbrock_parameters(int method) {',
              '  ChemRosenbrock p = {0};',
              '  p.max_steps = 1000; p.roundoff = DBL_EPSILON;',
              '  p.factor_min = 0.2; p.factor_max = 6.0; p.rejection_factor = 0.1; p.safety = 0.9;',
              '  p.init_iterations = 10; p.init_backtracks = 24; p.init_tolerance = 0.1;',
              '  p.init_backtrack_factor = 0.5; p.init_sufficient_decrease = 1e-4;',
              '  switch (method) {']
    for method, name in zip([2, 3, 4, 5, 6], names):
        body = re.search(r'inline RosenbrockSolverParameters RosenbrockSolverParameters::' + name + r'RosenbrockParameters\(\)\s*\{(.*?)return parameters;', text, re.S).group(1)
        body = re.sub(r'//[^\n]*', '', body)
        result.append(f'    case {method}: {{')
        for stmt in body.split(';'):
            stmt = stmt.strip()
            if not stmt or stmt == 'RosenbrockSolverParameters parameters':
                continue
            stmt = stmt.replace('Real ', 'double ').replace('std::numbers::sqrt2', 'sqrt(2.0)')
            stmt = stmt.replace('parameters.', 'p.').replace('estimator_of_local_order_', 'order').replace('new_function_evaluation_', 'new_function')
            stmt = re.sub(r'\b(stages|a|c|m|e|alpha|gamma)_', r'\1', stmt)
            stmt = stmt.replace('true', '1').replace('false', '0')
            fill = re.fullmatch(r'p\.(\w+)\.fill\(([^)]+)\)', stmt)
            if fill:
                field, val = fill.groups()
                result.append(f'      for (size_t j = 0; j < sizeof p.{field}/sizeof *p.{field}; ++j) p.{field}[j] = {val};')
            else:
                result.append('      ' + stmt + ';')
        result.append('      break; }')
    result += ['    default: p.stages = 0; break;', '  }', '  return p;', '}']
    return '\n'.join(result) + '\n'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source')
    parser.add_argument('output')
    args = parser.parse_args()
    Path(args.output).write_text(translate(args.source))
