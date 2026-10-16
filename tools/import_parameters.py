"""Translate MICM's five Rosenbrock coefficient tables into reviewable C data.

This maintainer tool is not needed by the build. The generated source is committed.
"""
from pathlib import Path
import re
import sys

source = Path(sys.argv[1]).read_text()
methods = [("TwoStageRosenbrock", 2), ("ThreeStageRosenbrock", 3),
           ("FourStageRosenbrock", 4), ("FourStageDifferentialAlgebraicRosenbrock", 43),
           ("SixStageDifferentialAlgebraicRosenbrock", 6)]
names = {"new_function_evaluation": "new_function", "estimator_of_local_order": "estimator_order"}
out = ['/* Rosenbrock coefficient tables ported from MICM.',
       ' * Copyright (C) 2023-2026 University Corporation for Atmospheric Research',
       ' * SPDX-License-Identifier: Apache-2.0 */', '#include "chem.h"', '#include <float.h>', '#include <math.h>',
       'ch_rosenbrock_parameters ch_rosenbrock_default(int method) {',
       '  ch_rosenbrock_parameters p = {0};',
       '  p.max_steps = 1000; p.round_off = DBL_EPSILON;',
       '  p.factor_min = 0.2; p.factor_max = 6.0; p.rejection_factor = 0.1; p.safety_factor = 0.9;',
       '  p.init_max_iterations = 10; p.init_tolerance = 0.1; p.init_max_backtracks = 24;',
       '  p.init_backtrack_factor = 0.5; p.init_sufficient_decrease = 1e-4;',
       '  p.controller = ch_default_controller;', '  switch (method) {']
for name, number in methods:
    body = re.search(rf'::{name}Parameters\(\)\s*\{{(.*?)return parameters;', source, re.S)[1]
    out.append(f'    case {number}: {{')
    if number == 2:
        out.append('      double g = 1.0 + 1.0 / sqrt(2.0);')
    for line in body.splitlines():
        line = line.strip()
        if not line.startswith('parameters.'):
            continue
        line = re.sub(r'parameters\.([a-z_]+)_', lambda m: 'p.' + names.get(m[1], m[1]), line)
        if '.fill(' in line:
            if line == 'p.new_function.fill(true);':
                out.append('      for (int i = 0; i < 6; ++i) p.new_function[i] = 1;')
            continue
        out.append('      ' + line.replace('true', '1').replace('false', '0'))
    out.append('      break;\n    }')
out.extend(['    default: p.stages = 0; break;', '  }', '  return p;', '}'])
Path(sys.argv[2]).write_text('\n'.join(out) + '\n')
