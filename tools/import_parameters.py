"""Refresh the C parameter tables from the supplied numerical-method specification.

This is a development-time transcription helper, not part of a normal build.
"""
from pathlib import Path
import re
import sys

source = Path(sys.argv[1]).read_text()
methods = [
    ('TwoStageRosenbrock', 'CHEM_ROS2'),
    ('ThreeStageRosenbrock', 'CHEM_ROS3'),
    ('FourStageRosenbrock', 'CHEM_ROS4'),
    ('FourStageDifferentialAlgebraicRosenbrock', 'CHEM_RODA3'),
    ('SixStageDifferentialAlgebraicRosenbrock', 'CHEM_RODA4'),
]
output = ['/* Parameter sets transcribed from MICM/UCAR, Apache-2.0. */',
          '#include "chem.h"', '#include <float.h>', '#include <math.h>', '',
          'chem_ros_parameters chem_ros_parameters_default(chem_ros_method method)', '{',
          '    chem_ros_parameters p = {0};',
          '    p.max_steps=1000; p.round_off=DBL_EPSILON;',
          '    p.factor_min=0.2; p.factor_max=6; p.rejection_factor=0.1; p.safety_factor=0.9;',
          '    p.init_max_iterations=10; p.init_tolerance=0.1; p.init_max_backtracks=24;',
          '    p.init_backtrack_factor=0.5; p.init_sufficient_decrease=1e-4;',
          '    switch (method) {']
for name, enum in methods:
    body = source.split('::' + name + 'Parameters()')[1].split('return parameters;')[0]
    output.append(f'    case {enum}: {{')
    for line in body.splitlines():
        line = line.strip()
        if line.startswith('Real g ='):
            output.append('        double g = 1.0 + 1.0 / sqrt(2.0);')
        if not line.startswith('parameters.'):
            continue
        line = re.sub(r'parameters\.(\w+)_', r'p.\1', line)
        line = line.replace('new_function_evaluation', 'new_function').replace('estimator_of_local_order','order')
        line = line.replace('true','1').replace('false','0')
        if '.fill(' in line:
            if '.fill(1)' in line:
                output.append('        for (size_t i=0; i<6; ++i) p.new_function[i]=1;')
            continue
        output.append('        ' + line)
    output.extend(['        break;', '    }'])
output.extend(['    }', '    return p;', '}', ''])
Path(sys.argv[2]).write_text('\n'.join(output))
