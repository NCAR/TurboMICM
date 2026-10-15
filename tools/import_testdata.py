#!/usr/bin/env python3
"""Retain the published reference vectors used by MICM's analytical cases.

One-time fixture translation only: normal builds never open MICM sources.
"""
from pathlib import Path
import re
source = Path('../micm/test/integration/analytical_policy.hpp').read_text()
output = ['/* Reference vectors from MICM analytical_policy.hpp; see NOTICE. */']
for name, cols in [('Robertson',3),('Oregonator',5),('Hires',8),('E5',6)]:
    body = source[source.index('void TestAnalytical'+name+'('):]
    start = re.search(r'analytical_concentrations\s*=\s*\{',body).end()-1
    depth=0
    for end in range(start,len(body)):
        depth += (body[end]=='{') - (body[end]=='}')
        if depth==0:
            break
    literal=body[start:end+1]
    output.append(f'static const double expected_{name.lower()}[][{cols}] = '+literal+';')
Path('tests/analytical_expected.h').write_text('\n\n'.join(output)+'\n')
