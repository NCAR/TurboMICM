"""Build a small adversarial mechanism for regeneration and algebra tests.

Start with the reader's extra-rate fixture, change a fixed parameter, then add
new species, repeated reactants, an emission, a parameterized reactant, and a
shared user parameter. The same generated builder feeds the generic oracle.
"""
import json
from pathlib import Path
import sys

source = json.loads(Path(sys.argv[1]).read_text())
data = source['camp-data']
data.extend([{'type': 'CHEM_SPEC', 'name': 'extra gas'},
             {'type': 'CHEM_SPEC', 'name': 'density', 'parameterization': 'air density'}])
mechanism = next(r for r in data if r['type'] == 'MECHANISM')
mechanism['reactions'][0]['k0_A'] = .017
names = [r['name'] for r in data if r['type'] == 'CHEM_SPEC' and not r.get('parameterization')]
mechanism['reactions'].extend([
    {'type': 'ARRHENIUS', 'A': .03, 'B': .2, 'C': -14, 'rate units': 'mol m-3',
     'reactants': [names[0], names[0], 'density'],
     'products': [{'species': 'extra gas', 'yield': 1.25}, {'species': names[0], 'yield': .1}]},
    {'type': 'USER_DEFINED', 'label': 'shared', 'scaling': .7,
     'reactants': {'extra gas': {}}, 'products': {names[0]: {'yield': .5}}},
    {'type': 'USER_DEFINED', 'label': 'shared', 'scaling': .01,
     'reactants': {}, 'products': {'extra gas': {}}},
])
Path(sys.argv[2]).write_text(json.dumps(source, indent=2) + '\n')
