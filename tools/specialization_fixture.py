"""Build-time regression: add a reaction and change a TS1 kinetic parameter."""
import json
import sys
from pathlib import Path

data = json.loads(Path(sys.argv[1]).read_text())
mechanism = next(r for r in data['camp-data'] if r['type'] == 'MECHANISM')
first = next(r for r in mechanism['reactions'] if r['type'] == 'ARRHENIUS')
first['A'] *= 1.03125
mechanism['reactions'].append(dict(type='USER_DEFINED', reactants={'O3': {}}, products={'O': {}},
                                  **{'MUSICA name': 'added rate', 'scaling factor': .75}))
data['initial state']['rate parameters']['added rate'] = 1e-8
Path(sys.argv[2]).write_text(json.dumps(data,indent=2)+'\n')
