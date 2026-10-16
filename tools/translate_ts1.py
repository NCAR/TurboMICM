"""Transcribe the supplied TS1 C++ data into the configuration format.

No C++ is compiled or executed. Reject an unrecognized declaration rather than
silently omitting a species, reaction, parameter, or initial condition.
"""
from collections import Counter
import json
from pathlib import Path
import re
import sys


def translate(path):
    text = Path(path).read_text()
    species = []
    for name, body, diffusion in re.findall(
        r'auto species = micm::Species\("([^"]+)"\);(.*?)phase_species.emplace_back\(species(?:,\s*([^)]*))?\);',
        text, re.S,
    ):
        sp = {'name': name, 'type': 'CHEM_SPEC'}
        if 'SetThirdBody()' in body:
            sp['parameterization'] = 'air density'
        weight = re.search(r'MOLECULAR_WEIGHT, micm::Real\{\s*([^}]+)', body)
        if weight:
            sp['molecular weight [kg mol-1]'] = float(weight[1])
        if diffusion:
            sp['diffusion coefficient [m2 s-1]'] = float(diffusion)
        species.append(sp)
    types = {'Arrhenius': 'ARRHENIUS', 'Troe': 'TROE', 'UserDefined': 'USER_DEFINED', 'Surface': 'SURFACE'}
    reactions = []
    for block in text.split('processes.push_back(')[1:]:
        block = block.split('.Build());', 1)[0]
        reactants = re.search(r'\.SetReactants\(\s*\{(.*?)\}\)', block, re.S)
        products = re.search(r'\.SetProducts\(\s*\{(.*?)\}\)', block, re.S)
        rate = re.search(r'\.SetRateConstant\(\s*micm::(\w+)RateConstantParameters\{(.*?)\}\)', block, re.S)
        if not (reactants and products and rate and rate[1] in types):
            raise ValueError('Unrecognized reaction: ' + block[:100])
        rx = {'type': types[rate[1]],
              'reactants': re.findall(r'SpeciesByName\(gas_phase, "([^"]+)"\)', reactants[1]),
              'products': [{'species': n, 'yield': float(y)} for n, y in re.findall(
                  r'StoichSpecies\{\s*SpeciesByName\(gas_phase, "([^"]+)"\),\s*([^}]+)\}', products[1])],
              'concentration units': 'mol m-3'}
        for key, value in re.findall(r'\.(\w+)_\s*=\s*("[^"]*"|[-+\d.eE]+)', rate[2]):
            rx[{'label':'MUSICA name','scaling_factor':'scaling factor',
                'reaction_probability':'reaction probability'}.get(key, key)] = json.loads(value) if value.startswith('"') else float(value)
        if rx['type'] == 'SURFACE':
            rx['gas-phase reactant'] = rx['reactants'][0]
        reactions.append(rx)
    initial = dict((n,float(v)) for n,v in re.findall(r'concentration\("([^"]+)",\s*([^)]+)\);',text))
    params = dict((n,float(v)) for n,v in re.findall(r'parameter\("([^"]+)",\s*([^)]+)\);',text))
    for sp in species:
        if sp['name'] in initial:
            sp['initial concentration [mol m-3]'] = initial.pop(sp['name'])
    if initial or len(species)!=210 or len(reactions)!=547:
        raise ValueError(f'Incomplete translation: {len(species)} species, {len(reactions)} reactions, unmatched {initial}')
    conditions = {}
    for field in ('temperature', 'pressure'):
        matches = re.findall(r'conditions_\[c\]\.'+field+r'_\s*=\s*([^;]+);',text)
        if len(matches)!=1:
            raise ValueError('Missing default condition: ' + field)
        conditions[field] = float(matches[0])
    return {'camp-data': species + [{'name':'TS1', 'type':'MECHANISM', 'reactions':reactions}],
            'initial conditions':conditions, 'initial rate parameters':params}


if __name__ == '__main__':
    result = translate(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(result,indent=2)+'\n')
    print(f"Translated {len(result['camp-data'])-1} species and {len(result['camp-data'][-1]['reactions'])} reactions")
