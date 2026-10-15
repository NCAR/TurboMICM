"""Translate the supplied TS1 C++ data declarations into CAMP JSON.

This is a data importer, not a C++ compiler. It accepts the declaration forms
in the supplied input and rejects incomplete extraction. No source is executed.
"""
import argparse
import json
from pathlib import Path
import re

NUMBER = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'

def translate(path):
    source = Path(path).read_text()
    species_section, process_section = source.split('inline std::vector<micm::Process> CreateProcesses', 1)
    process_section, init_section = process_section.split('void InitState', 1)
    species = []
    for name, body in re.findall(r'auto species = micm::Species\("([^"]+)"\);(.*?)phase_species.emplace_back\(species[^;]*;', species_section, re.S):
        entry = {'name': name, 'type': 'CHEM_SPEC'}
        weight = re.search(r'MOLECULAR_WEIGHT, micm::Real\{\s*(' + NUMBER + ')', body)
        if weight:
            entry['molecular weight [kg mol-1]'] = float(weight[1])
        if 'SetThirdBody' in body:
            entry['parameterization'] = 'air density'
        species.append(entry)
    # The diffusion value is on emplace_back, outside the body above.
    for name, diffusion in re.findall(r'auto species = micm::Species\("([^"]+)"\);(?:(?!auto species).)*?phase_species.emplace_back\(species,\s*(' + NUMBER + r')\)', species_section, re.S):
        next(s for s in species if s['name'] == name)['diffusion coefficient [m2 s-1]'] = float(diffusion)
    reactions = []
    names = r'SpeciesByName\(gas_phase, "([^"]+)"\)'
    types = {'Arrhenius': 'ARRHENIUS', 'Troe': 'TROE', 'UserDefined': 'USER_DEFINED', 'Surface': 'SURFACE'}
    for chunk in process_section.split('processes.push_back(')[1:]:
        chunk = chunk.split('.Build());', 1)[0]
        lhs = re.search(r'\.SetReactants\(\s*\{(.*?)\}\)', chunk, re.S)
        rhs = re.search(r'\.SetProducts\(\s*\{(.*?)\}\)', chunk, re.S)
        rate = re.search(r'\.SetRateConstant\(\s*micm::(\w+)RateConstantParameters\{(.*?)\}\)', chunk, re.S)
        if not (lhs and rhs and rate):
            raise ValueError('Unrecognized reaction declaration')
        entry = {'type': types[rate[1]], 'reactants': re.findall(names, lhs[1]),
                 'products': [{'species': n, 'yield': float(v)} for n, v in re.findall(names + r',\s*(' + NUMBER + ')', rhs[1])]}
        for key, value in re.findall(r'\.(\w+)_\s*=\s*(' + NUMBER + ')', rate[2]):
            key = {'scaling_factor': 'scaling factor', 'reaction_probability': 'reaction probability'}.get(key, key)
            entry[key] = float(value)
        label = re.search(r'\.label_\s*=\s*"([^"]+)"', rate[2])
        if label:
            entry['MUSICA name'] = label[1]
        if entry['type'] == 'SURFACE':
            phase = re.search(r'PhaseSpeciesByName\(gas_phase, "([^"]+)"\)', rate[2])
            if not phase or entry['reactants'] != [phase[1]]:
                raise ValueError('Unexpected surface reactants')
            entry['gas-phase reactant'] = phase[1]
            entry['gas-phase products'] = entry.pop('products')
            del entry['reactants']
        reactions.append(entry)
    initial = {}
    for label, value in re.findall(r'(?:concentration|parameter)\("([^"]+)",\s*(' + NUMBER + r')\)', init_section):
        if label in initial:
            raise ValueError(f'Duplicate initial value {label}')
        initial[label] = float(value)
    for key in ('temperature', 'pressure'):
        initial[key] = float(re.search(r'\.' + key + r'_\s*=\s*(' + NUMBER + ')', init_section)[1])
    if len(species) != 210 or len(reactions) != 547:
        raise ValueError(f'Incomplete translation: {len(species)} species, {len(reactions)} reactions')
    species_names = {s['name'] for s in species if not s.get('parameterization')}
    if not species_names <= initial.keys():
        raise ValueError('Missing initial species concentration')
    return {'camp-data': species + [{'name': 'TS1', 'type': 'MECHANISM', 'rate units': 'mol m-3',
                                    'reactions': reactions, 'initial state': initial}]}

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=Path)
    p.add_argument('output', type=Path)
    args = p.parse_args()
    args.output.write_text(json.dumps(translate(args.source), indent=2) + '\n')
