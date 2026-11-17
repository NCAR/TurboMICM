#!/usr/bin/env python3
"""Translate the generated MICM TS1 C++ description to the configuration format.

This deliberately recognizes the source's restricted builder syntax, not arbitrary
C++. Counts and unconsumed constructs are checked so a changed source fails loudly.
"""
import argparse
import json
import re
from collections import Counter
from pathlib import Path

NUMBER = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'


def translate(source):
    phase = source.split('inline micm::Phase CreateGasPhase()', 1)[1].split('return micm::Phase', 1)[0]
    species = []
    for block in re.findall(r'\{\s*auto species = micm::Species\("([^"]+)"\);(.*?)\n\s*\}', phase, re.S):
        name, body = block
        item = {'type': 'CHEM_SPEC', 'name': name}
        weight = re.search(r'MOLECULAR_WEIGHT, micm::Real\{\s*(' + NUMBER + r')\s*\}', body)
        if weight:
            item['molecular weight [kg mol-1]'] = float(weight.group(1))
        diffusion = re.search(r'phase_species\.emplace_back\(species,\s*(' + NUMBER + r')\)', body)
        if diffusion:
            item['diffusion coefficient [m2 s-1]'] = float(diffusion.group(1))
        if 'SetThirdBody()' in body:
            item['third body'] = True
        species.append(item)
    if len(species) != 210:
        raise ValueError(f'expected 210 species, found {len(species)}')

    section = source.split('inline std::vector<micm::Process> CreateProcesses', 1)[1].split('return processes;', 1)[0]
    blocks = re.findall(r'processes\.push_back\((.*?)\.Build\(\)\);', section, re.S)
    if len(blocks) != 547 or section.count('processes.push_back(') != len(blocks):
        raise ValueError(f'expected 547 reactions, found {len(blocks)}')
    reactions = []
    for block in blocks:
        react = re.search(r'\.SetReactants\(\s*\{(.*?)\}\s*\)', block, re.S)
        prod = re.search(r'\.SetProducts\(\s*\{(.*?)\}\s*\)', block, re.S)
        rate = re.search(r'\.SetRateConstant\(\s*micm::(\w+)RateConstantParameters\s*\{(.*?)\}\s*\)', block, re.S)
        if not (react and prod and rate):
            raise ValueError('unrecognized reaction: ' + block[:200])
        names = re.findall(r'SpeciesByName\(gas_phase,\s*"([^"]+)"\)', react.group(1))
        quantities = Counter(names)
        reactants = {n: ({} if q == 1 else {'qty': q}) for n, q in quantities.items()}
        pairs = re.findall(r'StoichSpecies\{\s*SpeciesByName\(gas_phase,\s*"([^"]+)"\),\s*(' + NUMBER + r')\s*\}', prod.group(1))
        if len(pairs) != prod.group(1).count('StoichSpecies{') or len(dict(pairs)) != len(pairs):
            raise ValueError('unrecognized or duplicate products: ' + prod.group(1))
        products = {n: {'yield': float(v)} for n, v in pairs}
        kind = {'Arrhenius': 'ARRHENIUS', 'Troe': 'TROE', 'UserDefined': 'USER_DEFINED', 'Surface': 'SURFACE'}[rate.group(1)]
        item = {'type': kind, 'reactants': reactants, 'products': products}
        attrs = rate.group(2)
        for key, value in re.findall(r'\.(\w+)_\s*=\s*("[^"]*"|' + NUMBER + r'|PhaseSpeciesByName\(gas_phase,\s*"[^"]+"\))', attrs):
            if key == 'phase_species':
                item['gas-phase reactant'] = re.search(r'"([^"]+)"', value).group(1)
            elif key == 'label':
                item['MUSICA name'] = json.loads(value)
            else:
                item[{'scaling_factor': 'scaling factor', 'reaction_probability': 'reaction probability'}.get(key, key)] = float(value)
        if len(re.findall(r'\.\w+_\s*=', attrs)) != len(re.findall(r'\.(\w+)_\s*=\s*("[^"]*"|' + NUMBER + r'|PhaseSpeciesByName\(gas_phase,\s*"[^"]+"\))', attrs)):
            raise ValueError('unrecognized rate parameters: ' + attrs)
        reactions.append(item)
    defaults = {}
    for kind, name, value in re.findall(r'\b(concentration|parameter)\("([^"]+)",\s*(' + NUMBER + r')\);', source):
        if name in defaults:
            raise ValueError('duplicate default: ' + name)
        defaults[name] = float(value)
    environment = {}
    for field in ('temperature', 'pressure'):
        matches = re.findall(r'state\.conditions_\[c\]\.' + field + r'_\s*=\s*(' + NUMBER + r')\s*;', source)
        if len(matches) != 1:
            raise ValueError(f'expected one initial {field}, found {len(matches)}')
        environment[field] = float(matches[0])
    return species, reactions, defaults, environment


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    species, reactions, defaults, environment = translate(args.source.read_text())
    args.output.mkdir(parents=True, exist_ok=True)
    for name, obj in [('config.json', {'camp-files': ['species.json', 'reactions.json'], 'rate units': 'mol m-3', 'initial conditions': defaults, 'initial environment': environment}),
                      ('species.json', {'camp-data': species}),
                      ('reactions.json', {'camp-data': [{'type': 'MECHANISM', 'name': 'TS1', 'reactions': reactions}]})]:
        (args.output / name).write_text(json.dumps(obj, indent=2) + '\n')
    print(f'translated {len(species)} species, {len(reactions)} reactions, {len(defaults)} defaults')


if __name__ == '__main__':
    main()
