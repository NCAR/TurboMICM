"""Translate the supplied TS1 C++ data declarations to the configuration format.

This reads declarations as text; it never compiles or executes reference code.
Unknown constructs fail the translation instead of being silently dropped.
"""
import argparse
import json
import re
from pathlib import Path

NUMBER = r'[-+0-9.eE]+'


def translate(path):
    text = Path(path).read_text()
    species = []
    for match in re.finditer(r'auto species = micm::Species\("([^"]+)"\);(.*?)phase_species.emplace_back\(species(?:,\s*(' + NUMBER + r'))?\);', text, re.S):
        name, body, diffusion = match.groups()
        record = {'name': name, 'type': 'CHEM_SPEC'}
        weight = re.search(r'MOLECULAR_WEIGHT, micm::Real\{\s*(' + NUMBER + ')', body)
        if weight:
            record['molecular weight [kg mol-1]'] = float(weight[1])
        if diffusion:
            record['diffusion coefficient [m2 s-1]'] = float(diffusion)
        if 'species.SetThirdBody();' in body:
            record['parameterization'] = 'third body'
        species.append(record)
    assert len(species) == 210, f'expected 210 species, got {len(species)}'
    reactions = []
    for block in re.findall(r'processes.push_back\((.*?)\.Build\(\)\);', text, re.S):
        react = re.search(r'\.SetReactants\(\s*\{(.*?)\}\)', block, re.S)
        prod = re.search(r'\.SetProducts\(\s*\{(.*?)\}\)', block, re.S)
        rate = re.search(r'\.SetRateConstant\(\s*micm::(\w+)RateConstantParameters\{(.*?)\}\)', block, re.S)
        assert react and prod and rate, block
        kinds = {'Arrhenius': 'ARRHENIUS', 'UserDefined': 'USER_DEFINED', 'Surface': 'SURFACE', 'Troe': 'TROE'}
        rtype = kinds[rate[1]]
        rs = re.findall(r'SpeciesByName\(gas_phase, "([^"]+)"\)', react[1])
        ps = [{'species': name, 'yield': float(value)} for name, value in re.findall(
            r'StoichSpecies\{ SpeciesByName\(gas_phase, "([^"]+)"\),\s*(' + NUMBER + r')\s*\}', prod[1])]
        record = {'type': rtype, 'reactant sequence': rs, 'product sequence': ps}
        for key, value in re.findall(r'\.(\w+)_\s*=\s*(' + NUMBER + r')', rate[2]):
            record[{'scaling_factor': 'scaling factor', 'reaction_probability': 'reaction probability'}.get(key, key)] = float(value)
        label = re.search(r'\.label_\s*=\s*"([^"]+)"', rate[2])
        if label:
            record['MUSICA name'] = label[1]
        if rtype == 'SURFACE':
            sp = re.search(r'PhaseSpeciesByName\(gas_phase, "([^"]+)"\)', rate[2])
            assert sp and rs == [sp[1]]
            record['gas-phase reactant'] = sp[1]
        reactions.append(record)
    assert len(reactions) == 547, f'expected 547 reactions, got {len(reactions)}'
    initial = {
        'temperature': float(re.search(r'temperature_ = (' + NUMBER + r');', text)[1]),
        'pressure': float(re.search(r'pressure_ = (' + NUMBER + r');', text)[1]),
        'concentrations': {name: float(v) for name, v in re.findall(r'concentration\("([^"]+)", (' + NUMBER + r')\);', text)},
        'rate parameters': {name: float(v) for name, v in re.findall(r'parameter\("([^"]+)", (' + NUMBER + r')\);', text)},
    }
    return {'concentration units':'mol m-3','camp-data': species + [{'name': 'TS1', 'type': 'MECHANISM', 'reactions': reactions}], 'initial state': initial}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source')
    parser.add_argument('output')
    args = parser.parse_args()
    Path(args.output).write_text(json.dumps(translate(args.source), indent=2) + '\n')
