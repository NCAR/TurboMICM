"""Translate the supplied TS1 C++ data declarations into a mechanism configuration.

No C++ is compiled or executed. Unexpected declarations fail rather than being
silently skipped. The result is read by the ordinary configuration compiler.
"""
import argparse
import json
from pathlib import Path
import re


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    text = args.source.read_text()
    head, tail = text.split('inline std::vector<micm::Process> CreateProcesses', 1)
    species = []
    for name, body, diffusion in re.findall(r'auto species = micm::Species\("([^"]+)"\);(.*?)phase_species.emplace_back\(species(.*?)\);', head, re.S):
        item = {'name': name, 'type': 'CHEM_SPEC'}
        weight = re.search(r'MOLECULAR_WEIGHT, micm::Real\{\s*([^}]+)\}', body)
        if weight:
            item['molecular weight [kg mol-1]'] = float(weight[1])
        if diffusion.strip():
            item['diffusion coefficient [m2 s-1]'] = float(diffusion.strip().removeprefix(','))
        if 'species.SetThirdBody();' in body:
            item['parameterization'] = 'air density'
        species.append(item)
    reactions = []
    pattern = r'processes.push_back\((.*?)\.Build\(\)\);'
    kinds = {'Arrhenius': 'ARRHENIUS', 'UserDefined': 'USER_DEFINED', 'Surface': 'SURFACE', 'Troe': 'TROE'}
    for block in re.findall(pattern, tail, re.S):
        reactants = re.search(r'\.SetReactants\(\s*\{(.*?)\}\s*\)', block, re.S)[1]
        products = re.search(r'\.SetProducts\(\s*\{(.*?)\}\s*\)', block, re.S)[1]
        kind, data = re.search(r'\.SetRateConstant\(\s*micm::(\w+)RateConstantParameters\{(.*?)\}\)', block, re.S).groups()
        rs = re.findall(r'SpeciesByName\(gas_phase, "([^"]+)"\)', reactants)
        ps = [(name, float(y)) for name, y in re.findall(r'SpeciesByName\(gas_phase, "([^"]+)"\),\s*([^}]+)\}', products)]
        item = {'type': kinds[kind], 'reactants': rs,
                'products': [{'species': name, 'yield': y} for name, y in ps]}
        for key, value in re.findall(r'\.(\w+)_\s*=\s*([^,}]+)', data):
            if key == 'label':
                item['MUSICA name'] = json.loads(value.strip())
            elif key == 'phase_species':
                continue
            else:
                item[{'reaction_probability': 'reaction probability', 'scaling_factor': 'scaling factor'}.get(key,key)] = float(value)
        if kind == 'Surface':
            assert len(rs) == 1
            item['gas-phase reactant'] = rs[0]
            item['gas-phase products'] = item.pop('products')
            del item['reactants']
        reactions.append(item)
    assert len(species) == 210, len(species)
    assert len(reactions) == 547, len(reactions)
    initial = {'temperature': 287.45, 'pressure': 101319.9,
               'concentrations': {n: float(v) for n,v in re.findall(r'concentration\("([^"]+)",\s*([^)]+)\)', text)},
               'parameters': {n: float(v) for n,v in re.findall(r'parameter\("([^"]+)",\s*([^)]+)\)', text)}}
    result = {'camp-data': species + [{'name': 'TS1', 'type': 'MECHANISM',
                                     'rate units': 'mol m-3', 'reactions': reactions}],
              'initial state': initial}
    args.destination.parent.mkdir(parents=True, exist_ok=True)
    args.destination.write_text(json.dumps(result, indent=2) + '\n')
    print(f'Translated {len(species)} species and {len(reactions)} reactions')


if __name__ == '__main__':
    main()
