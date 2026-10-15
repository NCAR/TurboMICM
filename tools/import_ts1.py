"""Translate the read-only TS1 C++ definition into an ordered CAMP configuration.

No C++ is compiled or executed. The generated JSON is the benchmark's build input;
normal library/benchmark execution is independent of the reference tree.
"""
import argparse
import json
import pathlib
import re

NUMBER = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'


def translate(source):
    species = []
    phase = source.split('inline std::vector<micm::Process> CreateProcesses')[0]
    for block in re.findall(r'auto species = micm::Species\("([^"]+)"\);(.*?)phase_species.emplace_back\(species(.*?)\);', phase, re.S):
        name, properties, diffusion = block
        item = {'name': name, 'type': 'CHEM_SPEC'}
        weight = re.search(r'MOLECULAR_WEIGHT, micm::Real\{\s*(' + NUMBER + r')\s*\}', properties)
        if weight:
            item['molecular weight [kg mol-1]'] = float(weight[1])
        if diffusion.strip():
            item['diffusion coefficient [m2 s-1]'] = float(diffusion.strip().removeprefix(',').strip())
        if 'SetThirdBody' in properties:
            item['parameterization'] = 'air-density'
        species.append(item)
    reactions = []
    for block in source.split('processes.push_back(')[1:]:
        block = block.split('.Build());', 1)[0]
        reactant_block = re.search(r'\.SetReactants\(\s*\{(.*?)\}\)', block, re.S)
        product_block = re.search(r'\.SetProducts\(\s*\{(.*?)\}\)', block, re.S)
        rate = re.search(r'\.SetRateConstant\(\s*micm::(\w+)RateConstantParameters\{(.*?)\}\)', block, re.S)
        if not rate:
            raise ValueError('Unrecognized reaction: ' + block[:120])
        item = {'type': {'UserDefined': 'USER_DEFINED'}.get(rate[1], rate[1].upper())}
        item['reactants'] = re.findall(r'SpeciesByName\(gas_phase, "([^"]+)"\)', reactant_block[1]) if reactant_block else []
        item['products'] = [{'species': name, 'yield': float(value)} for name, value in
                            re.findall(r'StoichSpecies\{\s*SpeciesByName\(gas_phase, "([^"]+)"\),\s*(' + NUMBER + r')\s*\}', product_block[1])] if product_block else []
        for key, value in re.findall(r'\.(\w+)_\s*=\s*(' + NUMBER + r')', rate[2]):
            item[key] = float(value)
        label = re.search(r'\.label_ = "([^"]+)"', rate[2])
        if label:
            item['MUSICA name'] = label[1]
        if rate[1] == 'Surface':
            item['gas-phase reactant'] = re.search(r'PhaseSpeciesByName\(gas_phase, "([^"]+)"\)', rate[2])[1]
            item['reaction probability'] = item.pop('reaction_probability')
        if not item['reactants']:
            raise ValueError('Missing reactants in reaction ' + str(len(reactions)))
        reactions.append(item)
    concentrations = {name: float(v) for name, v in re.findall(r'concentration\("([^"]+)",\s*(' + NUMBER + r')\);', source)}
    custom = {name: float(v) for name, v in re.findall(r'parameter\("([^"]+)",\s*(' + NUMBER + r')\);', source)}
    assert len(species) == 210, len(species)
    assert len(reactions) == 547, len(reactions)
    assert len(concentrations) == 209, len(concentrations)
    return {'concentration units': 'mol m-3',
            'initial state': {'temperature': 287.45, 'pressure': 101319.9, 'concentrations': concentrations, 'custom': custom},
            'camp-data': species + [{'name': 'TS1', 'type': 'MECHANISM', 'reactions': reactions}]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=pathlib.Path)
    parser.add_argument('output', type=pathlib.Path)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(translate(args.source.read_text()), indent=2) + '\n')


if __name__ == '__main__':
    main()
