"""Translate the specified C++ mechanism data into extended CAMP JSON.

This reads data declarations only. It never compiles or executes the input.
Unknown declarations fail rather than silently dropping reactions.
"""
import argparse
import json
from pathlib import Path
import re

NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"


def translate(path):
    text = Path(path).read_text()
    concentrations = dict((n, float(v)) for n, v in re.findall(
        rf'concentration\("([^"]+)",\s*({NUM})\)', text))
    parameters = dict((n, float(v)) for n, v in re.findall(
        rf'parameter\("([^"]+)",\s*({NUM})\)', text))
    species = []
    for name, body in re.findall(r'auto species = micm::Species\("([^"]+)"\);(.*?)\n      }', text, re.S):
        s = {'name': name, 'type': 'CHEM_SPEC', 'initial concentration': concentrations.get(name, 0)}
        mw = re.search(rf'MOLECULAR_WEIGHT, micm::Real{{\s*({NUM})\s*}}', body)
        if mw:
            s['molecular weight [kg mol-1]'] = float(mw[1])
        diff = re.search(rf'emplace_back\(species,\s*({NUM})\)', body)
        if diff:
            s['diffusion coefficient [m2 s-1]'] = float(diff[1])
        if 'SetThirdBody()' in body:
            s['third body'] = True
        species.append(s)
    reactions = []
    types = {'Arrhenius': 'ARRHENIUS', 'Troe': 'TROE', 'UserDefined': 'USER_DEFINED', 'Surface': 'SURFACE'}
    for block in text.split('processes.push_back(')[1:]:
        block = block.split('.Build());', 1)[0]
        reactants_text = re.search(r'\.SetReactants\(\s*{(.*?)}\)', block, re.S)[1]
        reactants = re.findall(r'SpeciesByName\(gas_phase, "([^"]+)"\)', reactants_text)
        # Lists preserve repeated reactants and the source's multiplication order.
        products_text = block.split('.SetProducts({', 1)[1].split('.SetRateConstant', 1)[0]
        products = [{'species': n, 'yield': float(v)} for n, v in re.findall(
            rf'StoichSpecies{{\s*SpeciesByName\(gas_phase, "([^"]+)"\),\s*({NUM})\s*}}', products_text)]
        rate_type, rate_text = re.search(r'micm::(\w+)RateConstantParameters\{(.*?)\}\)', block, re.S).groups()
        r = {'type': types[rate_type], 'reactants': reactants, 'products': products}
        for key, val in re.findall(rf'\.(\w+)_\s*=\s*({NUM})', rate_text):
            r[key.replace('reaction_probability', 'reaction probability')] = float(val)
        label = re.search(r'\.label_\s*=\s*"([^"]+)"', rate_text)
        if label:
            r['MUSICA name'] = label[1]
        if rate_type == 'Surface':
            r['gas-phase reactant'] = reactants[0]
            r['gas-phase products'] = products
            del r['reactants'], r['products']
        reactions.append(r)
    if len(species) != 210 or len(reactions) != 547:
        raise ValueError(f'unexpected TS1 counts: {len(species)} species, {len(reactions)} reactions')
    return {'camp-data': species + [{'name': 'TS1', 'type': 'MECHANISM', 'reactions': reactions,
        'initial conditions': {'temperature': 287.45, 'pressure': 101319.9},
        'initial rate parameters': parameters}]}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source')
    p.add_argument('output')
    a = p.parse_args()
    Path(a.output).write_text(json.dumps(translate(a.source), indent=2) + '\n')
