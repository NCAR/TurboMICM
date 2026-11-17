"""Translate the checked-in TS1 C++ mechanism into CAMP-style JSON input.

The source is generated C++ with a deliberately regular shape. Refuse to emit
partial output if any expected declaration cannot be translated.
"""
import argparse
import collections
import json
from pathlib import Path
import re

NUMBER = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
SPECIES = r'SpeciesByName\(gas_phase,\s*"([^"]+)"\)'


def one(pattern, text, what):
    matches = re.findall(pattern, text, re.S)
    if len(matches) != 1:
        raise ValueError(f'expected one {what}, found {len(matches)}')
    return matches[0]


def translate(source):
    phase, remainder = source.split('inline std::vector<micm::Process> CreateProcesses', 1)
    process_text, initialization = remainder.split('static void InitState', 1)
    concentrations = dict((k, float(v)) for k, v in re.findall(r'concentration\("([^"]+)",\s*(' + NUMBER + r')\)', initialization))
    defaults = dict((k, float(v)) for k, v in re.findall(r'parameter\("([^"]+)",\s*(' + NUMBER + r')\)', initialization))
    conditions = {
        'temperature [K]': float(one(r'conditions_\[c\]\.temperature_\s*=\s*(' + NUMBER + r')\s*;', initialization, 'default temperature')),
        'pressure [Pa]': float(one(r'conditions_\[c\]\.pressure_\s*=\s*(' + NUMBER + r')\s*;', initialization, 'default pressure')),
    }
    species = []
    blocks = re.findall(r'auto species = micm::Species\("([^"]+)"\);(.*?)phase_species\.emplace_back\(species(?:,\s*(' + NUMBER + r'))?\);', phase, re.S)
    for name, body, diffusion in blocks:
        entry = {'name': name, 'type': 'CHEM_SPEC'}
        if 'SetThirdBody()' in body:
            entry['third body'] = True
        else:
            entry['molecular weight [kg mol-1]'] = float(one(r'MOLECULAR_WEIGHT,\s*micm::Real\{\s*(' + NUMBER + r')\s*\}', body, 'molecular weight'))
        if diffusion:
            entry['diffusion coefficient [m2 s-1]'] = float(diffusion)
        if name in concentrations:
            entry['initial concentration [mol m-3]'] = concentrations.pop(name)
        species.append(entry)
    if len(species) != 210 or concentrations:
        raise ValueError(f'species mismatch: {len(species)}; missing {list(concentrations)}')
    reactions = []
    chunks = process_text.split('processes.push_back(')[1:]
    for i, chunk in enumerate(chunks):
        chunk = chunk.split('.Build());', 1)[0]
        reactant_text = one(r'\.SetReactants\(\s*\{(.*?)\}\s*\)', chunk, 'reactants')
        product_text = one(r'\.SetProducts\(\s*\{(.*?)\}\s*\)', chunk, 'products')
        reactants = collections.Counter(re.findall(SPECIES, reactant_text))
        products = re.findall(r'micm::StoichSpecies\{\s*' + SPECIES + r'\s*,\s*(' + NUMBER + r')\s*\}', product_text)
        if len(products) != len(re.findall('StoichSpecies', product_text)) or len(products) != len(dict(products)):
            raise ValueError(f'cannot parse products of reaction {i}')
        rate_type, body = one(r'\.SetRateConstant\(\s*micm::(\w+)\s*\{(.*?)\}\s*\)', chunk, 'rate')
        fields = dict(re.findall(r'\.(\w+)_\s*=\s*(' + NUMBER + r'|"[^"]+"|PhaseSpeciesByName\(gas_phase,\s*"[^"]+"\))', body))
        r = {'reactants': {n: {'qty': q} for n, q in reactants.items()},
             'products': {n: {'yield': float(y)} for n, y in products},
             'rate coefficient units': 'SI'}
        if rate_type in ('ArrheniusRateConstantParameters', 'TroeRateConstantParameters'):
            r['type'] = 'ARRHENIUS' if rate_type.startswith('Arrhenius') else 'TROE'
            expected = {'A','B','C','D','E'} if r['type'] == 'ARRHENIUS' else {'k0_A','k0_B','k0_C','kinf_A','kinf_B','kinf_C','Fc','N'}
            if fields.keys() != expected:
                raise ValueError(f'reaction {i} {rate_type} fields: {list(fields)}')
            r.update({key: float(value) for key, value in fields.items()})
        elif rate_type == 'UserDefinedRateConstantParameters':
            if fields.keys() != {'label', 'scaling_factor'}:
                raise ValueError(f'unexpected user rate fields in reaction {i}')
            r.update({'type': 'USER_DEFINED', 'MUSICA name': fields['label'].strip('"'), 'scaling factor': float(fields['scaling_factor'])})
        elif rate_type == 'SurfaceRateConstantParameters':
            if fields.keys() != {'label','phase_species','reaction_probability'}:
                raise ValueError(f'unexpected surface fields in reaction {i}')
            name = one(r'PhaseSpeciesByName\(gas_phase,\s*"([^"]+)"\)', fields['phase_species'], 'surface species')
            r.update({'type': 'SURFACE', 'gas-phase reactant': name, 'gas-phase products': r.pop('products'),
                      'surface parameter label': fields['label'].strip('"'), 'reaction probability': float(fields['reaction_probability'])})
            if reactants != {name: 1}:
                raise ValueError(f'unexpected surface reactants in reaction {i}')
            r.pop('reactants')
        else:
            raise ValueError(f'unsupported source rate: {rate_type}')
        reactions.append(r)
    if len(reactions) != 547:
        raise ValueError(f'expected 547 reactions, got {len(reactions)}')
    return species, reactions, defaults, conditions


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('source', type=Path)
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    species, reactions, defaults, conditions = translate(args.source.read_text(encoding='utf-8'))
    args.directory.mkdir(parents=True, exist_ok=True)
    files = {'config.json': {'camp-files': ['species.json', 'reactions.json'], 'custom parameter defaults': defaults,
                             **conditions},
             'species.json': {'camp-data': species},
             'reactions.json': {'camp-data': [{'name': 'TS1', 'type': 'MECHANISM', 'reactions': reactions}]}}
    for name, content in files.items():
        (args.directory/name).write_text(json.dumps(content, indent=2) + '\n', encoding='utf-8')

if __name__ == '__main__':
    main()
