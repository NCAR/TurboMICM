"""Read the CAMP-style mechanism files used by the MICM tutorials.

The result retains the input ordering. It is independent of the code generator so
it can also be used by validation tools. Invalid references fail at build time.
"""
import json
from pathlib import Path

import yaml


class MechanismError(ValueError):
    pass


def read_document(path):
    path = Path(path)
    with path.open(encoding='utf-8') as stream:
        if path.suffix.lower() == '.json':
            data = json.load(stream)
        elif path.suffix.lower() in ('.yaml', '.yml'):
            data = yaml.safe_load(stream)
        else:
            raise MechanismError(f'{path}: expected JSON or YAML')
    if not isinstance(data, dict):
        raise MechanismError(f'{path}: expected a mapping')
    return data


def read_mechanism(path):
    path = Path(path)
    config = read_document(path)
    files = config.get('camp-files')
    if not isinstance(files, list) or not files:
        raise MechanismError('camp-files must be a nonempty list')
    species, reactions = [], []
    for filename in files:
        doc = read_document(path.parent / filename)
        for item in doc.get('camp-data', []):
            kind = item.get('type')
            if kind == 'CHEM_SPEC':
                species.append(item)
            elif kind == 'MECHANISM':
                reactions.extend(item.get('reactions', []))
            else:
                raise MechanismError(f'{filename}: unsupported camp-data type {kind!r}')
    names = [item.get('name') for item in species]
    if not names or any(not isinstance(n, str) or not n for n in names) or len(set(names)) != len(names):
        raise MechanismError('species names must be nonempty and unique')
    allowed = {'ARRHENIUS', 'BRANCHED', 'SURFACE', 'TERNARY_CHEMICAL_ACTIVATION', 'TROE', 'TUNNELING',
               'PHOTOLYSIS', 'USER_DEFINED', 'FIRST_ORDER_LOSS', 'EMISSION'}
    for i, reaction in enumerate(reactions):
        kind = reaction.get('type')
        if kind not in allowed:
            raise MechanismError(f'reaction {i}: unsupported type {kind!r}')
        for key in ('reactants', 'products', 'alkoxy products', 'nitrate products', 'gas-phase products'):
            terms = reaction.get(key, {})
            if not isinstance(terms, dict):
                raise MechanismError(f'reaction {i}: {key} must be a mapping')
            for name, props in terms.items():
                if name not in names or not isinstance(props, dict):
                    raise MechanismError(f'reaction {i}: invalid {key} species {name!r}')
                qty = props.get('qty', 1)
                if key == 'reactants' and (not isinstance(qty, int) or isinstance(qty, bool) or qty < 1):
                    raise MechanismError(f'reaction {i}: reactant qty must be a positive integer')
        for key in ('species', 'gas-phase reactant'):
            if key in reaction and reaction[key] not in names:
                raise MechanismError(f'reaction {i}: unknown {key} {reaction[key]!r}')
    return {'species': species, 'reactions': reactions, 'initial conditions': config.get('initial conditions', {})}
