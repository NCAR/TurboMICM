"""Read MICM YAML/JSON mechanisms and emit C model constructors.

The reader normalizes file inclusion, rate defaults, stoichiometry and units.
Code generation is deliberately separate from parsing, so tests can also build
the normalized mechanism directly through the public C API.
"""
from dataclasses import dataclass, field
import argparse
import copy
import json
import math
from pathlib import Path
import re
import yaml

AVOGADRO_CM = 6.02214076e23 * 1e-6
RATE_NAMES = ['ARRHENIUS', 'TROE', 'TERNARY_CHEMICAL_ACTIVATION',
              'TERNARY_CHEMICAL_ACTIVATION_JPL19', 'TUNNELING', 'BRANCHED',
              'SURFACE', 'USER_DEFINED', 'TAYLOR_SERIES', 'REVERSIBLE']
FIELDS = {
    'ARRHENIUS': [('A',1),('B',0),('C',0),('D',300),('E',0)],
    'TROE': [('k0_A',1),('k0_B',0),('k0_C',0),('kinf_A',1),('kinf_B',0),('kinf_C',0),('Fc',.6),('N',1)],
    'TUNNELING': [('A',1),('B',0),('C',0)],
    'BRANCHED': [('X',0),('Y',0),('a0',0),('n',0),('branch',0)],
    'SURFACE': [('diffusion',0),('speed factor',0),('reaction probability',1)],
    'USER_DEFINED': [('scaling factor',1)],
    'REVERSIBLE': [('A',1),('C',0),('k_r',0)],
}
FIELDS['TERNARY_CHEMICAL_ACTIVATION'] = FIELDS['TROE']
FIELDS['TERNARY_CHEMICAL_ACTIVATION_JPL19'] = [
    (prefix+key, default) for prefix in ('k0_', 'kinf_', 'kint_')
    for key,default in [('A',1),('B',0),('C',0),('D',298)]] + [('Fc',.6),('N',1)]
FIELDS['TAYLOR_SERIES'] = FIELDS['ARRHENIUS']


def number(value):
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f'non-finite numeric value: {value}')
    return result


def terms(value, reactants=False):
    if isinstance(value, list):
        if reactants:
            return list(value)
        return [(p['species'],number(p.get('yield',1))) for p in value]
    if not isinstance(value,dict):
        raise ValueError('reaction terms must be an object or an ordered list')
    if reactants:
        result = []
        for name, props in value.items():
            count = props.get('qty',1)
            if int(count) != count or count < 1:
                raise ValueError('reactant qty must be a positive integer')
            result.extend([name]*int(count))
        return result
    return [(name,number(props.get('yield',1))) for name,props in value.items()]


@dataclass
class Mechanism:
    species: list = field(default_factory=list)
    reactions: list = field(default_factory=list)
    parameters: list = field(default_factory=list)
    initial: dict = field(default_factory=dict)


def read_config(path):
    data, initial, stack = [], {}, set()
    def visit(path):
        path = path.resolve()
        if path in stack:
            raise ValueError(f'cyclic camp-files inclusion: {path.name}')
        stack.add(path)
        document = json.loads(path.read_text()) if path.suffix.lower() == '.json' else yaml.safe_load(path.read_text())
        if not isinstance(document,dict):
            raise ValueError('configuration document must be a mapping')
        for child in document.get('camp-files',[]):
            visit(path.parent/child)
        data.extend(document.get('camp-data',[]))
        initial.update(document.get('initial state',{}))
        stack.remove(path)
    visit(Path(path))
    model = Mechanism(initial=initial)
    for item in data:
        if item['type'] == 'CHEM_SPEC':
            name = item['name']
            if any(s['name'] == name for s in model.species):
                raise ValueError(f'duplicate species {name}')
            if item.get('parameterization') not in (None, 'air density'):
                raise ValueError('unknown species parameterization')
            model.species.append(copy.deepcopy(item))
        elif item['type'] != 'MECHANISM':
            raise ValueError(f'unknown camp-data type {item["type"]}')
    known = {s['name']: s for s in model.species}
    def parameter(name):
        if name not in model.parameters:
            model.parameters.append(name)
        return model.parameters.index(name)
    for mechanism in data:
        if mechanism['type'] != 'MECHANISM':
            continue
        units = mechanism.get('rate units','molecule cm-3')
        if units not in ('molecule cm-3','mol m-3'):
            raise ValueError(f'unknown rate units {units}')
        for original in mechanism.get('reactions',[]):
            parts = [copy.deepcopy(original)]
            if original['type'] == 'BRANCHED' and 'branch' not in original:
                parts = []
                for i,key in enumerate(('alkoxy products','nitrate products')):
                    r = copy.deepcopy(original)
                    r['branch'] = i
                    r['products'] = r[key]
                    parts.append(r)
            for r in parts:
                kind = r['type']
                if kind == 'SURFACE':
                    r['reactants'] = [r['gas-phase reactant']]
                    r['products'] = r['gas-phase products']
                if kind in ('PHOTOLYSIS','FIRST_ORDER_LOSS','EMISSION'):
                    if kind == 'FIRST_ORDER_LOSS':
                        r['reactants'],r['products'] = [r['species']],[]
                    if kind == 'EMISSION':
                        r['reactants'],r['products'] = [],[{'species':r['species']}]
                    kind = 'USER_DEFINED'
                if kind not in FIELDS:
                    raise ValueError(f'unsupported rate type {kind}')
                reactants, products = terms(r.get('reactants',{}),True), terms(r.get('products',{}))
                for name in reactants + [name for name,_ in products]:
                    if name not in known:
                        raise ValueError(f'unknown species {name}')
                values = [number(r.get(k,default)) for k,default in FIELDS[kind]]
                coefficients = 0
                if kind == 'TAYLOR_SERIES':
                    coeffs = r.get('coefficients',[1])
                    if not 1 <= len(coeffs) <= 16:
                        raise ValueError('Taylor series needs 1 to 16 coefficients')
                    coefficients = len(coeffs)
                    values.extend(map(number,coeffs))
                if units == 'molecule cm-3':
                    order = len(reactants)
                    scale = AVOGADRO_CM**(order-1)
                    if kind in ('ARRHENIUS','TUNNELING','BRANCHED','TAYLOR_SERIES','REVERSIBLE'):
                        values[0] *= scale
                    elif kind == 'TROE':
                        values[0] *= AVOGADRO_CM**order
                        values[3] *= scale
                    elif kind == 'TERNARY_CHEMICAL_ACTIVATION':
                        values[0] *= scale
                        values[3] *= AVOGADRO_CM**(order-2)
                    elif kind == 'TERNARY_CHEMICAL_ACTIVATION_JPL19':
                        values[0] *= AVOGADRO_CM**order
                        values[4] *= scale
                        values[8] *= scale
                index = -1
                if kind == 'USER_DEFINED':
                    index = parameter(r['MUSICA name'])
                elif kind == 'SURFACE':
                    s = known[reactants[0]]
                    for key in ('molecular weight [kg mol-1]','diffusion coefficient [m2 s-1]'):
                        if key not in s:
                            raise ValueError(f'surface species {s["name"]} lacks {key}')
                    label = r['MUSICA name']
                    index = parameter(label+'.effective radius [m]')
                    if parameter(label+'.particle number concentration [# m-3]') != index+1:
                        raise ValueError('surface parameters must be adjacent')
                model.reactions.append(dict(type=RATE_NAMES.index(kind), reactants=reactants,
                    products=products, p=values, custom_index=index, coefficient_count=coefficients))
    if not model.species:
        raise ValueError('mechanism has no species')
    for key,allowed in [('concentrations',known),('parameters',model.parameters)]:
        for name,value in initial.get(key,{}).items():
            if name not in allowed:
                raise ValueError(f'unknown initial {key} name: {name}')
            number(value)
    return model


def c_array(type_, values):
    return f'({type_}[]){{' + ','.join(map(str,values)) + '}' if values else 'NULL'


def generate(model, name):
    if not re.fullmatch(r'[a-zA-Z_][a-zA-Z0-9_]*',name):
        raise ValueError('invalid C identifier')
    q = json.dumps
    names = {s['name']:i for i,s in enumerate(model.species)}
    lines = ['/* Generated by tools/configure.py. */', '#include "chem.h"',
             f'ch_model *ch_build_{name}(void) {{', '  ch_model *m = ch_model_create();',
             '  if (!m) return NULL;']
    for s in model.species:
        lines.append(f'  if (ch_model_add_species(m,{q(s["name"])},{number(s.get("absolute tolerance",1e-3))},{int(s.get("parameterization") == "air density")}) < 0) goto fail;')
        for key in ('molecular weight [kg mol-1]','diffusion coefficient [m2 s-1]'):
            if key in s:
                lines.append(f'  if (ch_model_set_species_property(m,{names[s["name"]]},{q(key)},{number(s[key])})) goto fail;')
    for label in model.parameters:
        initial = number(model.initial.get('parameters',{}).get(label,0))
        lines.append(f'  if (ch_model_add_parameter(m,{q(label)},{initial}) < 0) goto fail;')
    for r in model.reactions:
        ri = c_array('int',[names[x] for x in r['reactants']])
        pi = c_array('int',[names[x] for x,_ in r['products']])
        yi = c_array('double',[y for _,y in r['products']])
        rate = '{.type='+str(r['type'])+',.p={'+','.join(map(str,r['p']))+'},.custom_index='+str(r['custom_index'])+',.coefficient_count='+str(r['coefficient_count'])+'}'
        lines.append(f'  if (ch_model_add_reaction(m,{len(r["reactants"])},{ri},{len(r["products"])},{pi},{yi},&(ch_rate){rate}) < 0) goto fail;')
    lines.extend(['  if (ch_model_finalize(m,1,1)) goto fail;', '  return m;', 'fail:', '  ch_model_free(m); return NULL;', '}',
                  f'int ch_init_{name}(ch_state *s) {{', '  if (!s) return CH_INVALID;'])
    for label,value in model.initial.get('concentrations',{}).items():
        lines.append(f'  if (ch_state_set(s,{q(label)},(double[]){{{number(value)}}},1)) return CH_INVALID;')
    for label,value in model.initial.get('parameters',{}).items():
        lines.append(f'  if (ch_state_set_parameter(s,{q(label)},(double[]){{{number(value)}}},1)) return CH_INVALID;')
    t,p = number(model.initial.get('temperature',0)),number(model.initial.get('pressure',0))
    lines.extend(['  ch_conditions *c = ch_state_conditions(s);',
                  '  for (size_t i = 0; i < ch_state_cells(s); ++i) {',
                  f'    c[i].temperature = {t}; c[i].pressure = {p};',
                  '    ch_ideal_air_density(&c[i]);', '  }', '  return CH_OK;', '}'])
    return '\n'.join(lines)+'\n'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--name',default='mechanism')
    args = parser.parse_args()
    model = read_config(args.config)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(generate(model,args.name))
