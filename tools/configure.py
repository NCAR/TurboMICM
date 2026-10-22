"""Strict build-time reader for CAMP-style JSON/YAML mechanism files.

The normalized representation is independent of the C emitter. Dict iteration
preserves input order; ordered term arrays additionally preserve duplicate terms.
"""
from copy import deepcopy
import json
import math
from pathlib import Path
import yaml

MOLECULES_PER_MOLE_M3 = 1e-6 * 6.02214076e23


def _document(path):
    with Path(path).open() as f:
        value = json.load(f) if Path(path).suffix.lower()=='.json' else yaml.safe_load(f)
    if not isinstance(value,dict):
        raise ValueError(f'{path}: expected a mapping')
    return value


def read(path):
    """Read a configuration, recursively resolving camp-files relative to its file."""
    entries, options, active = [], {}, set()
    def visit(path):
        path=Path(path).resolve()
        if path in active:
            raise ValueError(f'Configuration include cycle at {path.name}')
        active.add(path)
        doc=_document(path)
        for key in doc:
            if key not in ('camp-files','camp-data','initial conditions','initial rate parameters'):
                raise ValueError(f'{path.name}: unknown document key {key}')
        for name in doc.get('camp-files',[]):
            visit(path.parent/name)
        entries.extend(doc.get('camp-data',[]))
        for key in ('initial conditions','initial rate parameters'):
            options.setdefault(key,{}).update(doc.get(key,{}))
        active.remove(path)
    visit(path)
    species=[deepcopy(e) for e in entries if e['type']=='CHEM_SPEC']
    names={s['name']:i for i,s in enumerate(species)}
    if len(names)!=len(species) or not species:
        raise ValueError('Species must be nonempty and have unique names')
    reactions=[]
    def terms(value, products=False):
        if isinstance(value,dict):
            terms=[]
            for name,info in value.items():
                if products:
                    terms.append((name,float(info.get('yield',1))))
                else:
                    qty=info.get('qty',1)
                    if not isinstance(qty,int) or qty<1:
                        raise ValueError('Reactant qty must be a positive integer')
                    terms.extend([name]*qty)
            return terms
        if isinstance(value,list):
            return [(e['species'],float(e.get('yield',1))) for e in value] if products else list(value)
        raise ValueError('Reaction terms must be a mapping or ordered array')
    for entry in entries:
        if entry['type']=='CHEM_SPEC':
            continue
        if entry['type']!='MECHANISM':
            raise ValueError('Unknown camp-data type: '+entry['type'])
        for raw in entry['reactions']:
            raw=deepcopy(raw)
            kind=raw['type']
            if kind=='SURFACE':
                raw.setdefault('reactants',{raw['gas-phase reactant']:{}})
                raw.setdefault('products',raw.get('gas-phase products',{}))
            elif kind=='EMISSION':
                raw['reactants']={}; raw['products']={raw['species']:{}}
            elif kind=='FIRST_ORDER_LOSS':
                raw['reactants']={raw['species']:{}}; raw['products']={}
            if kind=='BRANCHED':
                branches=[('alkoxy','alkoxy products'),('nitrate','nitrate products')]
            else:
                branches=[(None,'products')]
            for branch,product_key in branches:
                rx=deepcopy(raw)
                rx['reactants']=terms(raw.get('reactants',{}))
                rx['products']=terms(raw.get(product_key,{}),True)
                rx['branch']=branch
                for name in rx['reactants']+[p[0] for p in rx['products']]:
                    if name not in names:
                        raise ValueError('Unknown species: '+name)
                rx['reactants']=[names[n] for n in rx['reactants']]
                rx['products']=[(names[n],v) for n,v in rx['products']]
                units=rx.get('concentration units','molecule cm-3')
                if units not in ('mol m-3','molecule cm-3'):
                    raise ValueError('Unknown concentration units: '+units)
                order=len(rx['reactants'])
                if units=='molecule cm-3':
                    factor=MOLECULES_PER_MOLE_M3
                    if kind in ('ARRHENIUS','TUNNELING','TAYLOR_SERIES','REVERSIBLE'):
                        rx['A']=float(rx.get('A',1))*factor**(order-1)
                    elif kind=='TROE':
                        rx['k0_A']=float(rx.get('k0_A',1))*factor**order
                        rx['kinf_A']=float(rx.get('kinf_A',1))*factor**(order-1)
                    elif kind=='TERNARY_CHEMICAL_ACTIVATION':
                        rx['k0_A']=float(rx.get('k0_A',1))*factor**(order-1)
                        rx['kinf_A']=float(rx.get('kinf_A',1))*factor**(order-2)
                    elif kind=='TERNARY_CHEMICAL_ACTIVATION_JPL19':
                        rx['k0_A']=float(rx.get('k0_A',1))*factor**order
                        rx['kinf_A']=float(rx.get('kinf_A',1))*factor**(order-1)
                        rx['kint_A']=float(rx.get('kint_A',1))*factor**(order-1)
                reactions.append(rx)
    return {'species':species,'reactions':reactions,**options}


def number(value):
    value=float(value)
    if not math.isfinite(value):
        raise ValueError('Configuration numbers must be finite')
    return repr(value)


def rate_initializer(rx, species):
    kind=rx['type']
    def fields(defaults):
        return ', '.join('.'+k+'='+number(rx.get(k,v)) for k,v in defaults.items())
    label=''
    if kind in ('ARRHENIUS','TAYLOR_SERIES'):
        p=fields(dict(A=1,B=0,C=0,D=300,E=0))
        if kind=='ARRHENIUS': typ, member='ARRHENIUS','arrhenius'
        else:
            coeff=rx.get('coefficients',[1])
            if not 1<=len(coeff)<=16: raise ValueError('Taylor series requires 1 to 16 coefficients')
            p=f'.arrhenius={{{p}}}, .n={len(coeff)}, .coefficients={{{",".join(map(number,coeff))}}}'
            typ,member='TAYLOR','taylor'
    elif kind in ('TROE','TERNARY_CHEMICAL_ACTIVATION'):
        typ='TROE' if kind=='TROE' else 'TERNARY'; member='falloff'
        p=fields(dict(k0_A=1,k0_B=0,k0_C=0,kinf_A=1,kinf_B=0,kinf_C=0,Fc=.6,N=1))
    elif kind=='TERNARY_CHEMICAL_ACTIVATION_JPL19':
        typ,member='TERNARY_JPL19','jpl19'
        p=fields(dict(k0_A=1,k0_B=0,k0_C=0,k0_D=298,kinf_A=1,kinf_B=0,kinf_C=0,kinf_D=298,
                      kint_A=1,kint_B=0,kint_C=0,kint_D=298,Fc=.6,N=1))
    elif kind=='TUNNELING':
        typ,member='TUNNELING','tunneling'; p=fields(dict(A=1,B=0,C=0))
    elif kind=='BRANCHED':
        typ,member='BRANCHED','branched'
        p=fields(dict(X=1,Y=0,a0=0))+', .n=(size_t)'+number(rx.get('n',0))+f", .nitrate={int(rx['branch']=='nitrate')}"
    elif kind=='REVERSIBLE':
        typ,member='REVERSIBLE','reversible'; p=fields(dict(A=1,C=0,k_r=0))
    elif kind in ('USER_DEFINED','PHOTOLYSIS','EMISSION','FIRST_ORDER_LOSS'):
        typ,member='USER','user'; p='.scaling='+number(rx.get('scaling factor',1))
        label=rx['MUSICA name']
    elif kind=='SURFACE':
        typ,member='SURFACE','surface'; label=rx['MUSICA name']
        sp=species[rx['reactants'][0]]
        p=('.diffusion='+number(sp['diffusion coefficient [m2 s-1]'])+
           ', .molecular_weight='+number(sp['molecular weight [kg mol-1]'])+
           ', .probability='+number(rx.get('reaction probability',1)))
    else:
        raise ValueError('Unsupported rate type: '+kind)
    return '{.type=CHEM_'+typ+', .label='+json.dumps(label)+', .p.'+member+'={'+p+'}}'
