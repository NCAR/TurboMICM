"""Additional configuration covering every configurable rate and repeated terms.

This supplements the unchanged Tier 1 configuration mechanisms and TS1. Rates
are deliberately small enough for trajectories under both tolerance settings.
"""
import json
import pathlib

species = [{'type': 'CHEM_SPEC', 'name': s} for s in ['A', 'B', 'C', 'M', 'isolated']]
species[0].update({'diffusion coefficient [m2 s-1]': 1e-5,
                   'molecular weight [kg mol-1]': 0.05})
species[3]['parameterization'] = 'air-density'
kinds = ['ARRHENIUS', 'TROE', 'TERNARY_CHEMICAL_ACTIVATION',
         'TERNARY_CHEMICAL_ACTIVATION_JPL19', 'TUNNELING', 'BRANCHED',
         'TAYLOR_SERIES', 'REVERSIBLE', 'USER_DEFINED', 'SURFACE',
         'PHOTOLYSIS', 'FIRST_ORDER_LOSS', 'EMISSION']
reactions = []
for i, kind in enumerate(kinds):
    r = dict(type=kind, reactants=['A', 'B', 'A', 'M'],
             products=[{'species': 'C', 'yield': 0.73}, {'species': 'C', 'yield': 0.2},
                       {'species': 'A', 'yield': 0.1}],
             A=1e-6, B=0.12, C=3.7, D=290, E=1e-8,
             k0_A=1e-7, kinf_A=1e-5, kint_A=1e-8,
             k0_B=0.1, kinf_B=-0.4, kint_C=2.1,
             X=1e-6, Y=100, a0=0.1, n=2,
             coefficients=[0.1, 0.002, -1e-6], reverse=0.05,
             label='parameter_' + str(i))
    if kind == 'BRANCHED':
        r['alkoxy products'] = r['products']
        r['nitrate products'] = {'B': {}}
    if kind == 'SURFACE':
        r['gas-phase reactant'] = 'A'
        r['gas-phase products'] = {'B': {}}
    if kind in ('EMISSION', 'FIRST_ORDER_LOSS'):
        r['species'] = 'B'
    reactions.append(r)
data = {'concentration units': 'mol m-3', 'camp-data': species + [
    {'type': 'MECHANISM', 'name': 'all configurable rates', 'reactions': reactions}]}
path = pathlib.Path('generated/extended.json')
path.parent.mkdir(exist_ok=True)
path.write_text(json.dumps(data, indent=2) + '\n')
