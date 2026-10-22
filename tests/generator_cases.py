"""Small independent configurations covering every reader rate family and edge.

These are test inputs, never an alternate solver or expected-result generator.
The native test compares the emitted chemistry to the generic library.
"""
import json
from pathlib import Path
import sys


def fixtures(destination):
    species = [dict(type='CHEM_SPEC',name=n) for n in ('A','B','C','unused')]
    species[0].update({'molecular weight [kg mol-1]':.025, 'diffusion coefficient [m2 s-1]':2.3e-5})
    species += [dict(type='CHEM_SPEC',name='M',parameterization='air density')]
    rates = [
        dict(type='ARRHENIUS',A=.02,B=-.5,C=50,D=298,E=1e-8),
        dict(type='TROE',k0_A=1e-4,k0_B=-1,kinf_A=.1,Fc=.6),
        dict(type='TERNARY_CHEMICAL_ACTIVATION',k0_A=1e-4,kinf_A=.1),
        dict(type='TERNARY_CHEMICAL_ACTIVATION_JPL19',k0_A=1e-4,kinf_A=.1,kint_A=.01),
        dict(type='TUNNELING',A=.01,B=-12,C=100),
        dict(type='BRANCHED',X=.01,Y=100,a0=.2,n=1,branch=0),
        dict(type='BRANCHED',X=.01,Y=100,a0=.2,n=1,branch=1),
        dict(type='TAYLOR_SERIES',A=.01,B=.2,C=-100,coefficients=[.1,1e-4,1e-8]),
        dict(type='REVERSIBLE',A=.01,C=100,k_r=.2),
        {'type':'USER_DEFINED','MUSICA name':'custom','scaling factor':.5},
    ]
    reactions = [dict(rate,reactants=['A'],products=[{'species':'B','yield':.75}]) for rate in rates]
    reactions += [
        {'type':'SURFACE','gas-phase reactant':'A','gas-phase products':{'C':{}},'MUSICA name':'surface','reaction probability':.3},
        {'type':'ARRHENIUS','A':.02,'reactants':['A','B','A','M','M'],
         'products':[{'species':'C','yield':.3},{'species':'B','yield':1},{'species':'C','yield':.2},{'species':'M'}]},
        {'type':'EMISSION','species':'A','MUSICA name':'source'},
        {'type':'FIRST_ORDER_LOSS','species':'C','MUSICA name':'loss'},
        {'type':'PHOTOLYSIS','reactants':{'B':{}},'products':{'C':{}},'MUSICA name':'jtest'},
    ]
    document = {'camp-data':species+[{'type':'MECHANISM','rate units':'mol m-3','reactions':reactions}],
                'initial state':{'concentrations':{'A':1e-6,'B':2e-6,'C':1e-6},
                                 'parameters':{'source':1e-10,'loss':.001,'custom':.002,'jtest':.001}}}
    destination.mkdir(parents=True,exist_ok=True)
    (destination/'all_rates.json').write_text(json.dumps(document,indent=2)+'\n')
    (destination/'empty.json').write_text(json.dumps({'camp-data':species})+'\n')


if __name__ == '__main__':
    fixtures(Path(sys.argv[1]))
