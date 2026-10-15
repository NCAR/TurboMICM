"""Closed-form chain and fast-equilibrium integration case ports.

Every retained source invocation keeps its method, LU algorithm, group length,
cell conditions, 999 solve intervals, and absolute assertion bound. C++ sparse
policy variants map to the common C sparse plan; CSR/CSC storage utilities have
separate comprehensive case ports.
"""
from functools import partial
import json
import math
from pathlib import Path
from chemffi import *

TEMPERATURES = [272.5,254.7,312.6]
PRESSURES = [101253.3,100672.5,101319.8]
GAS = 1.380649e-23 * 6.02214076e23


def rates(kind,stiff):
    if kind=='Troe':return rate(1,[4e-11,0,0,1,0,0,.6,1]),rate(1,[1.2e-3,1.6,3,136,5,24,.9,.8])
    if kind=='TernaryChemicalActivation':return rate(2,[4e-5,0,0,1,0,0,.6,1]),rate(2,[1.2e-3,1.6,3,136,5,24,.9,.8])
    if kind=='Photolysis':return rate(9),rate(9)
    if kind=='Tunneling':return rate(4,[4e-3,0,0]),rate(4,[1.2e-4,167,1e8])
    if kind=='Arrhenius':return rate(0,[4e-3,0,50,300,0]),rate(0,[1.2e-4,1.6 if stiff else 7,75,50,.5])
    if kind=='Branched':return rate(5,[1e-4,204.3,1e-3,2]),rate(6,[1e-4,204.3,1e-3,2])
    raise AssertionError(kind)


def expected_rates(kind,stiff,t,p,air):
    if kind in ['Troe','TernaryChemicalActivation']:
        k0 = 4e-11 if kind=='Troe' else 4e-5
        def falloff(low,high,fc,n):
            ratio=low*air/high
            return low*(air if kind=='Troe' else 1)/(1+ratio)*fc**(n/(n+math.log10(ratio)**2))
        return falloff(k0,1,.6,1),falloff(1.2e-3*math.exp(3/t)*(t/300)**1.6,136*math.exp(24/t)*(t/300)**5,.9,.8)
    if kind=='Photolysis':return .002,.003
    if kind=='Tunneling':return .004,1.2e-4*math.exp(-167/t+1e8/t**3)
    if kind=='Arrhenius':return .004*math.exp(50/t),1.2e-4*math.exp(75/t)*(t/50)**(1.6 if stiff else 7)*(1+.5*p)
    if kind=='Branched':
        a=2e-22*math.exp(2)*2.45e19;b=.43*(293/298)**-8
        z=a/(1+a/b)*.41**(1/(1+math.log10(a/b)**2))*(1-1e-3)/1e-3
        a=2e-22*math.exp(2)*(air*6.02214076e23*1e-6);b=.43*(t/298)**-8
        av=a/(1+a/b)*.41**(1/(1+math.log10(a/b)**2))
        base=1e-4*math.exp(-204.3/t)
        return base*z/(z+av),base*av/(z+av)
    raise AssertionError(kind)


def execute(kind,stiff,method,group,algorithm,tolerance):
    for start in range(0,3,group or 3):
        count=min(group or 3,3-start)
        first,second=rates(kind,stiff)
        names=['A1','A2','B','C'] if stiff else ['A','B','C']
        photo=kind=='Photolysis'
        labels=['photoA1B','photoA2B','photoB'] if stiff and photo else (['photoA','photoB'] if photo else [])
        if stiff:
            first2=Rate.from_buffer_copy(first)
            if photo:first.parameter,first2.parameter,second.parameter=0,1,2
            reactions=[([0],[(2,1.)],first),([1],[(2,1.)],first2),([2],[(3,1.)],second),
                       ([0],[(1,1.)],rate(0,[4e10,0,0,300,0])),([1],[(0,1.)],rate(0,[.9*4e10,0,0,300,0]))]
        else:
            if photo:first.parameter,second.parameter=0,1
            reactions=[([0],[(1,1.)],first),([1],[(2,1.)],second)]
        with Mechanism(names,reactions,parameters=labels,cells=count) as m:
            assert lib.chem_set_linear_solver(m.state,getattr(lib,'chem_linear_'+algorithm)())==0
            initial=[1.-c/3. for c in range(start,start+count)]
            m.set_y([[.5*a,.5*a,0.,0.] if stiff else [a,0.,0.] for a in initial])
            kk=[]
            for c in range(count):
                t,p=TEMPERATURES[start+c],PRESSURES[start+c]
                lib.chem_temperatures(m.state)[c]=t;lib.chem_pressures(m.state)[c]=p
                lib.chem_air_densities(m.state)[c]=p/(GAS*t)
                kk.append(expected_rates(kind,stiff,t,p,p/(GAS*t)))
                if photo:
                    for j,v in enumerate([.002,.002,.003] if stiff else [.002,.003]):m.custom[j*count+c]=v
            parameters=lib.chem_rosenbrock_parameters(method) if method else lib.chem_backward_euler_parameters()
            for step in range(1,1000):
                assert lib.chem_update_rates(m.state)==0
                if not stiff:
                    actual=lib.chem_rate_constants(m.state)
                    for c,(k1,k2) in enumerate(kk):
                        assert abs(actual[c]-k1)<=1e-9*abs(k1)
                        assert abs(actual[count+c]-k2)<=1e-9*abs(k2)
                if method:status,stats=m.solve(1.,parameters=parameters)
                else:
                    stats=Stats();status=lib.chem_solve_backward_euler(m.state,C.byref(parameters),1.,C.byref(stats))
                assert status==0,(kind,stiff,method,group,algorithm,step,status,stats.final_time)
                for c,row in enumerate(m.rows()):
                    k1,k2=kk[c]
                    a=initial[c]
                    expected=[a*math.exp(-k1*step),a*(k1/(k2-k1))*(math.exp(-k1*step)-math.exp(-k2*step)),
                              a*(1+(k1*math.exp(-k2*step)-k2*math.exp(-k1*step))/(k2-k1))]
                    actual=[row[0]+row[1],row[2],row[3]] if stiff else row
                    for j,(value,reference) in enumerate(zip(actual,expected)):
                        assert abs(value-reference)<=tolerance,(kind,stiff,method,group,algorithm,step,start+c,j,value,reference,tolerance)


def check(name,variants):
    if name in ['SurfaceRxn','Robertson','E5','HIRES','Oregonator']:
        from test_stiff_analytical import check as stiff_check
        return stiff_check(name,variants)
    stiff=name.endswith('SuperStiffButAnalytical')
    kind=name.removesuffix('SuperStiffButAnalytical')
    # Several source template policies have the same C runtime representation.
    # Execute every distinct numerical setting with its tightest source bound.
    settings={}
    for variant in variants:
        key=(variant['method'],variant['group'],variant['algorithm'])
        settings[key]=min(settings.get(key,math.inf),variant['tolerance'])
    for (method,group,algorithm),tolerance in settings.items():execute(kind,stiff,method,group,algorithm,tolerance)


CASES={key:partial(check,key.split('.')[-1],variants) for key,variants in
       json.loads((Path(__file__).parent/'data/analytical_variants.json').read_text()).items()}
