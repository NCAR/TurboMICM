"""Compare generated mechanisms with independently assembled C API mechanisms."""
import ctypes as C
from pathlib import Path
import subprocess
import sys

sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'tools'))
from configure import normalize, read_documents, emit_c, ConfigurationError
from translate_ts1 import translate
from chemffi import *

ROOT=Path(__file__).resolve().parent.parent
CASES={}


def make_from_data(data,owner):
    m=object.__new__(Mechanism)
    m.data=data;m.n=data.species_count;m.nc=1;m.keep=[owner]
    m.handle,m.state=P(),P()
    assert lib.chem_mechanism_create(C.byref(data),C.byref(m.handle))==0
    assert lib.chem_state_create(m.handle,1,C.byref(m.state))==0
    m.y=lib.chem_concentrations(m.state);m.custom=lib.chem_parameters(m.state)
    return m


def by_hand(user=False):
    conversion=1e-6*6.02214076e23
    surface=rate(10,[2.3e2,.025,.9]);surface.parameter=0
    reactions=[([0],[(1,1)],rate(0,[2.15e-4,0,110,300,0])),
        ([1],[(2,1)],rate(5,[1.2,204.3,1e-3,2])),
        ([1],[(3,1)],rate(6,[1.2,204.3,1e-3,2])),
        ([2],[(4,1)],surface),
        ([3],[(5,2)],rate(2,[1.2,2.3,302.3,2.6/conversion,-3.1,402.1,.9,1.2])),
        ([4,4],[(6,1)],rate(1,[1.8e-30*conversion*conversion,-3,0,2.8e-11*conversion,0,0,.6,1])),
        ([5],[(6,1)],rate(4,[1.2,2.3,302.3]))]
    labels=['C.effective radius [m]','C.particle number concentration [# m-3]']
    if user:
        for rr,pp,label in [([2],[(6,1)],'my photolysis rate'),([],[(0,1)],'my emission rate'),([1],[],'my loss rate')]:
            r=rate(9);r.parameter=len(labels);labels.append(label);reactions.append((rr,pp,r))
    return Mechanism(list('ABCDEFG'),reactions=reactions,parameters=labels)


def names(m):return [m.data.species_names[i].decode() for i in range(m.n)]

def concentrations(m):return dict(zip(names(m),m.rows()[0]))

def setup(m,user):
    lib.chem_temperatures(m.state)[0]=287.45;lib.chem_pressures(m.state)[0]=101319.9
    for i,name in enumerate(names(m)):m.y[i]=1 if name=='A' else 0
    for j in range(m.data.parameter_count):
        label=m.data.parameter_names[j].decode()
        if label.endswith('.effective radius [m]'):m.custom[j]=1e-7
        elif label.endswith('.particle number concentration [# m-3]'):m.custom[j]=2.5e6
        elif label=='my emission rate':m.custom[j]=1e-20
        elif label=='my loss rate':m.custom[j]=1e-20*1e-3
    assert lib.chem_calculate_air_density(m.state)==0


def advance(m,photo,user):
    if user:
        for j in range(m.data.parameter_count):
            if m.data.parameter_names[j]==b'my photolysis rate':m.custom[j]=photo
    assert lib.chem_update_rates(m.state)==0
    elapsed=0
    while elapsed<500:
        status,stats=m.solve(500-elapsed)
        assert status in [0,3] and stats.final_time>0,(status,stats.final_time)
        elapsed+=stats.final_time


def configuration_equivalence(user):
    owner=C.CDLL(str(ROOT/'build/configs.so'))
    data=Data.in_dll(owner,'user_data' if user else 'no_user_data')
    with make_from_data(data,owner) as configured, by_hand(user) as hand:
        setup(configured,user);setup(hand,user)
        for step in range(10):
            advance(configured,1e-10*1.5**step,user);advance(hand,1e-10*1.5**step,user)
            a,b=concentrations(configured),concentrations(hand)
            for name in a:
                assert abs(a[name]-b[name])<=1e-30+1e-12*abs(b[name]),(step,name,a[name],b[name])


def program(user):
    result=subprocess.run([ROOT/'bin/rate_tutorial']+(['--user'] if user else []),capture_output=True,text=True,check=True)
    expected=(ROOT/'tests/data'/('user_output.txt' if user else 'no_user_output.txt')).read_text()
    assert result.stdout==expected,'tutorial output differs from reference'
    # The build-time reader is independently checked alongside the program case.
    configuration_equivalence(user)


CASES['rate_constants_no_user_defined_example_by_hand']=lambda:program(False)
CASES['rate_constants_user_defined_example_by_hand']=lambda:program(True)


def verify_reader():
    import json
    assert translate(ROOT.parent/'micm/test/benchmark/ts1_mechanism.hpp')==json.loads((ROOT/'configs/ts1.json').read_text())
    assert len(normalize(ROOT/'configs/ts1.json')['reactions'])==547
    m=normalize(ROOT/'configs/robertson/config.json')
    assert len(m['species'])==3 and len(m['reactions'])==3
    # JSON/YAML must have the same data semantics, including reaction ordering.
    directory=ROOT/'build/config-reader-test';directory.mkdir(exist_ok=True)
    source=read_documents(ROOT/'configs/rate_constants_no_user_defined/config.yaml')
    path=directory/'same.json';path.write_text(json.dumps({'camp-data':source}))
    assert normalize(path)==normalize(ROOT/'configs/rate_constants_no_user_defined/config.yaml')
    (directory/'cycle.json').write_text('{"camp-files":["cycle.json"]}')
    try:read_documents(directory/'cycle.json')
    except ConfigurationError:pass
    else:raise AssertionError('include cycle was accepted')


def robertson_equivalence():
    owner=C.CDLL(str(ROOT/'build/configs.so'))
    data=Data.in_dll(owner,'robertson_data')
    reactions=[]
    for index,(rr,pp) in enumerate([([0],[(1,1)]),([1,1],[(1,1),(2,1)]),([1,2],[(0,1),(2,1)])]):
        r=rate(9);r.parameter=index;reactions.append((rr,pp,r))
    with make_from_data(data,owner) as configured, Mechanism(['A','B','C'],reactions=reactions,parameters=['r1','r2','r3']) as hand:
        for m in [configured,hand]:
            for i,name in enumerate(names(m)):m.y[i]=1 if name=='A' else 0
            for i in range(3):m.custom[i]=[.04,3e7,1e4][i]
            assert lib.chem_update_rates(m.state)==0
        for step in range(10):
            for m in [configured,hand]:
                elapsed=0
                while elapsed<200:
                    status,stats=m.solve(200-elapsed)
                    assert status in [0,3] and stats.final_time>0
                    elapsed+=stats.final_time
            a,b=concentrations(configured),concentrations(hand)
            for name in a:assert abs(a[name]-b[name])<=1e-30+1e-12*abs(b[name]),(step,name,a[name],b[name])


if __name__=='__main__':
    verify_reader();program(False);program(True);robertson_equivalence()
    print('configuration checks PASS')
