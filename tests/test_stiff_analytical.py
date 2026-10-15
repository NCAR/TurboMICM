"""Surface, Robertson, HIRES, E5, and Oregonator reference quantities."""
import json
import math
from pathlib import Path
from chemffi import *

ROOT=Path(__file__).parent
REFERENCES=json.loads((ROOT/'data/analytical_reference.json').read_text())


def user_system(names,reactions,values,atol,rtol):
    rx=[]
    for i,(rr,pp) in enumerate(reactions):
        r=rate(9);r.parameter=i;rx.append((rr,pp,r))
    m=Mechanism(names,rx,parameters=[f'r{i+1}' for i in range(len(rx))],atols=atol)
    lib.chem_relative_tolerance(m.state,rtol)
    for i,value in enumerate(values):m.custom[i]=value
    assert lib.chem_update_rates(m.state)==0
    return m


def solve(m,method,dt,parameters):
    if method:return m.solve(dt,parameters=parameters)
    stats=Stats()
    status=lib.chem_solve_backward_euler(m.state,C.byref(parameters),dt,C.byref(stats))
    return status,stats


def advance(m,method,dt,parameters,substeps=1,strict=False):
    elapsed=0.;limit=dt/substeps
    while elapsed<dt:
        status,stats=solve(m,method,min(limit,dt-elapsed),parameters)
        assert status in ([0,3] if strict else [0,3,7]),(status,elapsed,dt)
        advanced=elapsed+stats.final_time
        assert advanced>elapsed,(status,elapsed,dt,stats.final_time)
        elapsed=advanced


def surface(method,tolerance):
    gmd,gsd=1e-6,.1
    radius=gmd/2*math.exp(5*math.log(gsd)*math.log(gsd)/2)
    number=6/(math.pi*gmd**3*math.exp(9/2*math.log(gsd)*math.log(gsd)))*(.002/1000+.003/1000)
    t=272.5
    speed=math.sqrt(8*(1.380649e-23*6.02214076e23)/(math.pi*.04607)*t)
    k=4*number*math.pi*radius*radius/(radius/.95e-5+4/(speed*.02))
    r=rate(10,[.95e-5,.04607,.02]);r.parameter=0
    with Mechanism(['foo','bar','baz'],[([0],[(1,1.),(2,.4)],r)],parameters=['radius','number']) as m:
        m.set_y([[1.,0.,0.]])
        lib.chem_temperatures(m.state)[0]=t;lib.chem_pressures(m.state)[0]=101253.3
        m.custom[0],m.custom[1]=radius,number
        parameters=lib.chem_rosenbrock_parameters(method) if method else lib.chem_backward_euler_parameters()
        dt=.1/k
        for step in range(1,11):
            assert lib.chem_update_rates(m.state)==0
            status,stats=solve(m,method,dt,parameters)
            assert status==0
            assert abs(lib.chem_rate_constants(m.state)[0]-k)<=1e-8
            foo=math.exp(-k*(step*dt))
            for actual,reference in zip(m.rows()[0],[foo,1-foo,.4*(1-foo)]):
                assert abs(actual-reference)<=tolerance*abs(reference),(method,step,actual,reference,tolerance)


def robertson(method,algorithm,tolerance):
    rx=[([0],[(1,1.)]),([1,1],[(1,1.),(2,1.)]),([1,2],[(0,1.),(2,1.)])]
    with user_system(['A','B','C'],rx,[.04,3e7,1e4],[1e-12]*3,1e-10) as m:
        lib.chem_set_linear_solver(m.state,getattr(lib,'chem_linear_'+algorithm)())
        m.set_y([[1.,0.,0.]])
        lib.chem_temperatures(m.state)[0]=272.5;lib.chem_pressures(m.state)[0]=101253.3;lib.chem_air_densities(m.state)[0]=1e6
        params=lib.chem_rosenbrock_parameters(method)
        current,target=0.,1.
        for expected in REFERENCES['Robertson'][1:]:
            advance(m,method,target-current,params)
            for j,(a,b) in enumerate(zip(m.rows()[0],expected)):
                assert abs(a-b)<.3e-3 or abs(a-b)<tolerance*abs(b),(method,algorithm,target,j,a,b)
            current,target=target,target*10


def hires(method,algorithm,tolerance):
    rx=[([0],[(1,1.)]),([1],[(0,.43/8.75),(3,8.32/8.75)]),([2],[(0,8.32/10.03),(3,1.71/10.03)]),
        ([],[(0,1.)]),([3],[(2,.43/1.12),(5,.69/1.12)]),([4],[(2,.035/1.745),(5,1.71/1.745)]),
        ([5],[(4,1.)]),([6],[(4,.43/1.81),(5,.69/1.81),(7,1.)]),([5,7],[(6,1.)])]
    with user_system([f'Y{i}' for i in range(8)],rx,[1.71,8.75,10.03,.0007,1.12,1.745,.43,1.81,280.],[1e-8]*8,1e-6) as m:
        lib.chem_set_linear_solver(m.state,getattr(lib,'chem_linear_'+algorithm)())
        m.set_y([REFERENCES['Hires'][0]])
        params=lib.chem_rosenbrock_parameters(method) if method else lib.chem_backward_euler_parameters()
        dt=321.8122
        for expected in REFERENCES['Hires'][1:]:
            advance(m,method,dt,params)
            for j,(a,b) in enumerate(zip(m.rows()[0],expected)):
                assert abs(a-b)<=tolerance,(method,algorithm,j,a,b,tolerance)
            dt+=100


def e5(method,algorithm,tolerance):
    rx=[([0],[(1,1.),(2,1.)]),([1,2],[(4,1.)]),([0,2],[(3,1.)]),([3],[(2,1.),(5,1.)])]
    atol=[1e-7,1e-17,1e-17,1e-17,1e-7,1e-7]
    with user_system([f'A{i+1}' for i in range(6)],rx,[7.89e-10,1.13e9,1.1e7,1.13e3],atol,1e-13) as m:
        lib.chem_set_linear_solver(m.state,getattr(lib,'chem_linear_'+algorithm)())
        m.set_y([REFERENCES['E5'][0]])
        params=lib.chem_rosenbrock_parameters(method)
        current,target=0.,10.
        for expected in REFERENCES['E5'][1:]:
            advance(m,method,target-current,params)
            for j,(a,b) in enumerate(zip(m.rows()[0][:4],expected[:4])):
                assert abs(a-b)<(1e-6 if j==0 else 1e-13) or abs(a-b)<tolerance*abs(b),(method,algorithm,target,j,a,b)
            current,target=target,target*100


def oregonator(method,algorithm,tolerance,substeps):
    s,q,w=77.27,8.375e-6,.161
    k2,k5=1.6e9,1.
    alpha=k5/(s*w*k2);eta=s*k5/(w*k2);rho=k5*k5/(w*w*k2);tau=w/k5
    rx=[([1],[(0,1.)]),([0,1],[(3,1.)]),([0],[(2,1.),(0,2.)]),([0,0],[(4,1.)]),([2],[(1,1.)])]
    values=[k5/(s*w),k2,s*k5/w,s*s*q*k2/2,k5]
    atol=[alpha*1e-11,eta*1e-11,rho*1e-11,eta,eta]
    with user_system(['X','Y','Z','P','Q'],rx,values,atol,1e-9) as m:
        lib.chem_set_linear_solver(m.state,getattr(lib,'chem_linear_'+algorithm)())
        m.set_y([[alpha,2*eta,3*rho,0.,0.]])
        params=lib.chem_rosenbrock_parameters(method) if method else lib.chem_backward_euler_parameters()
        for step,expected in enumerate(REFERENCES['Oregonator'][1:],1):
            advance(m,method,30*tau,params,substeps=substeps,strict=True)
            for j,(a,b) in enumerate(zip(m.rows()[0][:3],[expected[0]*alpha,expected[1]*eta,expected[2]*rho])):
                assert abs(a-b)<tolerance*abs(b),(method,algorithm,step,j,a,b,tolerance)


def check(name,variants):
    settings={}
    for variant in variants:
        # These source helpers create one logical cell. Runtime groups therefore
        # have identical storage; preserve distinct methods/algorithms/bounds.
        key=(variant['method'],variant['algorithm'],variant.get('substeps',1))
        settings[key]=min(settings.get(key,math.inf),variant['tolerance'])
    for (method,algorithm,substeps),tolerance in settings.items():
        if name=='SurfaceRxn':surface(method,tolerance)
        elif name=='Robertson':robertson(method,algorithm,tolerance)
        elif name=='HIRES':hires(method,algorithm,tolerance)
        elif name=='E5':e5(method,algorithm,tolerance)
        elif name=='Oregonator':oregonator(method,algorithm,tolerance,substeps)
        else:raise AssertionError(name)
