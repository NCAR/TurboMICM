"""Closed-form reaction chains, with the reference cells, intervals, and bounds."""
import math
import json
from api import *
from harness import case, equal, near

TABLES = json.loads((ROOT/'tests/analytical-tables.json').read_text())

TEMP = [272.5,254.7,312.6]
PRESSURE = [101253.3,100672.5,101319.8]
KINDS = ('Troe','Photolysis','TernaryChemicalActivation','Tunneling','Arrhenius','Branched')


def chain_rates(kind,stiff):
    if kind == 'Troe': return [rate(TROE,k0_A=4e-11),rate(TROE,k0_A=1.2e-3,k0_B=1.6,k0_C=3,kinf_A=136,kinf_B=5,kinf_C=24,Fc=.9,N=.8)]
    if kind == 'Photolysis': return [rate(USER,label='photoA'),rate(USER,label='photoB')]
    if kind == 'TernaryChemicalActivation': return [rate(TERNARY,k0_A=4e-5),rate(TERNARY,k0_A=1.2e-3,k0_B=1.6,k0_C=3,kinf_A=136,kinf_B=5,kinf_C=24,Fc=.9,N=.8)]
    if kind == 'Tunneling': return [rate(TUNNELING,A=4e-3),rate(TUNNELING,A=1.2e-4,B=167,C=1e8)]
    if kind == 'Arrhenius': return [rate(A=4e-3,C=50),rate(A=1.2e-4,B=1.6 if stiff else 7,C=75,D=50,E=.5)]
    return [rate(BRANCHED,branch=branch,X=1e-4,Y=204.3,a0=1e-3,n=2) for branch in (0,1)]


def independent_rate(r,t,p,d):
    if r.type in (TROE,TERNARY):
        k0 = r.k0_A*math.exp(r.k0_C/t)*(t/300)**r.k0_B
        ki = r.kinf_A*math.exp(r.kinf_C/t)*(t/300)**r.kinf_B
        ratio = k0*d/ki
        return k0*(d if r.type == TROE else 1)/(1+ratio)*r.Fc**(r.N/(r.N+math.log10(ratio)**2))
    if r.type == USER: return 2e-3 if r.label == b'photoA' else 3e-3
    if r.type == TUNNELING: return r.A*math.exp(-r.B/t+r.C/t**3)
    if r.type == ARRHENIUS: return r.A*math.exp(r.C/t)*(t/r.D)**r.B*(1+r.E*p)
    a,b = 2e-22*math.exp(2)*2.45e19,.43*(293/298)**-8
    z = a/(1+a/b)*.41**(1/(1+math.log10(a/b)**2))*(1-1e-3)/1e-3
    a,b = 2e-22*math.exp(2)*(d*6.02214076e23*1e-6),.43*(t/298)**-8
    capital_a = a/(1+a/b)*.41**(1/(1+math.log10(a/b)**2))
    return 1e-4*math.exp(-204.3/t)*((z if not r.branch else capital_a)/(z+capital_a))


def chain_tolerance(kind,stiff,method):
    if method == 0:
        return {'Troe':1e-5 if stiff else 1e-6,'Photolysis':1e-3,'TernaryChemicalActivation':1e-2 if stiff else 1e-5,
                'Tunneling':1e-3,'Arrhenius':1e-3,'Branched':1e-2 if stiff else 1e-5}[kind]
    if stiff:
        if kind == 'Arrhenius': return {2:1e-4,3:2e-5,4:2e-5,5:2e-5,6:1e-5}[method]
        return {'Troe':1e-5,'Photolysis':2e-5,'TernaryChemicalActivation':2e-3,'Tunneling':1e-4,'Branched':2e-3}[kind]
    return {'Troe':1e-10,'Photolysis':2e-6,'TernaryChemicalActivation':1e-8,'Tunneling':2e-5 if method == 2 else 1e-8,
            'Arrhenius':4e-6 if method == 2 else 1e-9,'Branched':1e-10 if method == 2 else 1e-13}[kind]


@case(*[family+':AnalyticalExamples.'+kind+suffix for family in ('analytical_rosenbrock_integration','analytical_backward_euler')
        for kind in KINDS for suffix in ('','SuperStiffButAnalytical')])
def chains(identifier):
    family,name = identifier.split(':'); name = name.split('.')[-1]
    stiff = name.endswith('SuperStiffButAnalytical'); kind = name.removesuffix('SuperStiffButAnalytical')
    rates = chain_rates(kind,stiff)
    for method in ((0,) if family == 'analytical_backward_euler' else (2,3,4,5,6)):
        for doolittle in (False,True):
            m = Model(['A1','A2','B','C'] if stiff else ['A','B','C'])
            if stiff:
                m.reaction(['A1'],[('B',1)],rates[0]); m.reaction(['A2'],[('B',1)],rates[0]); m.reaction(['B'],[('C',1)],rates[1])
                m.reaction(['A1'],[('A2',1)],rate(A=4e10)); m.reaction(['A2'],[('A1',1)],rate(A=.9*4e10))
            else:
                m.reaction(['A'],[('B',1)],rates[0]); m.reaction(['B'],[('C',1)],rates[1])
            m.finalize(); s = State(m,3)
            try:
                for i in range(3):
                    initial = 1-i/3
                    if stiff:
                        s.data.variables[m.var('A1')*3+i] = s.data.variables[m.var('A2')*3+i] = .5*initial
                    else: s.data.variables[m.var('A')*3+i] = initial
                    s.data.conditions[i] = Conditions(TEMP[i],PRESSURE[i],PRESSURE[i]/((1.380649e-23*6.02214076e23)*TEMP[i]))
                if kind == 'Photolysis': s.set('photoA',.002,custom=True); s.set('photoB',.003,custom=True)
                s.rates()
                p = lib.chem_euler_parameters() if method == 0 else lib.chem_solver(method)
                if doolittle: p.factor = C.cast(lib.chem_lu_doolittle,P)
                k = [[independent_rate(r,TEMP[i],PRESSURE[i],s.data.conditions[i].air_density) for r in rates] for i in range(3)]
                tolerance = chain_tolerance(kind,stiff,method)
                if not stiff:
                    for i in range(3):
                        for j in range(2): near(s.data.rates[j*3+i],k[i][j],abs(s.data.rates[j*3+i])*1e-9)
                for time in range(1,1000):
                    s.advance(1,parameters=p)
                    for i,(k1,k2) in enumerate(k):
                        initial = 1-i/3; e1,e2 = math.exp(-k1*time),math.exp(-k2*time)
                        expected = [initial*e1,initial*(k1/(k2-k1))*(e1-e2),initial*(1+(k1*e2-k2*e1)/(k2-k1))]
                        values = [s.get('A1',i)+s.get('A2',i) if stiff else s.get('A',i),s.get('B',i),s.get('C',i)]
                        for actual,reference in zip(values,expected): near(actual,reference,tolerance)
            finally: s.close(); m.close()


@case('analytical_rosenbrock_integration:AnalyticalExamples.SurfaceRxn','analytical_backward_euler:AnalyticalExamples.SurfaceRxn')
def surface(identifier):
    radius = 1e-6/2*math.exp(5*math.log(.1)*math.log(.1)/2)
    number = 6/(math.pi*(1e-6)**3*math.exp(9/2*math.log(.1)*math.log(.1)))*(.002/1000+.003/1000)
    speed = math.sqrt(8*(1.380649e-23*6.02214076e23)/(math.pi*.04607)*272.5)
    k = 4*number*math.pi*radius*radius/(radius/.95e-5+4/(speed*.02)); dt = .1/k
    for method,tolerance in ([(0,.05)] if identifier.startswith('analytical_backward') else [(2,1e-2),(3,1e-5),(4,1e-6),(5,1e-5),(6,1e-7)]):
        m = Model(['foo','bar','baz']); m.reaction(['foo'],[('bar',1),('baz',.4)],rate(SURFACE,label='foo',diffusion=.95e-5,molecular_weight=.04607,probability=.02)); m.finalize(); s = State(m)
        try:
            s.set('foo',1); s.set('foo.effective radius [m]',radius,custom=True); s.set('foo.particle number concentration [# m-3]',number,custom=True)
            s.data.conditions[0] = Conditions(272.5,101253.3,0); s.rates(); near(s.data.rates[0],k,1e-8)
            for step in range(1,11):
                s.advance(dt,method=method); remaining = math.exp(-k*(step*dt))
                for name,expected in [('foo',remaining),('bar',1-remaining),('baz',.4*(1-remaining))]: near(s.get(name),expected,tolerance*abs(expected))
        finally: s.close(); m.close()


def stiff_model(kind):
    if kind == 'Robertson':
        names = ['A','B','C']
        reactions = [(['A'],[('B',1)],.04),(['B','B'],[('B',1),('C',1)],3e7),(['B','C'],[('A',1),('C',1)],1e4)]
        initial,rtol,atol = [1,0,0],1e-10,[1e-12]*3
    elif kind == 'E5':
        names = ['A'+str(i) for i in range(1,7)]
        reactions = [(['A1'],[('A2',1),('A3',1)],7.89e-10),(['A2','A3'],[('A5',1)],1.13e9),
                     (['A1','A3'],[('A4',1)],1.1e7),(['A4'],[('A3',1),('A6',1)],1.13e3)]
        initial,rtol,atol = [1.76e-3,0,0,0,0,0],1e-13,[1e-7,1e-17,1e-17,1e-17,1e-7,1e-7]
    elif kind == 'Hires':
        names = ['Y'+str(i) for i in range(8)]
        reactions = [(['Y0'],[('Y1',1)],1.71),(['Y1'],[('Y0',.43/8.75),('Y3',8.32/8.75)],8.75),
          (['Y2'],[('Y0',8.32/10.03),('Y3',1.71/10.03)],10.03),([],[('Y0',1)],.0007),
          (['Y3'],[('Y2',.43/1.12),('Y5',.69/1.12)],1.12),(['Y4'],[('Y2',.035/1.745),('Y5',1.71/1.745)],1.745),
          (['Y5'],[('Y4',1)],.43),(['Y6'],[('Y4',.43/1.81),('Y5',.69/1.81),('Y7',1)],1.81),(['Y5','Y7'],[('Y6',1)],280)]
        initial,rtol,atol = [1,0,0,0,0,0,0,.0057],1e-6,[1e-8]*8
    else:
        names = ['X','Y','Z','P','Q']; ss,q,w,k2 = 77.27,8.375e-6,.161,1.6e9
        alpha,eta,rho = 1/(ss*w*k2),ss/(w*k2),1/(w*w*k2)
        reactions = [(['Y'],[('X',1)],1/(ss*w)),(['X','Y'],[('P',1)],k2),(['X'],[('Z',1),('X',2)],ss/w),
                     (['X','X'],[('Q',1)],ss*ss*q*k2/2),(['Z'],[('Y',1)],1)]
        initial,rtol,atol = [alpha,2*eta,3*rho,0,0],1e-9,[alpha*1e-11,eta*1e-11,rho*1e-11,eta,eta]
    m = Model(names)
    for i,(reactants,products,k) in enumerate(reactions): m.reaction(reactants,products,rate(USER,label='r'+str(i+1)))
    m.finalize(); s = State(m); s.tolerances(rtol,atol)
    for name,value in zip(names,initial): s.set(name,value)
    s.data.conditions[0] = Conditions(272.5,101253.3,1e6)
    for i,(_,_,k) in enumerate(reactions): s.set('r'+str(i+1),k,custom=True)
    s.rates()
    return m,s,names


@case(*['analytical_rosenbrock_integration:AnalyticalExamples.'+name for name in ('Robertson','E5','Oregonator','HIRES')],
      'analytical_backward_euler:AnalyticalExamples.HIRES','analytical_backward_euler:AnalyticalExamples.Oregonator')
def stiff_problems(identifier):
    kind = identifier.split('.')[-1].replace('HIRES','Hires')
    backward = identifier.startswith('analytical_backward')
    for method in ((0,) if backward else (2,3,4,5,6)):
        for doolittle in (False,True):
            m,s,names = stiff_model(kind)
            try:
                p = lib.chem_euler_parameters() if backward else lib.chem_solver(method)
                if doolittle: p.factor = C.cast(lib.chem_lu_doolittle,P)
                previous = 0
                table = TABLES[kind]
                for i,row in enumerate(table[1:]):
                    if kind == 'Robertson': target = 10**i; dt = target-previous; previous = target
                    elif kind == 'E5': target = 10*100**i; dt = target-previous; previous = target
                    elif kind == 'Hires': dt = 321.8122+100*i
                    else: dt = 30*.161
                    s.advance(dt,parameters=p,substeps=18000 if backward and kind == 'Oregonator' else 1)
                    for j,expected in enumerate(row[:3 if kind == 'Oregonator' else 4 if kind == 'E5' else len(names)]):
                        if kind == 'Oregonator': expected *= [1/(77.27*.161*1.6e9),77.27/(.161*1.6e9),1/(.161*.161*1.6e9)][j]
                        actual = s.get(names[j]); error = abs(actual-expected)
                        if kind == 'Robertson': bound = max(.0003,1e-6*abs(expected))
                        elif kind == 'E5': bound = max(1e-6 if j == 0 else 1e-13,1e-10*abs(expected))
                        elif kind == 'Hires': bound = .1 if backward else {2:1e-6,3:1e-7,4:1e-7,5:1e-6,6:1e-6}[method]
                        else: bound = (5e-3 if backward else 1e-6)*abs(expected)
                        assert error < bound,(kind,method,doolittle,i,j,actual,expected,bound)
            finally: s.close(); m.close()
