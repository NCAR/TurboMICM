"""Equilibrium evolution, conservation, and algebraic error regressions."""
import math
import sys
from api import *
from harness import case, equal, near


def keq(k,dh,t): return k*math.exp(dh/(1.380649e-23*6.02214076e23)*(1/t-1/298.15))


def equilibrium(m,reactant,product,k,name,dh=0,power=1):
    m.constraint(product,[(reactant,power),(product,1)],k,equilibrium_reactants=1,name=name,delta_enthalpy=dh)


def parameter(s,name): return s.data.custom[lib.chem_custom_id(s.model.p,name.encode())*s.data.stride]


@case(*['equilibrium_constraint:EquilibriumIntegration.'+name for name in ('SetConstraintsAPIWorks','SetConstraintsAPIMultipleConstraints','DAESolveWithConstraint','DAESolveWithConstraintAndReorderState','DAESolveWithTwoCoupledConstraints','DAESolveWithNonUnitStoichiometry')])
def equilibrium_cases(identifier):
    name = identifier.split('.')[-1]; metadata = name.startswith('SetConstraints'); multiple = name == 'SetConstraintsAPIMultipleConstraints'
    coupled,nonunit = 'TwoCoupled' in name,'NonUnit' in name
    names = list('ABCDEF' if multiple else 'ABCD' if coupled else 'ABC'); m = Model(names)
    if multiple:
        m.reaction(['A'],[('B',1)],rate(A=.5)); m.reaction(['D'],[('E',1)],rate(A=.2))
        constraints = [('B','C',.034,'B_C_eq',-2400),('E','F',.012,'E_F_eq',-2000)]
        initial,t = [1,.1,.034*.1,.5,.05,.012*.05],310
    elif nonunit:
        m.reaction(['C'],[('A',1)],rate(A=.5)); constraints = [('A','B',10,'A2_B_eq',-2400)]
        initial,t = [.1,.1,1],298
    else:
        m.reaction(['A'],[('B',1)],rate(A=.1 if metadata else 1))
        constraints = [('B','C',.034 if metadata else 3 if coupled else 2,'B_C_eq',-2400)]
        if coupled: constraints.append(('B','D',5,'B_D_eq',-2400))
        initial = ([1,.1,.034*.1] if metadata else [1]+[0]*(len(names)-1))
        t = 300 if metadata or coupled else 400 if 'Reorder' in name else 270
    for r,p,k,label,dh in constraints: equilibrium(m,r,p,k,label,dh,2 if nonunit else 1)
    m.finalize(reorder='Reorder' in name); s = State(m)
    try:
        equal(m.data.nvariables,len(names)); equal(m.data.nconstraints,len(constraints))
        for n in names:
            assert m.var(n) != MAX
            equal(m.data.algebraic[m.var(n)],int(n in [q[1] for q in constraints]))
        for n,x in zip(names,initial): s.set(n,x)
        s.data.conditions[0] = Conditions(t,101325,0); s.rates()
        for _,_,k,label,dh in constraints:
            expected = keq(k,dh,t)
            near(parameter(s,label),expected,abs(expected)*100*sys.float_info.epsilon)
        if not metadata:
            time,stop = 0,.05 if nonunit else .1
            while time < stop:
                s.rates(); s.solve(.001)
                for r,p,k,label,dh in constraints: near(parameter(s,label)*s.get(r)**(2 if nonunit else 1)-s.get(p),0,1e-5 if nonunit else 1e-6)
                time += .001
            if not nonunit: near(s.get('A')+s.get('B'),1,.01)
    finally: s.close(); m.close()


@case('linear_constraint:DAESolveWithConstraint.TerminatorAndRobertson')
def terminator_robertson(_):
    m = Model(['Cl2','Cl','A','B','C'])
    for i in range(2): m.data.species[i].absolute_tolerance = 1e-20
    m.reaction(['Cl2'],[('Cl',2)],rate(USER,label='terminator_k1')); m.reaction(['Cl','Cl'],[('Cl2',1)],rate(A=1))
    m.reaction(['A'],[('B',1)],rate(USER,label='robertson_r1')); m.reaction(['B','B'],[('B',1),('C',1)],rate(USER,label='robertson_r2'))
    m.reaction(['B','C'],[('A',1),('C',1)],rate(USER,label='robertson_r3')); m.constraint('C',[(n,1) for n in 'ABC'],1)
    m.finalize(); s = State(m)
    try:
        s.tolerances(1e-8,[1e-20,1e-20,1e-3,1e-3,1e-3]); s.set('Cl',1.2e-6); s.set('Cl2',1.8e-10); s.set('A',1)
        for n,k in [('robertson_r1',.04),('robertson_r2',3e7),('robertson_r3',1e4)]: s.set(n,k,custom=True)
        s.data.conditions[0] = Conditions(298,101300,42); s.rates(); chlorine = 2*1.8e-10+1.2e-6
        for i in range(12):
            s.advance(10**i,method=5)
            near(sum(s.get(n) for n in 'ABC'),1,1e-10)
            near(2*s.get('Cl2')+s.get('Cl'),chlorine,abs(chlorine)*1e-6)
            for n in ['Cl','Cl2','A','B','C']: assert math.isfinite(s.get(n)) and s.get(n) >= 0,(i,n,s.get(n))
            if i == 0: assert s.get('A') < .99 and s.get('B') > 0 and s.get('C') > 1e-3
        assert s.get('C') > .95 and s.get('A') < .05 and s.get('Cl2') > .4*chlorine and s.get('Cl') < .1*chlorine
    finally: s.close(); m.close()


def conservation_model(cascade=False,k=1e3,k1=10,k2=0):
    names = ['A_gas','A_aq']+(['B_aq'] if cascade else [])+['P']
    m = Model(names); m.reaction(['B_aq' if cascade else 'A_aq'],[('P',1)],rate(A=k))
    equilibrium(m,'A_gas','A_aq',k1,'eq1')
    if cascade: equilibrium(m,'A_aq','B_aq',k2,'eq2')
    m.constraint('A_gas',[(n,1) for n in names[1:]+names[:1]],1e-6); m.finalize()
    return m


@case(*['dae_constraint_overshoot:DAEConstraintOvershoot.'+n for n in ('AlgebraicVariableStaysNonNegative','EquilibriumPlusConservation','AllRosenbrockOrdersConstrained')])
def overshoot(identifier):
    name = identifier.split('.')[-1]; eq = name == 'EquilibriumPlusConservation'
    for method in ((3,4,5,6) if name == 'AllRosenbrockOrdersConstrained' else (5,)):
        if eq: m = conservation_model()
        else:
            m = Model(['A','B','C']); m.reaction(['A'],[('B',1)],rate(A=1e4)); m.constraint('C',[(n,1) for n in 'ABC'],1e-6); m.finalize()
        s = State(m)
        try:
            if eq:
                s.set('A_gas',1e-6/11); s.set('A_aq',10e-6/11); s.tolerances(1e-6,[1e-8,1e-8,1e-12])
            else: s.set('A',.9e-6); s.set('C',.1e-6); s.tolerances(1e-6,[1e-12]*3)
            s.data.conditions[0] = Conditions(298,101325,0); s.rates(); s.solve(30,method)
            near(sum(s.data.variables[:3]),1e-6,1e-12)
            if eq:
                assert s.get('A_gas') >= -1e-18 and s.get('A_aq') >= -1e-18
                near(s.get('P'),1e-6,1e-8)
            else:
                assert s.get('C') >= -1e-18 and s.get('A') < 1e-12
                near(s.get('B'),.9e-6,1e-10); near(s.get('C'),.1e-6,1e-10)
        finally: s.close(); m.close()


def cascade(balance_atol,k,k1,k2):
    m = conservation_model(True,k,k1,k2); s = State(m)
    try:
        s.tolerances(1e-6,[balance_atol,1e-8,1e-8,1e-12]); gas = 1e-6/(1+k1+k1*k2)
        s.set('A_gas',gas); s.set('A_aq',k1*gas); s.set('B_aq',k1*k2*gas)
        s.data.conditions[0] = Conditions(298,101325,0); s.rates(); stats = s.solve(30,5)
        return s.get('A_gas'),s.get('P'),stats.accepted
    finally: s.close(); m.close()


@case('dae_algebraic_error_insensitivity:DAEAlgebraicError.ErrorSensitiveToBalanceAtol','dae_algebraic_error_insensitivity:DAEAlgebraicError.AlgebraicVariableDoesNotOvershootDeeply')
def algebraic_error(identifier):
    if identifier.endswith('ErrorSensitiveToBalanceAtol'):
        loose,tight = cascade(1e-3,100,2,2),cascade(1e-8,100,2,2)
        assert tight[2] > loose[2],(tight,loose)
    else:
        gas,product,steps = cascade(1e-8,1e4,100,50)
        assert gas >= -1e-8; near(product,1e-6,1e-8); near(gas,0,1e-10)
