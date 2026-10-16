"""DAE integration: equilibrium, conservation, and algebraic error control."""
import math
from ffi import *
from harness import case,group
from cases_analytical import RGAS,add_user,robertson,advance


def arrhenius(m,src,dst,k):
    rate = rate_default(0); rate.p[0] = k; m.reaction([src],{dst:1},rate)


@group('equilibrium_constraint:EquilibriumIntegration.','''SetConstraintsAPIWorks SetConstraintsAPIMultipleConstraints
DAESolveWithConstraint DAESolveWithConstraintAndReorderState DAESolveWithTwoCoupledConstraints DAESolveWithNonUnitStoichiometry''')
def equilibrium(id_):
    name = id_.split('.')[-1]
    multiple = name == 'SetConstraintsAPIMultipleConstraints'; coupled = name.endswith('CoupledConstraints')
    nonlinear = name.endswith('Stoichiometry'); metadata = name.startswith('SetConstraints')
    names = list('ABCDEF' if multiple else 'ABCD' if coupled else 'ABC')
    temperature = 310 if multiple else 300 if metadata or coupled else 298 if nonlinear else 400 if 'Reorder' in name else 270
    with Model(names) as m:
        arrhenius(m,'C' if nonlinear else 'A','A' if nonlinear else 'B',.5 if nonlinear or multiple else .1 if metadata else 1.)
        if multiple: arrhenius(m,'D','E',.2)
        specifications = [('A2_B_eq','B',{'A':2},10.,-2400)] if nonlinear else [('B_C_eq','C',{'B':1},.034 if metadata else 3. if coupled else 2.,-2400)]
        if multiple: specifications.append(('E_F_eq','F',{'E':1},.012,-2000))
        if coupled: specifications.append(('B_D_eq','D',{'B':1},5.,-2400))
        for label,row,reactants,k,enthalpy in specifications:
            m.equilibrium(row,reactants,{row:1},k,enthalpy,name=label)
        m.finalize(reorder='Reorder' in name)
        assert model_species_count(m.handle) == len(names)
        assert model_constraint_count(m.handle) == len(specifications)
        for n in names:
            assert bool(model_is_algebraic(m.handle,m.indices[n])) == (n in [s[1] for s in specifications])
        with State(m) as s:
            s.conditions(temperature); s.set('A',.1 if nonlinear else 1)
            if nonlinear: s.set('B',.1); s.set('C',1)
            if metadata:
                s.set('B',.1); s.set('C',.0034)
                if multiple: s.set('D',.5); s.set('E',.05); s.set('F',.0006)
            s.update()
            constants = {}
            for label,row,rs,k,h in specifications:
                idx = model_parameter_index(m.handle,label.encode()); assert idx >= 0
                expected = k*math.exp((h/RGAS)*(1/temperature-1/298.15))
                near(s.custom[idx,0],expected,rtol=100*np.finfo(float).eps)
                constants[label] = expected
            if metadata: return
            for _ in range(50 if nonlinear else 100):
                s.update(); s.solve(.001)
                for label,row,rs,k,h in specifications:
                    target = constants[label]*math.prod(s.y[m.indices[n],0]**v for n,v in rs.items())
                    near(s.y[m.indices[row],0],target,rtol=0,atol=1e-5 if nonlinear else 1e-6)
            if not nonlinear: near(s.y[m.indices['A'],0]+s.y[m.indices['B'],0],1.,rtol=0,atol=.01)


@group('dae_constraint_overshoot:DAEConstraintOvershoot.','AlgebraicVariableStaysNonNegative EquilibriumPlusConservation AllRosenbrockOrdersConstrained')
def overshoot(id_):
    eq = id_.endswith('EquilibriumPlusConservation')
    methods = (3,4,43,6) if id_.endswith('AllRosenbrockOrdersConstrained') else (43,)
    for method in methods:
        with Model(['A','B','C']) as m:
            if eq:
                arrhenius(m,'B','C',1e3)
                m.equilibrium('B',{'A':1},{'B':1},10,name='eq')
                m.linear('A',{'B':1,'C':1,'A':1},1e-6)
            else:
                arrhenius(m,'A','B',1e4); m.linear('C',{'A':1,'B':1,'C':1},1e-6)
            m.finalize()
            with State(m) as s:
                s.conditions(298); s.y[:,0] = [1e-6/11,1e-5/11,0] if eq else [.9e-6,0,.1e-6]
                s.tolerance(1e-6,[1e-8,1e-8,1e-12] if eq else 1e-12); s.update(); s.solve(30,method)
                near(s.y.sum(),1e-6,rtol=0,atol=1e-12)
                if eq:
                    assert np.all(s.y[:2] >= -1e-18)
                    near(s.y[2,0],1e-6,rtol=0,atol=1e-8)
                else:
                    assert s.y[2,0] >= -1e-18 and s.y[0,0] < 1e-12
                    near(s.y[1:,0],[.9e-6,.1e-6],rtol=0,atol=1e-10)


def cascade(atol,k,k1,k2):
    with Model(['A_gas','A_aq','B_aq','P']) as m:
        arrhenius(m,'B_aq','P',k)
        m.equilibrium('A_aq',{'A_gas':1},{'A_aq':1},k1,name='eq1')
        m.equilibrium('B_aq',{'A_aq':1},{'B_aq':1},k2,name='eq2')
        m.linear('A_gas',{'A_aq':1,'B_aq':1,'P':1,'A_gas':1},1e-6); m.finalize()
        with State(m) as s:
            s.tolerance(1e-6,[atol,1e-8,1e-8,1e-12]); s.conditions(298)
            v = 1e-6/(1+k1+k1*k2); s.y[:,0] = [v,k1*v,k2*k1*v,0]
            s.update(); _,stats = s.solve(30,43)
            return s.y[:,0].copy(),stats


@group('dae_algebraic_error_insensitivity:DAEAlgebraicError.','ErrorSensitiveToBalanceAtol AlgebraicVariableDoesNotOvershootDeeply')
def sensitivity(id_):
    if id_.endswith('ErrorSensitiveToBalanceAtol'):
        _,loose = cascade(1e-3,100,2,2); _,tight = cascade(1e-8,100,2,2)
        assert tight.accepted > loose.accepted
    else:
        y,_ = cascade(1e-8,1e4,100,50)
        assert y[0] >= -1e-8
        near(y[3],1e-6,rtol=0,atol=1e-8); near(y[0],0.,rtol=0,atol=1e-10)


@case('linear_constraint:DAESolveWithConstraint.TerminatorAndRobertson')
def constrained_robertson(_):
    with Model() as m:
        m.species('Cl2',atol=1e-20); m.species('Cl',atol=1e-20)
        for n in 'ABC': m.species(n)
        r = rate_default(0); m.reaction(['Cl','Cl'],{'Cl2':1},r)
        add_user(m,['Cl2'],{'Cl':2},0.,'terminator_k1'); robertson(m)
        m.linear('C',{'A':1,'B':1,'C':1},1.); m.finalize()
        with State(m) as s:
            s.y[:,0] = [1.8e-10,1.2e-6,1,0,0]; s.tolerance(1e-8,[1e-20,1e-20,.001,.001,.001])
            s.conditions(298,101300,42); s.update(); chlorine = 2*1.8e-10+1.2e-6
            for i in range(12):
                advance(s,10.**i,43); y = s.y[:,0]
                near(sum(y[2:]),1.,rtol=0,atol=1e-10); near(2*y[0]+y[1],chlorine,rtol=1e-6)
                assert np.all(np.isfinite(y)) and np.all(y>=0)
                if i == 0: assert y[2] < .99 and y[3] > 0 and y[4] > 1e-3
            assert y[4] > .95 and y[2] < .05 and y[0] > .4*chlorine and y[1] < .1*chlorine
