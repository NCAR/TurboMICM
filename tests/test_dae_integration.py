"""Conservation, equilibrium integration, and algebraic-error regressions."""
from functools import partial
import math
from chemffi import *
from test_state import lib

GAS=1.380649e-23*6.02214076e23
EPS=2.220446049250313e-16


def arr(a):return rate(0,[a,0,0,300,0])


def eq(row,rr,pp,k,enthalpy=0,parameter=-1,name=None):
    c,rr,pp=equilibrium(row,rr,pp,k,enthalpy)
    c.parameter=parameter
    c.name=name.encode() if name else None
    return c,rr,pp


def setup(m,rows,temperature=298.,rtol=1e-6):
    m.set_y(rows)
    for c in range(m.nc):lib.chem_temperatures(m.state)[c]=temperature;lib.chem_pressures(m.state)[c]=101325.
    lib.chem_relative_tolerance(m.state,rtol)
    assert lib.chem_update_rates(m.state)==0


def integrate(m,dt,method,check=lambda row:None,partial_ok=False):
    advanced=0.;accepted=rejected=0
    while advanced<dt:
        status,stats=m.solve(dt-advanced,method=method)
        assert status in ([0,3] if partial_ok else [0]),(status,advanced,dt,stats.final_time)
        assert advanced+stats.final_time>advanced
        advanced+=stats.final_time;accepted+=stats.accepted;rejected+=stats.rejected
        for row in m.rows():check(row)
    return accepted,rejected


def cascade(atol,k,k1,k2):
    constraints=[eq(1,[(0,1.)],[(1,1.)],k1),eq(2,[(1,1.)],[(2,1.)],k2),
                 linear(0,[(1,1.),(2,1.),(3,1.),(0,1.)],1e-6)]
    with Mechanism(['A_gas','A_aq','B_aq','P'],[([2],[(3,1.)],arr(k))],constraints,atols=[atol,1e-8,1e-8,1e-12]) as m:
        ag=1e-6/(1+k1+k1*k2)
        setup(m,[[ag,k1*ag,k2*k1*ag,0.]])
        minimum=[math.inf]
        def check(row):minimum[0]=min(minimum[0],row[0])
        accepted,rejected=integrate(m,30.,-4,check)
        row=m.rows()[0]
        return row[0],row[3],minimum[0],accepted,rejected


def sensitivity():
    loose=cascade(1e-3,100.,2.,2.)
    tight=cascade(1e-8,100.,2.,2.)
    assert tight[3]>loose[3],(loose,tight)


def cascade_overshoot():
    ag,p,minimum,accepted,rejected=cascade(1e-8,1e4,100.,50.)
    assert minimum>=-1e-8 and abs(p-1e-6)<=1e-8 and abs(ag)<=1e-10,(ag,p,minimum)


def conservation(method=-4):
    with Mechanism(['A','B','C'],[([0],[(1,1.)],arr(1e4))],[linear(2,[(0,1.),(1,1.),(2,1.)],1e-6)],atols=[1e-12]*3) as m:
        setup(m,[[.9e-6,0.,.1e-6]])
        def check(row):
            assert abs(sum(row)-1e-6)<=1e-12
            assert row[2]>=-1e-18
        integrate(m,30.,method,check)
        a,b,c=m.rows()[0]
        assert a<1e-12 and abs(b-.9e-6)<=1e-10 and abs(c-.1e-6)<=1e-10


def all_orders():
    for method in [3,4,-4,6]:conservation(method)


def equilibrium_conservation():
    constraints=[eq(1,[(0,1.)],[(1,1.)],10.),linear(0,[(1,1.),(2,1.),(0,1.)],1e-6)]
    with Mechanism(['A_gas','A_aq','P'],[([1],[(2,1.)],arr(1e3))],constraints,atols=[1e-8,1e-8,1e-12]) as m:
        ag=1e-6/11
        setup(m,[[ag,10*ag,0.]])
        def check(row):
            assert abs(sum(row)-1e-6)<=1e-12
            assert row[0]>=-1e-18 and row[1]>=-1e-18
        integrate(m,30.,-4,check)
        row=m.rows()[0]
        assert abs(row[2]-1e-6)<=1e-8 and row[0]>=-1e-18


def constant(k,h,t):return k*math.exp((h/GAS)*(1/t-1/298.15))


def metadata(m,algebraic,labels):
    assert m.data.species_count==m.n and m.data.constraint_count==len(algebraic)
    assert [lib.chem_is_algebraic(m.handle,j) for j in range(m.n)]==[int(j in algebraic) for j in range(m.n)]
    for i,name in enumerate(m.data.species_names[:m.n]):assert lib.chem_species_index(m.handle,name)==i
    for label in labels:assert lib.chem_parameter_index(m.handle,label.encode())!=Z(-1).value


def api(multiple=False):
    if multiple:
        constraints=[eq(2,[(1,1.)],[(2,1.)],.034,-2400,0,'B_C_eq'),eq(5,[(4,1.)],[(5,1.)],.012,-2000,1,'E_F_eq')]
        with Mechanism(list('ABCDEF'),[([0],[(1,1.)],arr(.5)),([3],[(4,1.)],arr(.2))],constraints,parameters=['B_C_eq','E_F_eq']) as m:
            metadata(m,[2,5],['B_C_eq','E_F_eq'])
            setup(m,[[1.,.1,.034*.1,.5,.05,.012*.05]],310.)
            expected=[constant(.034,-2400,310.),constant(.012,-2000,310.)]
            assert abs(m.custom[0]-expected[0])<=1e-10
            assert abs(m.custom[1]-expected[1])<=abs(expected[1])*100*EPS
    else:
        constraints=[eq(2,[(1,1.)],[(2,1.)],.034,-2400,0,'B_C_eq')]
        with Mechanism(list('ABC'),[([0],[(1,1.)],arr(.1))],constraints,parameters=['B_C_eq']) as m:
            metadata(m,[2],['B_C_eq'])
            setup(m,[[1.,.1,.034*.1]],300.)
            expected=constant(.034,-2400,300.)
            assert abs(m.custom[0]-expected)<=abs(expected)*100*EPS


def equilibrium_integration(reorder=False,coupled=False,nonunit=False):
    if nonunit:
        names=list('ABC');rx=[([2],[(0,1.)],arr(.5))]
        constraints=[eq(1,[(0,2.)],[(1,1.)],10.,-2400,0,'A2_B_eq')]
        labels=['A2_B_eq'];initial=[.1,10*.1*.1,1.];temp=298.;total=.05
    else:
        names=list('ABCD') if coupled else list('ABC')
        rx=[([0],[(1,1.)],arr(1.))]
        constraints=[eq(2,[(1,1.)],[(2,1.)],3. if coupled else 2.,-2400,0,'B_C_eq')]
        labels=['B_C_eq']
        if coupled:
            constraints.append(eq(3,[(1,1.)],[(3,1.)],5.,-2400,1,'B_D_eq'));labels.append('B_D_eq')
        initial=[1.,0.,0.,0.] if coupled else [1.,0.,0.]
        temp=300. if coupled else (400. if reorder else 270.);total=.1
    if reorder:
        # Exercise the generic runtime reordering API and apply its permutation
        # to names, rates, constraints, and initial values consistently.
        from test_linear import lib as linear_lib
        entries=[(0,0),(1,0),(2,1),(2,2)]
        coords=(MatrixIndex*len(entries))(*[MatrixIndex(*ij) for ij in entries]);perm=(Z*3)()
        assert linear_lib.chem_diagonal_markowitz(3,len(entries),coords,perm)==0
        inv={old:new for new,old in enumerate(perm)}
        names=[names[i] for i in perm];initial=[initial[i] for i in perm]
        rx=[([inv[i] for i in rr],[(inv[i],y) for i,y in pp],r) for rr,pp,r in rx]
        constraints=[(c,[(inv[i],y) for i,y in rr],[(inv[i],y) for i,y in pp]) for c,rr,pp in constraints]
        for c,rr,pp in constraints:c.algebraic_species=inv[c.algebraic_species]
    with Mechanism(names,rx,constraints,parameters=labels) as m:
        setup(m,[initial],temp)
        for c,rr,pp in constraints:
            expected=constant(c.equilibrium_constant,c.enthalpy,temp)
            tol=abs(expected)*100*EPS if not (coupled or reorder or nonunit) else 1e-10
            assert abs(m.custom[c.parameter]-expected)<=tol
        time=0.
        while time<total:
            assert lib.chem_update_rates(m.state)==0
            status,stats=m.solve(.001)
            assert status==0,(reorder,coupled,nonunit,status)
            row=m.rows()[0]
            for c,rr,pp in constraints:
                left=m.custom[c.parameter]
                for i,y in rr:left*=row[i]**y
                right=1.
                for i,y in pp:right*=row[i]**y
                assert abs(left-right)<=(1e-5 if nonunit else 1e-6)
            time+=.001
        row=dict(zip(names,m.rows()[0]))
        if not nonunit:assert abs(row['A']+row['B']-1.)<=.01


def chlorine_robertson():
    photo=rate(9);photo.parameter=0
    rx=[([0],[(1,2.)],photo),([1,1],[(0,1.)],arr(1.))]
    for i,(rr,pp) in enumerate([([2],[(3,1.)]),([3,3],[(3,1.),(4,1.)]),([3,4],[(2,1.),(4,1.)])],1):
        r=rate(9);r.parameter=i;rx.append((rr,pp,r))
    constraints=[linear(4,[(2,1.),(3,1.),(4,1.)],1.)]
    with Mechanism(['Cl2','Cl','A','B','C'],rx,constraints,parameters=['terminator_k1','r1','r2','r3'],atols=[1e-20,1e-20,1e-3,1e-3,1e-3]) as m:
        for i,v in enumerate([0.,.04,3e7,1e4]):m.custom[i]=v
        setup(m,[[1.8e-10,1.2e-6,1.,0.,0.]],298.,1e-8)
        lib.chem_pressures(m.state)[0]=101300.;lib.chem_air_densities(m.state)[0]=42.
        chlorine=2*1.8e-10+1.2e-6
        dt=1.
        for step in range(12):
            integrate(m,dt,-4,partial_ok=True)
            cl2,cl,a,b,c=m.rows()[0]
            assert abs(a+b+c-1.)<=1e-10
            assert abs(2*cl2+cl-chlorine)<=1e-6*abs(chlorine)
            assert all(math.isfinite(v) and v>=0 for v in [cl2,cl,a,b,c])
            if step==0:assert a<.99 and b>0 and c>1e-3
            dt*=10
        assert c>.95 and a<.05 and cl2>.4*chlorine and cl<.1*chlorine


CASES={
    'DAEAlgebraicError.ErrorSensitiveToBalanceAtol':sensitivity,
    'DAEAlgebraicError.AlgebraicVariableDoesNotOvershootDeeply':cascade_overshoot,
    'DAEConstraintOvershoot.AlgebraicVariableStaysNonNegative':conservation,
    'DAEConstraintOvershoot.EquilibriumPlusConservation':equilibrium_conservation,
    'DAEConstraintOvershoot.AllRosenbrockOrdersConstrained':all_orders,
    'EquilibriumIntegration.SetConstraintsAPIWorks':api,
    'EquilibriumIntegration.SetConstraintsAPIMultipleConstraints':partial(api,True),
    'EquilibriumIntegration.DAESolveWithConstraint':equilibrium_integration,
    'EquilibriumIntegration.DAESolveWithConstraintAndReorderState':partial(equilibrium_integration,reorder=True),
    'EquilibriumIntegration.DAESolveWithTwoCoupledConstraints':partial(equilibrium_integration,coupled=True),
    'EquilibriumIntegration.DAESolveWithNonUnitStoichiometry':partial(equilibrium_integration,nonunit=True),
    'DAESolveWithConstraint.TerminatorAndRobertson':chlorine_robertson,
}
