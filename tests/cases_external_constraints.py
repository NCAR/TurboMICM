"""External process and constraint composition, dispatch, and Jacobian checks."""
import math
from contextlib import ExitStack
from ffi import *
from harness import case,group
from external import ExternalFixture
from cases_dae import arrhenius

PREFIX = 'external_model_constraints:ExternalModelConstraints.'
GAS,AQ = 'A_GAS','AEROSOL.A_AQ'


def transfer(m,k,total=None,solvent=False):
    return ExternalFixture(m,[([GAS,'S'] if solvent else [GAS],{GAS:-1,AQ:1},k)],
                            [(AQ,{GAS:1,AQ:1},total)] if total is not None else [])


@group(PREFIX,'''AddExternalModelWithConstraints AddExternalModelProcessOnly DAESolveEnforcesConservation
AlgebraicRowHoldsOnlyConstraintTerms CombinedBuiltInAndExternalConstraints AddExternalModelOnlyStandardRosenbrock
AddExternalModelConstraintsOnly MultiGridCell ProcessJacobianElementInAlgebraicRowSurvivesFiltering
ProcessJacobianElementInBuiltInAlgebraicRowSurvivesFiltering''')
def aerosol_constraints(id_):
    name = id_.split('.')[-1]; solvent = name == 'ProcessJacobianElementInAlgebraicRowSurvivesFiltering'
    builtin = 'BuiltInAlgebraic' in name; combined = name == 'CombinedBuiltInAndExternalConstraints'
    process_only = name in ('AddExternalModelProcessOnly','AddExternalModelOnlyStandardRosenbrock')
    k = 10. if name == 'AlgebraicRowHoldsOnlyConstraintTerms' else .1 if solvent or builtin or name in ('MultiGridCell','DAESolveEnforcesConservation') else .01
    with Model([GAS]+(['S'] if solvent else ['B','C'] if combined else ['B'] if builtin else [])+[AQ]) as m:
        external = transfer(m,k,None if process_only or builtin else 1.,solvent)
        if combined:
            arrhenius(m,GAS,'B',.05); m.equilibrium('C',{'B':1},{'C':1},5.)
        if builtin: m.equilibrium(AQ,{'B':1},{AQ:1},5.)
        m.finalize()
        assert model_constraint_count(m.handle) == (0 if process_only else 2 if combined else 1)
        assert not model_is_algebraic(m.handle,m.indices[GAS])
        assert bool(model_is_algebraic(m.handle,m.indices[AQ])) == (not process_only)
        if combined: assert model_is_algebraic(m.handle,m.indices['C'])
        if builtin:
            assert model_jacobian_index(m.handle,m.indices[AQ],m.indices[GAS]) >= 0
            return
        if name.startswith('AddExternalModel') and name != 'AddExternalModelConstraintsOnly':
            assert len(m.names) == 2
            return
        if combined: return
        nc = 3 if name == 'MultiGridCell' else 1
        with State(m,nc) as s:
            s.conditions(298)
            if name == 'AlgebraicRowHoldsOnlyConstraintTerms':
                s.set(GAS,.7); s.set(AQ,.2)
                near(s.forcing(True)[:,0],[-7.,-.1],rtol=0,atol=1e-12)
                near(s.jacobian(True)[0],[[-10.,0.],[1.,1.]],rtol=0,atol=1e-12)
                return
            initial = .8 if name == 'AddExternalModelConstraintsOnly' else .9-.2*np.arange(nc)
            s.set(GAS,initial); s.set(AQ,1-initial)
            if solvent:
                s.set('S',1.)
                assert model_jacobian_index(m.handle,m.indices[AQ],m.indices['S']) >= 0
            steps = 10 if solvent else 20 if name == 'DAESolveEnforcesConservation' else 1
            for _ in range(steps):
                s.solve(10 if steps > 1 else 50,43)
                near(s.y[m.indices[GAS]]+s.y[m.indices[AQ]],np.ones(nc),rtol=0,atol=1e-4)


def equilibrium_model(kind,extra=False,noop=None):
    m = Model(list('ABCD' if extra else 'ABC')); arrhenius(m,'A','B',.1)
    if kind == 'kinetic':
        arrhenius(m,'B','C',10); arrhenius(m,'C','B',2)
        if extra: arrhenius(m,'B','D',8); arrhenius(m,'D','B',8/3)
    elif kind == 'builtin': m.equilibrium('C',{'B':1},{'C':1},5.)
    else:
        if noop == 'first': ExternalFixture(m)
        ExternalFixture(m,constraints=[('C',{'B':5,'C':-1},0)])
        if extra: ExternalFixture(m,constraints=[('D',{'B':3,'D':-1},0)])
        if kind == 'conservative': ExternalFixture(m,constraints=[('B',{n:1 for n in m.ids},1.)])
        if noop == 'last': ExternalFixture(m)
    return m.finalize()


def trajectory(m,method=43,steps=200,dt=1.):
    with State(m) as s:
        s.conditions(298); s.set('A',1.)
        result = []
        for _ in range(steps):
            s.update(); s.solve(dt,method); result.append(s.y[:,0].copy())
        return np.array(result)


@group(PREFIX,'''KineticVsConservativeConstraintConvergence SimpleConstraintPreservesRatioNotMass
BuiltInVsExternalModelConstraintStepByStep MultiEquilibriumKineticVsComposedConstraints
ConstraintActiveWhenProcessOnlyModelAddedFirst''')
def compositions(id_):
    name = id_.split('.')[-1]
    if name == 'ConstraintActiveWhenProcessOnlyModelAddedFirst':
        with equilibrium_model('simple',noop='last') as a, equilibrium_model('simple',noop='first') as b:
            for m in (a,b): assert model_constraint_count(m.handle) == 1 and model_is_algebraic(m.handle,2)
            av,bv = trajectory(a,steps=20)[-1],trajectory(b,steps=20)[-1]
            assert av[1] > 0 and bv[1] > 0; near(av,bv,rtol=0)
            near(5*av[1]-av[2],0.,rtol=0,atol=1e-6)
        return
    if name == 'BuiltInVsExternalModelConstraintStepByStep':
        with equilibrium_model('builtin') as a, equilibrium_model('simple') as b:
            av,bv = trajectory(a,steps=100,dt=.5),trajectory(b,steps=100,dt=.5)
            near(av,bv,rtol=0,atol=1e-8)
            for v in (av,bv): near(5*v[:,1]-v[:,2],np.zeros(100),rtol=0,atol=1e-6)
        return
    extra = name == 'MultiEquilibriumKineticVsComposedConstraints'; simple = name == 'SimpleConstraintPreservesRatioNotMass'
    with equilibrium_model('kinetic',extra) as a, equilibrium_model('simple' if simple else 'conservative',extra) as b:
        av,bv = trajectory(a,3)[-1],trajectory(b)[-1]
        near(sum(av),1.,rtol=0,atol=1e-6)
        if simple:
            near(av[2]/av[1],5.,rtol=0,atol=1e-3); near(bv[2]/bv[1],5.,rtol=0,atol=1e-3)
            near(sum(bv),6.,rtol=0,atol=.1)
        else:
            expected = np.array([0,1,5,3] if extra else [0,1,5])/(9 if extra else 6)
            tolerance = np.array([1e-6]+[1e-3]*(len(expected)-1))
            near(av,expected,rtol=0,atol=tolerance); near(bv,expected,rtol=0,atol=tolerance)
            near(av,bv,rtol=0,atol=tolerance); near(sum(bv),1.,rtol=0,atol=1e-3)
            if extra:
                near(5*bv[1]-bv[2],0.,rtol=0,atol=1e-6); near(3*bv[1]-bv[3],0.,rtol=0,atol=1e-6)


@case(PREFIX+'TemperatureDependentConstraintParameter')
def temperature(_):
    with Model(list('ABC')) as m:
        arrhenius(m,'A','B',.1); m.parameter('C_K_eq')
        update = lambda t: 2*math.exp(3000*(1/298.15-1/t))
        ExternalFixture(m,constraints=[('C',{'B':'C_K_eq','C':-1},0)],updates={'C_K_eq':update}); m.finalize()
        for t in (298.15,350.):
            with State(m) as s:
                s.conditions(t); s.set('A',1.)
                for _ in range(200): s.update(); s.solve(1.,43)
                assert s.y[1,0] > 0
                near(s.y[2,0]/s.y[1,0],update(t),rtol=0,atol=1e-4)
                near(s.custom[0,0],update(t),rtol=0)
        assert update(350) > update(298.15)


@case(PREFIX+'BuilderCanBuildMultipleSolvers')
def reusable_builder(_):
    with Model(list('ABCD')) as original, ExitStack() as stack:
        arrhenius(original,'A','B',.1); original.equilibrium('D',{'A':1},{'D':1},2.)
        ExternalFixture(original,constraints=[('C',{'B':5,'C':-1},0)])
        results = []
        for _ in range(2):
            m = Model(); model_free(m.handle); m.handle = model_clone_builder(original.handle); assert m.handle
            stack.enter_context(m); m.finalize()
            assert model_constraint_count(m.handle) == 2
            with State(m) as s:
                s.y[:,0] = [1,0,0,2]; s.conditions(298)
                for _ in range(20): s.update(); s.solve(1.,43)
                results.append(s.y.copy())
        near(results[0],results[1],rtol=0)
        y = results[0][:,0]; near(2*y[0]-y[3],0.,rtol=0,atol=1e-6); near(5*y[1]-y[2],0.,rtol=0,atol=1e-6)


@group('external_model_constraints:ExternalModelFiniteDifferenceJacobian.',
       'ProcessForcingJacobian ConstraintResidualJacobian EquilibriumConstraintModelJacobian')
def external_fd(id_):
    name = id_.split('.')[-1]; process = name == 'ProcessForcingJacobian'; eq = name == 'EquilibriumConstraintModelJacobian'
    with Model(['A','B'] if eq else [GAS,AQ]) as m:
        e = ExternalFixture(m,constraints=[('B',{'A':2.5,'B':-1},0)]) if eq else transfer(m,3.5 if process else .1,None if process else 1.)
        m.finalize()
        with State(m,1 if eq else 2) as s:
            s.y[:] = [[3],[5]] if eq else [[.8,.3],[.2,.7]] if process else [[.6,.2],[.4,.8]]
            fn,jfn = (e.force_cb,e.jacobian_cb) if process else (e.constraint_cb,e.constraint_jac_cb)
            @ForcingFunction
            def callback(n,c,y,f,context): fn(s.handle,y,f,None)
            fd = np.zeros((2,2,s.cells)); require(finite_difference(2,s.cells,ptr(s.y),1e-8,callback,None,ptr(fd)))
            packed = np.zeros((model_jacobian_size(m.handle),s.cells)); jfn(s.handle,ptr(s.y),ptr(packed),None)
            negative = np.zeros_like(fd); pattern = np.zeros((2,2),np.uint8)
            for r,c in e.pairs:
                i,j = m.indices[r],m.indices[c]; pattern[i,j] = 1
                negative[i,j] = packed[model_jacobian_index(m.handle,i,j)]
            assert compare_jacobian(2,s.cells,pattern.ctypes.data_as(BP),ptr(negative),ptr(fd),1e-7,1e-7).passed
            if process or eq: assert check_sparsity(2,s.cells,pattern.ctypes.data_as(BP),ptr(fd),1e-6).passed
