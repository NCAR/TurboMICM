"""External constraint composition and process/constraint Jacobian filtering."""
import math
from api import *
from harness import case, equal, near, real_equal
from external_helpers import affine_constraint, custom, external, transfers
from dae_integration_cases import equilibrium
from constraint_cases import jac_value, verify_jacobian

AQ = 'AEROSOL.A_AQ'


def stub(m,k,total=None,solvent=False):
    if not solvent: transfers(m,[('A_GAS',AQ,float(k))])
    else:
        gi,ai,si = (m.var(n) for n in ['A_GAS',AQ,'S'])
        def forcing(mp,sp,y,f,ctx):
            s = sp.contents; st = s.stride
            for c in range(s.cells):
                value = k*y[gi*st+c]*y[si*st+c]; f[gi*st+c] -= value; f[ai*st+c] += value
        def jacobian(mp,sp,y,j,ctx):
            s = sp.contents; st,n = s.stride,mp.contents.nvariables
            for c in range(s.cells):
                for col,value in [(gi,k*y[si*st+c]),(si,k*y[gi*st+c])]:
                    j[mp.contents.index[gi*n+col]*st+c] += value; j[mp.contents.index[ai*n+col]*st+c] -= value
        external(m,forcing,jacobian)
    if total is not None: affine_constraint(m,AQ,[('A_GAS',1),(AQ,1)],total)


METADATA = ('AddExternalModelWithConstraints','AddExternalModelProcessOnly','AddExternalModelOnlyStandardRosenbrock',
            'CombinedBuiltInAndExternalConstraints','AddExternalModelConstraintsOnly','DAESolveEnforcesConservation',
            'AlgebraicRowHoldsOnlyConstraintTerms','MultiGridCell','ProcessJacobianElementInAlgebraicRowSurvivesFiltering',
            'ProcessJacobianElementInBuiltInAlgebraicRowSurvivesFiltering')


@case(*['external_model_constraints:ExternalModelConstraints.'+n for n in METADATA])
def stub_cases(identifier):
    name = identifier.split('.')[-1]; solvent = name == 'ProcessJacobianElementInAlgebraicRowSurvivesFiltering'
    builtin = name == 'ProcessJacobianElementInBuiltInAlgebraicRowSurvivesFiltering'
    combined = name == 'CombinedBuiltInAndExternalConstraints'
    process_only = name in ('AddExternalModelProcessOnly','AddExternalModelOnlyStandardRosenbrock') or builtin
    m = Model(['A_GAS']+(['S'] if solvent else ['B','C'] if combined else ['B'] if builtin else [])+[AQ])
    if combined:
        m.reaction(['A_GAS'],[('B',1)],rate(A=.05)); equilibrium(m,'B','C',5,'B_C_eq')
    if builtin: equilibrium(m,'B',AQ,5,'B_AQ_eq')
    k = 10 if name == 'AlgebraicRowHoldsOnlyConstraintTerms' else .1 if name in ('DAESolveEnforcesConservation','MultiGridCell') or solvent or builtin else .01
    stub(m,k,None if process_only else 1,solvent); m.finalize(); cells = 3 if name == 'MultiGridCell' else 1; s = State(m,cells)
    try:
        equal(m.data.nvariables,4 if combined else 3 if solvent or builtin else 2)
        equal(m.data.nconstraints,2 if combined else 0 if process_only and not builtin else 1)
        equal(m.data.algebraic[m.var(AQ)],int(not process_only or builtin)); equal(m.data.algebraic[m.var('A_GAS')],0)
        if name == 'AlgebraicRowHoldsOnlyConstraintTerms':
            s.set('A_GAS',.7); s.set(AQ,.2); f = arr([0]*2); lib.chem_forcing(s.p,s.data.variables,f)
            lib.chem_jacobian(s.p,s.data.variables,s.data.jacobian); gi,ai = m.var('A_GAS'),m.var(AQ)
            near(f[gi],-k*.7,1e-12); near(f[ai],.7+.2-1,1e-12)
            for r,c,expected in [(gi,gi,k),(ai,gi,-1),(ai,ai,-1)]: near(jac_value(m,s,s.data.jacobian,r,c),expected,1e-12)
        elif builtin:
            assert m.data.index[m.var(AQ)*m.data.nvariables+m.var('A_GAS')] >= 0
        elif name in ('DAESolveEnforcesConservation','MultiGridCell','AddExternalModelConstraintsOnly') or solvent:
            for c in range(cells):
                gas = .8 if name == 'AddExternalModelConstraintsOnly' else .9-.2*c
                s.data.variables[m.var('A_GAS')*cells+c] = gas; s.data.variables[m.var(AQ)*cells+c] = 1-gas
                s.data.conditions[c] = Conditions(298,101325,0)
            if solvent: s.set('S',1)
            for _ in range(20 if name == 'DAESolveEnforcesConservation' else 10 if solvent else 1):
                s.solve(10 if name == 'DAESolveEnforcesConservation' or solvent else 50,5)
                for c in range(cells): near(s.get('A_GAS',c)+s.get(AQ,c),1,1e-4)
    finally: s.close(); m.close()


def abc(kind='simple',multi=False,builtin=False,no_op_first=None,four=False):
    names = list('ABCD' if multi or four else 'ABC'); m = Model(names); m.reaction(['A'],[('B',1)],rate(A=.1))
    if kind == 'kinetic':
        for r,p,k in [('B','C',10),('C','B',2)]+([('B','D',8),('D','B',8/3)] if multi else []): m.reaction([r],[(p,1)],rate(A=k))
    else:
        if no_op_first: external(m)
        if builtin: equilibrium(m,'B','C',5,'B_C_eq')
        else: affine_constraint(m,'C',[('B',5),('C',-1)])
        if kind == 'conservative': affine_constraint(m,'B',[(n,1) for n in names],1)
        if multi: affine_constraint(m,'D',[('B',3),('D',-1)])
        if four: equilibrium(m,'A','D',2,'A_D_eq')
        if no_op_first is False: external(m)
    m.finalize(); return m


def initial_state(m):
    s = State(m); s.set('A',1)
    s.data.conditions[0] = Conditions(298,101325,0); return s


def run_abc(kind,multi=False):
    m = abc(kind,multi); s = initial_state(m)
    try:
        for _ in range(200): s.rates(); s.solve(1,3 if kind == 'kinetic' else 5)
        return [s.get(n) for n in ('ABCD' if multi else 'ABC')]
    finally: s.close(); m.close()


@case('external_model_constraints:ExternalModelConstraints.KineticVsConservativeConstraintConvergence',
      'external_model_constraints:ExternalModelConstraints.SimpleConstraintPreservesRatioNotMass',
      'external_model_constraints:ExternalModelConstraints.MultiEquilibriumKineticVsComposedConstraints')
def equilibrium_convergence(identifier):
    simple = identifier.endswith('SimpleConstraintPreservesRatioNotMass'); multi = identifier.endswith('MultiEquilibriumKineticVsComposedConstraints')
    kinetic = run_abc('kinetic',multi); constrained = run_abc('simple' if simple else 'conservative',multi)
    near(sum(kinetic),1,1e-6)
    if simple:
        near(kinetic[2]/kinetic[1],5,1e-3); near(constrained[2]/constrained[1],5,1e-3); near(sum(constrained),6,.1)
    else:
        denominator = 9 if multi else 6; expected = [0,1/denominator,5/denominator]+([3/denominator] if multi else [])
        for a,b,e in zip(kinetic,constrained,expected):
            near(a,e,1e-6 if e == 0 else 1e-3); near(b,e,1e-6 if e == 0 else 1e-3); near(a,b,1e-6 if e == 0 else 1e-3)
        near(sum(constrained),1,1e-3)
        if multi: near(5*constrained[1]-constrained[2],0,1e-6); near(3*constrained[1]-constrained[3],0,1e-6)


@case('external_model_constraints:ExternalModelConstraints.BuiltInVsExternalModelConstraintStepByStep',
      'external_model_constraints:ExternalModelConstraints.ConstraintActiveWhenProcessOnlyModelAddedFirst',
      'external_model_constraints:ExternalModelConstraints.BuilderCanBuildMultipleSolvers')
def equivalent_models(identifier):
    name = identifier.split('.')[-1]; multiple = name == 'BuilderCanBuildMultipleSolvers'; order = name == 'ConstraintActiveWhenProcessOnlyModelAddedFirst'
    models = [abc(four=True)] if multiple else [abc(no_op_first=False),abc(no_op_first=True)] if order else [abc(builtin=True),abc()]
    states = [initial_state(models[0]),initial_state(models[-1])]
    try:
        if multiple:
            for s in states: s.set('D',2); equal(s.model.data.nconstraints,2)
        for _ in range(20 if multiple or order else 100):
            for s in states: s.rates(); s.solve(1 if multiple or order else .5,5)
            values = [[s.get(n) for n in ('ABCD' if multiple else 'ABC')] for s in states]
            if multiple or order: equal(*values)
            else:
                for a,b in zip(*values): near(a,b,1e-8)
            for a,b,c,*rest in values:
                near(5*b-c,0,1e-6)
                if multiple: near(2*a-rest[0],0,1e-6)
        if order:
            for s in states: assert s.get('B') > 0
    finally:
        for s in states: s.close()
        for m in models: m.close()


@case('external_model_constraints:ExternalModelConstraints.TemperatureDependentConstraintParameter')
def temperature_constraint(_):
    m = Model(list('ABC')); m.reaction(['A'],[('B',1)],rate(A=.1)); ki = custom(m,'C_K_eq')
    affine_constraint(m,'C',[('B',0),('C',-1)],parameter_index=ki)
    def update(sp,context):
        s = sp.contents
        for c in range(s.cells): s.custom[ki*s.stride+c] = 2*math.exp(3000*(1/298.15-1/s.conditions[c].temperature))
    external(m,update=update); m.finalize()
    try:
        for t in (298.15,350):
            s = initial_state(m)
            try:
                s.data.conditions[0].temperature = t
                for _ in range(200): s.rates(); s.solve(1,5)
                expected = 2*math.exp(3000*(1/298.15-1/t))
                assert s.get('B') > 0; near(s.get('C')/s.get('B'),expected,1e-4)
                if t == 350: assert expected > 2
            finally: s.close()
    finally: m.close()


@case(*['external_model_constraints:ExternalModelFiniteDifferenceJacobian.'+n for n in ('ProcessForcingJacobian','ConstraintResidualJacobian','EquilibriumConstraintModelJacobian')])
def external_fd(identifier):
    name = identifier.split('.')[-1]
    if name == 'EquilibriumConstraintModelJacobian':
        m = Model(['A','B']); affine_constraint(m,'B',[('A',2.5),('B',-1)]); rows = [[3,5]]
    else:
        m = Model(['A_GAS',AQ]); rows = [[.8,.2],[.3,.7]] if name == 'ProcessForcingJacobian' else [[.6,.4],[.2,.8]]
        if name == 'ProcessForcingJacobian': stub(m,3.5)
        else: affine_constraint(m,AQ,[('A_GAS',1),(AQ,1)],1)
    m.finalize(); s = State(m,len(rows))
    try:
        for c,row in enumerate(rows):
            for j,x in enumerate(row): s.data.variables[j*len(rows)+c] = x
        verify_jacobian(s,constraints_only=name != 'ProcessForcingJacobian')
    finally: s.close(); m.close()
