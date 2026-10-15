"""C ownership, validation, and state indexing equivalents of the state cases.

Original assertions: Copyright (C) 2023-2026 University Corporation for Atmospheric Research.
SPDX-License-Identifier: Apache-2.0
"""
from api import *
from harness import case, equal, near

NAMES = ['foo','bar','baz','quz']


def fixture(labels=('quux','corge')):
    m = Model(NAMES)
    for label in labels:
        m.reaction([],[],rate(USER,label=label))
    for i in range(10-len(labels)):
        m.reaction([],[],rate())
    return m.finalize()


@case(*['state:State.'+name for name in ('DefaultConstructor','Constructor','CopyConstructor','CopyAssignmentOperator','MoveConstructor','MoveAssignmentOperator')])
def construction(identifier):
    name = identifier.split('.')[-1]
    if name == 'DefaultConstructor':
        states = [StateData() for _ in range(10)]
        for state in states:
            equal(state.cells,0)
        return
    m = fixture(); s = State(m,3); copied = None
    try:
        equal(s.data.cells,3); equal(m.data.nvariables,4); equal(m.data.nreactions,10); equal(m.data.ncustom,2)
        for c in range(3): equal((s.data.conditions[c].temperature,s.data.conditions[c].pressure,s.data.conditions[c].air_density),(0,0,0))
        for i,n in enumerate(NAMES):
            equal(m.var(n),i)
        for i,n in enumerate(('quux','corge')):
            equal(lib.chem_custom_id(m.p,n.encode()),i)
        s.tolerances(1e-5,[1e-10]*4)
        if name.startswith('Copy'):
            copied = State(m,pointer=lib.chem_state_copy(s.p)); t = copied
        else:
            t = s # C ownership transfer retains the allocated state address.
        equal(t.data.cells,3); equal(t.data.relative_tolerance,1e-5)
        for i in range(4):
            equal(t.data.absolute_tolerance[i],1e-10)
        for field,count in [('variables',12),('custom',6),('rates',30)]:
            equal(list(getattr(t.data,field)[:count]),list(getattr(s.data,field)[:count]))
        if copied:
            s.set('foo',[1,2,3]); equal(copied.get('foo'),0)
    finally:
        if copied: copied.close()
        s.close(); m.close()


@case(*['state:State.'+name for name in ('SettingSingleConcentrationWithInvalidArgumentsThowsException',
    'SetSingleConcentration','SetConcentrationByElementSingleValue','SetConcentrationByElementVector',
    'SettingConcentrationsWithInvalidArguementsThrowsException','SetConcentrations','SetStateWithSquareBracketOperator','SetMultiCellStateWithSquareBrackets')])
def concentrations(identifier):
    name = identifier.split('.')[-1]
    cells = 1 if name in ('SetConcentrationByElementSingleValue','SetStateWithSquareBracketOperator') else 3
    m = fixture(); s = State(m,cells)
    try:
        if name.startswith('Setting'):
            equal(lib.chem_state_set(s.p,b'foo',arr([1]),1),4)
            equal(lib.chem_state_set(s.p,b'foo',arr([1,2]),2),4)
            equal(lib.chem_state_set(s.p,b'not foo',arr([1]),1),3)
            equal(lib.chem_state_set(s.p,b'not foo',arr([1,2,3,4]),4),3)
            equal(lib.chem_state_set(s.p,b'FUU',arr([.1]),1),3)
        elif name == 'SetSingleConcentration':
            s.set('bar',[12,42,35.2])
            for c,v in enumerate([12,42,35.2]): equal(s.get('bar',c),v)
            t = State(m,1)
            try:
                t.set('bar',324.2); near(t.get('bar'),324.2,324.2*10*2.220446049250313e-16)
            finally: t.close()
        elif name == 'SetConcentrationByElementSingleValue':
            s.set('baz',42); s.set('quz',12); equal(s.get('baz'),42); equal(s.get('quz'),12)
        elif name == 'SetConcentrationByElementVector':
            s.set('foo',[12,42,35.2])
            for c,v in enumerate([12,42,35.2]): equal(s.get('foo',c),v)
        elif name == 'SetConcentrations':
            data = {'bar':[.2,.22,.222],'baz':[.3,.33,.333],'foo':[.9,.99,.999],'quz':[.4,.44,.444]}
            for n,values in data.items(): s.set(n,values)
            for c in range(3):
                for n in NAMES: equal(s.get(n,c),data[n][c])
        elif name == 'SetStateWithSquareBracketOperator':
            for n,v in zip(NAMES,[42,12,35.2,24.2]): s.set(n,v); equal(s.get(n),v)
            at = m.var('foo')*s.data.stride
            for op,expected in [(lambda x:x+8,50),(lambda x:x-20,30),(lambda x:x*2,60),(lambda x:x/3,20)]:
                s.data.variables[at] = op(s.data.variables[at]); equal(s.get('foo'),expected)
        else:
            for n,values in [('foo',[1,2,3]),('bar',[4,5,6])]:
                s.set(n,values)
                for c,v in enumerate(values): equal(s.get(n,c),v)
    finally: s.close(); m.close()


@case(*['state:State.'+name for name in ('SettingCustomRateParameterWithInvalidVectorSizeThrowsException',
    'SettingCustomRateParameterWithInvalidLabelThrowsException','SetCustomRateParameter','SetCustomRateParameters',
    'UnsafelySetCustomRateParameterOneCell','UnsafelySetCustomRateParameterMultiCell',
    'UnsafelySetCustomRateParameterCatchesTooFewGridCells','UnsafelySetCustomRateParameterCatchesTooParameters')])
def custom_parameters(identifier):
    name = identifier.split('.')[-1]
    labels = ['O1','O2','O3','AAA','BBB']
    if name in ('SettingCustomRateParameterWithInvalidLabelThrowsException','SetCustomRateParameter','UnsafelySetCustomRateParameterOneCell'): cells = 1
    elif 'Catches' in name: cells = 2
    else: cells = 3
    m = fixture(labels); s = State(m,cells)
    try:
        if 'InvalidVectorSize' in name:
            for n in labels: equal(lib.chem_state_set_custom(s.p,n.encode(),arr([.3]),1),4)
        elif 'InvalidLabel' in name:
            equal(lib.chem_state_set_custom(s.p,b'CCC',arr([.7]),1),3)
        elif name == 'SetCustomRateParameter':
            s.set('O2',42.3,custom=True); equal(s.data.custom[1],42.3)
        elif name == 'SetCustomRateParameters':
            data = {'O3':[.3,.33,.333],'O1':[.1,.11,.111],'O2':[.5,.55,.555],'BBB':[.7,.77,.777],'AAA':[.5,.55,.555]}
            for n,values in data.items(): s.set(n,values,custom=True)
            for c in range(cells):
                for i,n in enumerate(labels): equal(s.data.custom[i*s.data.stride+c],data[n][c])
        elif 'CatchesTooFewGridCells' in name:
            equal(lib.chem_state_custom_array(s.p,arr([.1,.2,.3,.4,.5]),1,5),4)
        elif 'CatchesTooParameters' in name:
            equal(lib.chem_state_custom_array(s.p,arr([.1,.2,.3]),1,3),4)
        else:
            values = [.1,.2,.3,.4,.5]
            equal(lib.chem_state_custom_array(s.p,arr(values*cells),cells,5),0)
            for c in range(cells):
                for i,v in enumerate(values): near(s.data.custom[i*s.data.stride+c],v,abs(v)*10*2.220446049250313e-16)
    finally: s.close(); m.close()


@case('solver_lambda_rate_constant:Solver.GetLambdaRateConstantByNameCanOverrideLambda')
def lambda_override(_):
    callback_type = C.CFUNCTYPE(D,C.POINTER(Conditions),DP,P)
    first = callback_type(lambda c,u,ctx: .001*c.contents.temperature)
    second = callback_type(lambda c,u,ctx: .002*c.contents.temperature)
    m = Model(['A','B']); m.reaction(['A'],[('B',1)],rate(LAMBDA,label='lambda_rc',function=C.cast(first,P))); m.finalize()
    s = State(m)
    try:
        r = lib.chem_model_lambda(m.p,b'lambda_rc'); assert r
        r.contents.function = C.cast(second,P)
        s.tolerances(1e-6,[1e-6]*2); s.set('A',1)
        s.data.conditions[0] = Conditions(300,101253.3,1e6); s.rates()
        near(s.data.rates[0],.6,1e-12)
        p = lib.chem_euler_parameters(); stats = Stats()
        equal(lib.chem_backward_euler(C.byref(p),1,s.p,C.byref(stats)),0); equal(stats.final_time,1)
    finally: s.close(); m.close()


@case('backward_euler:BackwardEuler.CanCallSolve')
def euler_call(_):
    m = Model(['A','B','C'])
    m.reaction(['A'],[('B',1)],rate(A=2.15e-11,C=110)); m.reaction(['B'],[('C',1)],rate(A=3.3e-11,C=55)); m.finalize()
    s = State(m)
    try:
        s.tolerances(1e-6,[1e-6]*3); s.set('A',1); s.data.conditions[0] = Conditions(272.5,101253.3,1e6); s.rates()
        p = lib.chem_euler_parameters(); stats = Stats()
        equal(lib.chem_backward_euler(C.byref(p),1,s.p,C.byref(stats)),0)
    finally: s.close(); m.close()


@case('backward_euler:BackwardEuler.IsConverged')
def euler_converged(_):
    for cells in (4,5):
        m = Model(['A','B','C']).finalize(); s = State(m,cells)
        try:
            s.tolerances(1e-3,[1e-6]*3)
            residual = arr([0]*(cells*3)); y = arr([0]*(cells*3))
            equal(lib.chem_euler_converged(s.p,1e-6,residual,y),1)
            residual[1*cells] = 1e-5; equal(lib.chem_euler_converged(s.p,1e-6,residual,y),0)
            equal(lib.chem_euler_converged(s.p,1e-4,residual,y),1)
            residual[2*cells+3] = 1e-3; equal(lib.chem_euler_converged(s.p,1e-4,residual,y),0)
            y[2*cells+3] = 10; equal(lib.chem_euler_converged(s.p,1e-4,residual,y),1)
            residual[2*cells+3] = .1; equal(lib.chem_euler_converged(s.p,1e-4,residual,y),0)
            s.data.absolute_tolerance[2] = 1; equal(lib.chem_euler_converged(s.p,1e-4,residual,y),1)
        finally: s.close(); m.close()
