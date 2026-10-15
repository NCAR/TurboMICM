"""State metadata, named updates, copy, and ownership-transfer case ports."""
from functools import partial
from chemffi import *


class NamedValues(C.Structure):
    _fields_ = [('name',C.c_char_p),('count',Z),('values',DP)]


bind('chem_cell_count',Z,P)
bind('chem_state_clone',I,P,C.POINTER(P))
bind('chem_state_assign',I,C.POINTER(P),P)
bind('chem_state_move',I,C.POINTER(P),C.POINTER(P))
bind('chem_state_mechanism',P,P)
bind('chem_state_relative_tolerance',D,P)
for kind in ['species','parameter']:
    bind('chem_'+kind+'_index',Z,P,C.c_char_p)
    bind('chem_'+kind+'_column',DP,P,C.c_char_p)
bind('chem_const_species_column',DP,P,C.c_char_p)
for kind in ['concentrations','parameters']:
    bind('chem_set_'+kind,I,P,Z,C.POINTER(NamedValues))
bind('chem_set_parameter_matrix',I,P,Z,Z,DP)
NAMES = ['foo','bar','baz','quz']
PARAMETERS = ['O1','O2','O3','AAA','BBB']


def state(cells=3, parameters=PARAMETERS, tolerances=False):
    m = Mechanism(NAMES, [([],[],rate(0)) for _ in range(10)], parameters=parameters,
                  cells=cells,atols=[1e-10]*4 if tolerances else None)
    if tolerances: lib.chem_relative_tolerance(m.state,1e-5)
    return m


def named(target,values,parameters=False,status=0):
    arrays = [(D*len(v))(*v) for v in values.values()]
    data = (NamedValues*len(values))(*[NamedValues(k.encode(),len(a),a) for k,a in zip(values,arrays)])
    fn = lib.chem_set_parameters if parameters else lib.chem_set_concentrations
    assert fn(target,len(values),data) == status


def assert_shape(m,handle,parameters=PARAMETERS):
    assert lib.chem_cell_count(handle) == m.nc
    owner = lib.chem_state_mechanism(handle)
    d = lib.chem_mechanism_data(owner).contents
    assert (d.species_count,d.reaction_count,d.parameter_count) == (4,10,len(parameters))
    assert [d.species_names[i].decode() for i in range(4)] == NAMES
    for i,name in enumerate(NAMES): assert lib.chem_species_index(owner,name.encode()) == i
    for i,name in enumerate(parameters): assert lib.chem_parameter_index(owner,name.encode()) == i
    assert lib.chem_species_index(owner,b'missing') == Z(-1).value
    assert lib.chem_parameter_index(owner,b'missing') == Z(-1).value


def default():
    for _ in range(11):
        empty = P()
        assert lib.chem_state_create(None,0,C.byref(empty)) == 0
        assert lib.chem_cell_count(empty) == 0
        assert not lib.chem_state_mechanism(empty)
        assert lib.chem_update_rates(empty) == 1
        parameters,stats = lib.chem_rosenbrock_parameters(3),Stats()
        assert lib.chem_solve(empty,C.byref(parameters),1.,C.byref(stats)) == 1
        parameters = lib.chem_backward_euler_parameters()
        assert lib.chem_solve_backward_euler(empty,C.byref(parameters),1.,C.byref(stats)) == 1
        lib.chem_state_free(empty)


def constructor():
    with state(parameters=['quux','corge']) as m:
        assert_shape(m,m.state,['quux','corge'])


def copy_or_move(move=False,assignment=False):
    with state(parameters=['quux','corge'],tolerances=True) as m:
        target = P()
        original = m.state.value
        # Retain the source assertions on zero-initialized arrays and add
        # nonzero data to make accidental pointer sharing observable.
        for fn,count in [(lib.chem_concentrations,12),(lib.chem_parameters,6),(lib.chem_rate_constants,30)]:
            data = C.cast(fn(m.state),DP)
            assert list(data[:count]) == [0]*count
            for i in range(count): data[i] = i+.25
        if assignment: assert lib.chem_state_create(None,0,C.byref(target)) == 0
        if move: assert lib.chem_state_move(C.byref(target),C.byref(m.state)) == 0
        elif assignment: assert lib.chem_state_assign(C.byref(target),m.state) == 0
        else: assert lib.chem_state_clone(m.state,C.byref(target)) == 0
        try:
            assert_shape(m,target,['quux','corge'])
            assert lib.chem_state_relative_tolerance(target) == 1e-5
            assert list(lib.chem_absolute_tolerances(target)[:4]) == [1e-10]*4
            for fn,count in [(lib.chem_concentrations,12),(lib.chem_parameters,6),(lib.chem_rate_constants,30)]:
                assert list(fn(target)[:count]) == [i+.25 for i in range(count)]
            if move:
                assert not m.state.value and target.value == original
            else:
                assert target.value != original
                lib.chem_concentrations(target)[0] = 42
                assert lib.chem_concentrations(m.state)[0] == .25
        finally: lib.chem_state_free(target)


def invalid_concentration(single=True):
    with state() as m:
        if single:
            for name,values in [('foo',[1.]),('foo',[1.,2.]),('not foo',[1.]),('not foo',[1.,2.,3.,4.])]:
                named(m.state,{name:values},status=1)
        else:
            named(m.state,{'FUU':[.1],'bar':[.2],'baz':[.3],'quz':[.4]},status=1)
        assert m.rows() == [[0]*4 for _ in range(3)]


def single_concentration(by_element=False,vector=False):
    if not by_element or vector:
        with state() as m:
            name = 'foo' if by_element else 'bar'
            named(m.state,{name:[12.,42.,35.2]})
            column = lib.chem_species_column(m.state,name.encode())
            assert list(column[:3]) == [12.,42.,35.2]
            assert [row[NAMES.index(name)] for row in m.rows()] == [12.,42.,35.2]
    if not vector:
        with state(cells=1) as m:
            data = {'baz':[42.],'quz':[12.]} if by_element else {'bar':[324.2]}
            named(m.state,data)
            for k,v in data.items(): assert lib.chem_species_column(m.state,k.encode())[0] == v[0]


def concentrations():
    with state() as m:
        named(m.state,{'bar':[.2,.22,.222],'baz':[.3,.33,.333], 'foo':[.9,.99,.999],'quz':[.4,.44,.444]})
        assert m.rows() == [[.9,.2,.3,.4],[.99,.22,.33,.44],[.999,.222,.333,.444]]


def alias(single=True):
    with state(cells=1 if single else 3) as m:
        if single:
            named(m.state,{'foo':[42.],'bar':[12.],'baz':[35.2],'quz':[24.2]})
            for name,v in [('foo',42.),('bar',12.),('baz',35.2),('quz',24.2)]:
                assert lib.chem_species_column(m.state,name.encode())[0] == v
                assert lib.chem_const_species_column(m.state,name.encode())[0] == v
            ptr = lib.chem_species_column(m.state,b'foo')
            for op,expected in [(lambda x:x+8,50.),(lambda x:x-20,30.),(lambda x:x*2,60.),(lambda x:x/3,20.)]:
                ptr[0] = op(ptr[0]);assert ptr[0] == expected
        else:
            named(m.state,{'foo':[1.,2.,3.],'bar':[4.,5.,6.]})
            for name,expected in [('foo',[1.,2.,3.]),('bar',[4.,5.,6.])]:
                assert list(lib.chem_species_column(m.state,name.encode())[:3]) == expected
                assert list(lib.chem_const_species_column(m.state,name.encode())[:3]) == expected
                assert [row[NAMES.index(name)] for row in m.rows()] == expected


def invalid_parameter(size=True):
    with state(cells=3 if size else 1) as m:
        data = {'O3':[.3],'O1':[.1],'O2':[.5],('BBB' if size else 'CCC'):[.7],'AAA':[.5]}
        named(m.state,data,parameters=True,status=1)
        assert list(m.custom[:m.nc*5]) == [0]*(m.nc*5)


def parameter(single=True):
    with state(cells=1 if single else 3) as m:
        if single:
            named(m.state,{'O2':[42.3]},parameters=True)
            assert m.custom[1] == 42.3
        else:
            data = {'O3':[.3,.33,.333],'O1':[.1,.11,.111],'O2':[.5,.55,.555],'BBB':[.7,.77,.777],'AAA':[.5,.55,.555]}
            named(m.state,data,parameters=True)
            expected = [[.1,.5,.3,.5,.7],[.11,.55,.33,.55,.77],[.111,.555,.333,.555,.777]]
            assert [[m.custom[j*3+i] for j in range(5)] for i in range(3)] == expected


def parameter_matrix(cells,rows=None,columns=5,status=0):
    rows = cells if rows is None else rows
    with state(cells=cells) as m:
        data = (D*(rows*columns))(*[.1*(j+1) for _ in range(rows) for j in range(columns)])
        assert lib.chem_set_parameter_matrix(m.state,rows,columns,data) == status
        if not status:
            for i in range(cells):
                for j in range(columns): assert abs(m.custom[j*cells+i] - [.1,.2,.3,.4,.5][j]) <= abs(data[j])*10*2.220446049250313e-16


CASES = {
    'State.DefaultConstructor':default,
    'State.Constructor':constructor,
    'State.CopyConstructor':partial(copy_or_move),
    'State.CopyAssignmentOperator':partial(copy_or_move,assignment=True),
    'State.MoveConstructor':partial(copy_or_move,move=True),
    'State.MoveAssignmentOperator':partial(copy_or_move,move=True,assignment=True),
    'State.SettingSingleConcentrationWithInvalidArgumentsThowsException':invalid_concentration,
    'State.SetSingleConcentration':single_concentration,
    'State.SetConcentrationByElementSingleValue':partial(single_concentration,True),
    'State.SetConcentrationByElementVector':partial(single_concentration,True,True),
    'State.SettingConcentrationsWithInvalidArguementsThrowsException':partial(invalid_concentration,False),
    'State.SetConcentrations':concentrations,
    'State.SetStateWithSquareBracketOperator':alias,
    'State.SetMultiCellStateWithSquareBrackets':partial(alias,False),
    'State.SettingCustomRateParameterWithInvalidVectorSizeThrowsException':invalid_parameter,
    'State.SettingCustomRateParameterWithInvalidLabelThrowsException':partial(invalid_parameter,False),
    'State.SetCustomRateParameter':parameter,
    'State.SetCustomRateParameters':partial(parameter,False),
    'State.UnsafelySetCustomRateParameterOneCell':partial(parameter_matrix,1),
    'State.UnsafelySetCustomRateParameterMultiCell':partial(parameter_matrix,3),
    'State.UnsafelySetCustomRateParameterCatchesTooFewGridCells':partial(parameter_matrix,2,1,status=1),
    'State.UnsafelySetCustomRateParameterCatchesTooParameters':partial(parameter_matrix,2,1,3,status=1),
}
