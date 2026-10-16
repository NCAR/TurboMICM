"""Species, phase, system, and state ownership/lookup equivalents in C."""
from ffi import *
from harness import case,group


def property_set(s,key,value,kind=0):
    data=C.c_char_p(value.encode()) if kind==1 else I(value) if kind in (2,3) else D(value)
    require(species_set_property(s,key.encode(),kind,data if kind==1 else C.byref(data)))


def property_get(s,key,kind=0):
    data=C.c_char_p() if kind==1 else I() if kind in (2,3) else D()
    require(species_get_property(s,key.encode(),kind,C.byref(data)))
    return data.value.decode() if kind==1 else data.value


@group('species:Species.','StringConstructor StringAndVectorConstructor GetProperty SetThirdBody CopyConstructor CopyAssignment')
def species(id_):
    s=species_create(b'thing');assert s
    try:
        assert species_name(s)==b'thing'
        assert not species_is_parameterized(s)
        require(species_parameterize(s,1,0,100,0,0))
        assert species_is_parameterized(s)
        assert species_parameter_value(s,C.byref(Conditions(12.5,0,0)))==1250
        require(species_parameterize(s,0,0,0,0,0))
        property_set(s,'name [units]',1)
        property_set(s,'name2 [units2]',2)
        assert species_property_count(s,0)==2
        assert property_get(s,'name [units]')==1 and property_get(s,'name2 [units2]')==2
        for kind in range(4):
            data=C.c_char_p() if kind==1 else D()
            assert species_get_property(s,b'not there',kind,C.byref(data))==-3
        assert species_get_property(s,b'name [units]',100,C.byref(D()))==-1
        property_set(s,'foo','bar',1);property_set(s,'baz',42,2);property_set(s,'qux',1,3)
        for enabled in (0,1):
            require(species_parameterize(s,enabled,15.4,0,0,0))
            copy=species_clone(s);assert copy
            try:
                assert species_name(copy)==b'thing'
                assert species_property_count(copy,0)==2
                assert property_get(copy,'name [units]')==1 and property_get(copy,'name2 [units2]')==2
                assert property_get(copy,'foo',1)=='bar' and property_get(copy,'baz',2)==42 and property_get(copy,'qux',3)==1
                assert species_is_parameterized(copy)==enabled
                if enabled:assert species_parameter_value(copy,C.byref(Conditions()))==15.4
                property_set(s,'foo','changed',1)
                assert property_get(copy,'foo',1)=='bar'
                property_set(s,'foo','bar',1)
            finally:species_free(copy)
    finally:species_free(s)
    s=species_create(b'M')
    try:
        require(species_parameterize(s,1,0,0,0,1))
        assert species_name(s)==b'M' and species_is_parameterized(s)
        assert species_parameter_value(s,C.byref(Conditions(0,0,42.4)))==42.4
    finally:species_free(s)


@group('phase:PhaseSpecies.','ConstructsWithoutDiffusionCoefficient ConstructsWithDiffusionCoefficient SetDiffusionCoefficient')
def phase_species(_):
    s=species_create(b'CO2');assert s
    try:
        assert species_name(s)==b'CO2'
        assert species_get_property(s,b'diffusion coefficient [m2 s-1]',0,C.byref(D()))==-3
        for value in (1.23e-5,2.5e-6):
            property_set(s,'diffusion coefficient [m2 s-1]',value)
            assert property_get(s,'diffusion coefficient [m2 s-1]')==value
    finally:species_free(s)


@group('phase:Phase.','Constructor ConstructorWithParameterizedSpecies UniqueNamesWithParameterizedSpecies')
@group('system:System.','ConstructorWithAllParameters ConstructorWithParameterizedSpecies')
def phase_system(_):
    for parameterized in (False,True):
        p=phase_create(b'gas' if parameterized else b'test_phase');assert p
        try:
            names=('foo','bar','baz') if parameterized else ('species1','species2')
            for name in names:
                s=species_create(name.encode());assert s
                if parameterized and name=='bar':require(species_parameterize(s,1,42,0,0,0))
                require(phase_add(p,s));species_free(s)
            assert phase_species_count(p)==len(names) and phase_state_size(p)==2
            assert phase_name(p)==(b'gas' if parameterized else b'test_phase')
            expected=['foo','baz'] if parameterized else list(names)
            for i,name in enumerate(expected):
                assert phase_state_name(p,i,0).decode()==name
                assert phase_state_name(p,i,1).decode()==phase_name(p).decode()+'.'+name
            with Model() as m:
                require(model_add_phase(m.handle,p,1));m.finalize()
                assert model_builder_species_count(m.handle)==len(names)
                assert m.names==expected
                reordered=[model_species_name(m.handle,i).decode() for i in (1,0)]
                assert reordered==expected[::-1]
        finally:phase_free(p)
    # The system-specific fixture uses a third, constant species named paramSpecies.
    with Model(['foo','bar']) as m:
        s=species_create(b'paramSpecies')
        require(species_parameterize(s,1,64.2,0,0,0))
        require(model_add_species_spec(m.handle,s));species_free(s)
        m.finalize();assert m.names==['foo','bar'] and model_builder_species_count(m.handle)==3


def state_model(params=('quux','corge')):
    m=Model(['foo','bar','baz','quz'])
    for p in params:m.parameter(p)
    for i in range(10):m.reaction(['foo'],{'bar':1},rate_default(0))
    return m.finalize()


@case('state:State.DefaultConstructor')
def default_state(_):
    # A zero-initialized C owning handle is the empty/default state.
    handles=(P*10)()
    for h in handles:
        assert state_cells(h)==0 and not state_values(h)
        state_free(h)


@group('state:State.','Constructor CopyConstructor CopyAssignmentOperator MoveConstructor MoveAssignmentOperator')
def ownership(_):
    with state_model() as m,State(m,3) as s:
        assert s.y.shape==(4,3) and s.custom.shape==(2,3) and s.rates.shape==(10,3)
        assert m.indices==dict(zip(('foo','bar','baz','quz'),range(4)))
        assert [model_parameter_index(m.handle,n.encode()) for n in ('quux','corge')]==[0,1]
        s.y[:]=np.arange(12).reshape(4,3);s.custom[:]=np.arange(6).reshape(2,3)
        s.rates[:]=np.arange(30).reshape(10,3);s.conditions();s.tolerance(1e-5,1e-10)
        for copy in (state_clone(s.handle),state_create(m.handle,3)):
            assert copy
            try:
                require(state_copy(copy,s.handle));require(state_copy(copy,copy))
                near(np.ctypeslib.as_array(state_values(copy),shape=(12,)),s.y.ravel(),rtol=0)
                near(np.ctypeslib.as_array(state_parameters(copy),shape=(6,)),s.custom.ravel(),rtol=0)
                near(np.ctypeslib.as_array(state_rates(copy),shape=(30,)),s.rates.ravel(),rtol=0)
                assert state_relative_tolerance(copy)==1e-5
                near(np.ctypeslib.as_array(state_absolute_tolerances(copy),shape=(4,)),[1e-10]*4,rtol=0)
                state_values(copy)[0]=-999;assert s.y[0,0]==0
                source=P(copy);dest=P(state_create(m.handle,1));state_move(C.byref(dest),C.byref(source))
                assert not source and dest.value==copy
                assert state_cells(dest)==3
            finally:state_free(copy)


@group('state:State.','SettingSingleConcentrationWithInvalidArgumentsThowsException SettingConcentrationsWithInvalidArguementsThrowsException')
def invalid_concentration(_):
    with state_model() as m,State(m,3) as s:
        for name,values in (('foo',[1]),('foo',[1,2]),('not foo',[1]),('not foo',[1,2,3,4]),('FUU',[.1])):
            assert state_set_exact(s.handle,name.encode(),doubles(values),len(values))<0
        near(s.y,np.zeros((4,3)),rtol=0)


@group('state:State.','SetSingleConcentration SetConcentrationByElementSingleValue SetConcentrationByElementVector SetConcentrations SetStateWithSquareBracketOperator SetMultiCellStateWithSquareBrackets')
def concentrations(_):
    with state_model() as m:
        with State(m,3) as s:
            s.set('bar',[12,42,35.2]);near(s.y[1],[12,42,35.2],rtol=0)
            s.set('foo',[12,42,35.2]);near(s.y[0],[12,42,35.2],rtol=0)
            values={'bar':[.2,.22,.222],'baz':[.3,.33,.333],'foo':[.9,.99,.999],'quz':[.4,.44,.444]}
            for name,v in values.items():require(state_set_exact(s.handle,name.encode(),doubles(v),3))
            near(s.y.T,[[.9,.2,.3,.4],[.99,.22,.33,.44],[.999,.222,.333,.444]],rtol=0)
            s.set('foo',[1,2,3]);s.set('bar',[4,5,6])
            near(s.y[:2],[[1,2,3],[4,5,6]],rtol=0)
        with State(m) as s:
            s.set('bar',324.2);assert s.y[1,0]==324.2
            s.set('baz',42);s.set('quz',12);near(s.y[2:,0],[42,12],rtol=0)
            for name,v in zip(m.names,(42,12,35.2,24.2)):s.set(name,v)
            near(s.y[:,0],[42,12,35.2,24.2],rtol=0)
            value=state_values(s.handle)
            value[0]+=8;assert value[0]==50
            value[0]-=20;assert value[0]==30
            value[0]*=2;assert value[0]==60
            value[0]/=3;assert value[0]==20


@group('state:State.','SettingCustomRateParameterWithInvalidVectorSizeThrowsException SettingCustomRateParameterWithInvalidLabelThrowsException')
def invalid_parameters(_):
    with state_model(('O1','O2','O3','AAA','BBB')) as m,State(m,3) as s:
        for name in ('O1','O2','O3','AAA','BBB'):
            assert state_set_parameter_exact(s.handle,name.encode(),doubles([.3]),1)==-4
        assert state_set_parameter_exact(s.handle,b'CCC',doubles([.7]*3),3)==-3


@group('state:State.','SetCustomRateParameter SetCustomRateParameters UnsafelySetCustomRateParameterOneCell UnsafelySetCustomRateParameterMultiCell UnsafelySetCustomRateParameterCatchesTooFewGridCells UnsafelySetCustomRateParameterCatchesTooParameters')
def parameters(_):
    with state_model(('O1','O2','O3','AAA','BBB')) as m:
        for cells in (1,2,3):
            with State(m,cells) as s:
                if cells==1:s.parameter('O2',42.3);assert s.custom[1,0]==42.3
                if cells==3:
                    for name,v in {'O3':[.3,.33,.333],'O1':[.1,.11,.111],'O2':[.5,.55,.555],'BBB':[.7,.77,.777],'AAA':[.5,.55,.555]}.items():
                        require(state_set_parameter_exact(s.handle,name.encode(),doubles(v),3))
                    near(s.custom.T,[[.1,.5,.3,.5,.7],[.11,.55,.33,.55,.77],[.111,.555,.333,.555,.777]],rtol=0)
                data=np.tile([.1,.2,.3,.4,.5],(cells,1))
                require(state_set_parameter_matrix(s.handle,ptr(data),cells,5))
                near(s.custom.T,data,rtol=0)
                assert state_set_parameter_matrix(s.handle,ptr(data),cells-1,5)==-4
                assert state_set_parameter_matrix(s.handle,ptr(data),cells,3)==-4
