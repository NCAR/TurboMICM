"""Species, phase, gas-system, and version case assertions."""
from functools import partial
from systemffi import *

BASIC = [('name [units]',0,1.),('name2 [units2]',0,2.)]


def string_constructor():
    with Species('thing') as s:
        assert s.name == 'thing'
        assert not s.parameterization.has_value
        s.parameterization = Parameterization(temperature=100,has_value=1)
        assert s.parameterization.has_value
        p = s.parameterization
        assert lib.chem_parameterization_evaluate(C.byref(p),12.5,0,0) == 1250.


def properties(errors=False):
    with Species('thing',BASIC) as s:
        assert s.name == 'thing'
        assert lib.chem_species_property_count(s.handle,0) == 2
        assert s.get('name [units]') == 1.
        assert s.get('name2 [units2]') == 2.
        if errors:
            for t in range(4):s.get('not there',t,status=9)
            s.get('name [units]',99,status=10)


def third_body():
    with Species('M') as s:
        lib.chem_species_third_body(s.handle)
        assert s.name == 'M' and s.parameterization.has_value
        p = s.parameterization
        assert lib.chem_parameterization_evaluate(C.byref(p),0,0,42.4) == 42.4


def copy(assignment=False):
    for parameterized in [False,True]:
        with Species('thing',BASIC) as s:
            s.set('foo',3,'bar');s.set('baz',1,42);s.set('qux',2,True)
            if parameterized:s.parameterization = Parameterization(constant=15.4,has_value=1)
            out = P()
            if assignment:
                assert lib.chem_species_create(b'old',0,None,C.byref(out)) == 0
                assert lib.chem_species_assign(C.byref(out),s.handle) == 0
            else:assert lib.chem_species_clone(s.handle,C.byref(out)) == 0
            obj = object.__new__(Species);obj.handle = out
            with obj:
                assert obj.name == 'thing'
                assert lib.chem_species_property_count(obj.handle,0) == 2
                assert obj.get('name [units]') == 1. and obj.get('name2 [units2]') == 2.
                assert obj.get('foo',3) == 'bar' and obj.get('baz',1) == 42 and obj.get('qux',2) == 1
                assert bool(obj.parameterization.has_value) == parameterized
                if parameterized:
                    p = obj.parameterization
                    assert lib.chem_parameterization_evaluate(C.byref(p),0,0,0) == 15.4
                s.set('foo',3,'changed')
                assert obj.get('foo',3) == 'bar'


def phase_species(mode=0):
    with Species('CO2') as s:
        member = PhaseSpecies(s.handle,1.23e-5 if mode==1 else 0,mode==1)
        assert lib.chem_species_name(member.species) == b'CO2'
        assert bool(member.has_diffusion) == (mode==1)
        if mode==2:member.diffusion,member.has_diffusion = 2.5e-6,1
        if mode:assert member.has_diffusion and member.diffusion == (1.23e-5 if mode==1 else 2.5e-6)


def phase(parameterized=False,names_only=False):
    names = ['foo','bar','baz'] if parameterized else ['species1','species2']
    species = [Species(n) for n in names]
    try:
        if parameterized:species[1].parameterization = Parameterization(constant=42,has_value=1)
        with Phase('gas' if parameterized else 'test_phase',[PhaseSpecies(s.handle,0,0) for s in species]) as p:
            assert p.name == ('gas' if parameterized else 'test_phase')
            assert p.count == len(names) and p.size == 2
            assert p.names() == (['gas.foo','gas.baz'] if parameterized else ['test_phase.species1','test_phase.species2'])
    finally:
        for s in species:lib.chem_species_free(s.handle)


def system(parameterized=False):
    species = [Species('foo'),Species('bar')]
    if parameterized:
        species.append(Species('paramSpecies'))
        species[-1].parameterization = Parameterization(constant=64.2,has_value=1)
    try:
        with Phase('gas',[PhaseSpecies(s.handle,0,0) for s in species]) as p:
            assert p.count == 2+parameterized and p.size == 2
            names = p.names(prefix=False)
            assert sorted(names) == ['bar','foo']
            if parameterized:assert p.names(prefix=False,permutation=[1,0]) == [names[1],names[0]]
    finally:
        for s in species:lib.chem_species_free(s.handle)


CASES = {
    'Species.StringConstructor':string_constructor,
    'Species.StringAndVectorConstructor':properties,
    'Species.GetProperty':partial(properties,True),
    'Species.SetThirdBody':third_body,
    'Species.CopyConstructor':copy,
    'Species.CopyAssignment':partial(copy,True),
    'PhaseSpecies.ConstructsWithoutDiffusionCoefficient':phase_species,
    'PhaseSpecies.ConstructsWithDiffusionCoefficient':partial(phase_species,1),
    'PhaseSpecies.SetDiffusionCoefficient':partial(phase_species,2),
    'Phase.Constructor':phase,
    'Phase.ConstructorWithParameterizedSpecies':partial(phase,True),
    'Phase.UniqueNamesWithParameterizedSpecies':partial(phase,True,True),
    'System.ConstructorWithAllParameters':system,
    'System.ConstructorWithParameterizedSpecies':partial(system,True),
}

bind('chem_version',C.c_char_p)
for component in ['major','minor','patch','tweak']:bind('chem_version_'+component,C.c_uint)


def version(component):
    value = getattr(lib,'chem_version'+('_'+component if component else ''))()
    # The source's version cases only call getters. Check consistency as well.
    parts = [lib.chem_version_major(),lib.chem_version_minor(),lib.chem_version_patch()]
    assert lib.chem_version().decode() == '.'.join(map(str,parts))
    assert value is not None


for name,component in [('FullVersion',''),('VersionMajor','major'),('VersionMinor','minor'),('VersionPatch','patch'),('VersionTweak','tweak')]:
    CASES['Version.'+name] = partial(version,component)
