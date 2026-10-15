"""Ports of species, phase, system, and version API behavior."""
import math
from api import *
from harness import case, equal


class Phase(C.Structure):
    _fields_ = [('name',C.c_char_p), ('count',Z), ('species',C.POINTER(Species))]


bind('chem_species_init',I,C.POINTER(Species),C.c_char_p)
bind('chem_species_destroy',None,C.POINTER(Species))
bind('chem_species_copy',I,C.POINTER(Species),C.POINTER(Species))
bind('chem_species_set_property',I,C.POINTER(Species),C.c_char_p,I,P)
bind('chem_species_get_property',I,C.POINTER(Species),C.c_char_p,I,P)
bind('chem_species_parameter',D,C.POINTER(Species),C.POINTER(Conditions))
bind('chem_phase_init',I,C.POINTER(Phase),C.c_char_p,C.POINTER(Species),Z)
bind('chem_phase_destroy',None,C.POINTER(Phase))
bind('chem_phase_state_size',Z,C.POINTER(Phase))
bind('chem_phase_names',I,C.POINTER(Phase),I,ZP,C.POINTER(C.POINTER(C.c_char_p)),ZP)
bind('chem_names_free',None,C.POINTER(C.c_char_p),Z)
bind('chem_version',C.c_char_p)
bind('chem_version_component',C.c_uint,C.c_uint)


def property_set(s,key,type_,value):
    obj = C.create_string_buffer(value.encode()) if type_ == 3 else (D(value) if type_ == 0 else I(value))
    equal(lib.chem_species_set_property(C.byref(s),key.encode(),type_,C.cast(C.byref(obj),P)),0)


def property_get(s,key,type_):
    obj = C.c_char_p() if type_ == 3 else (D() if type_ == 0 else I())
    equal(lib.chem_species_get_property(C.byref(s),key.encode(),type_,C.cast(C.byref(obj),P)),0)
    return obj.value


@case(*['species:Species.'+name for name in ('StringConstructor','StringAndVectorConstructor','GetProperty','SetThirdBody','CopyConstructor','CopyAssignment')])
def species(identifier):
    name = identifier.split('.')[-1]
    s,copy = Species(),Species()
    equal(lib.chem_species_init(C.byref(s),b'M' if name == 'SetThirdBody' else b'thing'),0)
    callback_type = C.CFUNCTYPE(D,C.POINTER(Conditions),DP,P)
    callback = callback_type(lambda c,u,ctx: 100*c.contents.temperature)
    constant = callback_type(lambda c,u,ctx: 15.4)
    try:
        if name == 'StringConstructor':
            equal(s.name,b'thing'); equal(s.parameterized,0)
            s.parameterized = 2; s.parameterization = C.cast(callback,P)
            equal(lib.chem_species_parameter(C.byref(s),C.byref(Conditions(12.5,0,0))),1250)
        elif name == 'SetThirdBody':
            s.parameterized = 1; equal(s.name,b'M'); equal(s.parameterized,1)
            equal(lib.chem_species_parameter(C.byref(s),C.byref(Conditions(0,0,42.4))),42.4)
        else:
            property_set(s,'name [units]',0,1.0); property_set(s,'name2 [units2]',0,2.0)
            equal(s.nproperties,2)
            equal(property_get(s,'name [units]',0),1.0); equal(property_get(s,'name2 [units2]',0),2.0)
            if name == 'GetProperty':
                value = D()
                for type_ in range(4): equal(lib.chem_species_get_property(C.byref(s),b'not there',type_,C.byref(value)),3)
                equal(lib.chem_species_get_property(C.byref(s),b'name [units]',4,C.byref(value)),1)
            if name.startswith('Copy'):
                property_set(s,'foo',3,'bar'); property_set(s,'baz',1,42); property_set(s,'qux',2,1)
                for parameterized in (False,True):
                    s.parameterized = 2 if parameterized else 0; s.parameterization = C.cast(constant,P)
                    equal(lib.chem_species_copy(C.byref(copy),C.byref(s)),0)
                    equal(copy.name,b'thing'); equal(copy.nproperties,5)
                    for key,type_,expected in [('name [units]',0,1.0),('name2 [units2]',0,2.0),('foo',3,b'bar'),('baz',1,42),('qux',2,1)]:
                        equal(property_get(copy,key,type_),expected)
                    equal(bool(copy.parameterized),parameterized)
                    if parameterized: equal(lib.chem_species_parameter(C.byref(copy),C.byref(Conditions())),15.4)
    finally:
        lib.chem_species_destroy(C.byref(copy)); lib.chem_species_destroy(C.byref(s))


@case(*['phase:PhaseSpecies.'+name for name in ('ConstructsWithoutDiffusionCoefficient','ConstructsWithDiffusionCoefficient','SetDiffusionCoefficient')])
def phase_species(identifier):
    s = Species(); equal(lib.chem_species_init(C.byref(s),b'CO2'),0)
    try:
        equal(s.name,b'CO2'); assert math.isnan(s.diffusion)
        if 'Without' not in identifier:
            s.diffusion = 2.5e-6 if 'SetDiffusion' in identifier else 1.23e-5
            equal(s.diffusion,2.5e-6 if 'SetDiffusion' in identifier else 1.23e-5)
    finally: lib.chem_species_destroy(C.byref(s))


@case('phase:Phase.Constructor','phase:Phase.ConstructorWithParameterizedSpecies','phase:Phase.UniqueNamesWithParameterizedSpecies',
      'system:System.ConstructorWithAllParameters','system:System.ConstructorWithParameterizedSpecies')
def phase_system(identifier):
    parameterized = 'Parameterized' in identifier
    plain = identifier == 'phase:Phase.Constructor'
    names = ['species1','species2'] if plain else ['foo','bar']
    if parameterized:
        names = ['foo','bar','paramSpecies'] if identifier.startswith('system') else ['foo','bar','baz']
    values = (Species*len(names))(); phase = Phase()
    try:
        for s,n in zip(values,names): equal(lib.chem_species_init(C.byref(s),n.encode()),0)
        if parameterized: values[2 if identifier.startswith('system') else 1].parameterized = 2
        label = b'test_phase' if plain else b'gas'
        equal(lib.chem_phase_init(C.byref(phase),label,values,len(values)),0)
        equal(phase.name,label); equal(phase.count,len(values)); equal(lib.chem_phase_state_size(C.byref(phase)),2)
        output = C.POINTER(C.c_char_p)(); count = Z()
        qualified = int(identifier.startswith('phase'))
        equal(lib.chem_phase_names(C.byref(phase),qualified,None,C.byref(output),C.byref(count)),0)
        expected = [n for i,n in enumerate(names) if not values[i].parameterized]
        prefix = label.decode()+'.' if qualified else ''
        equal([output[i].decode() for i in range(count.value)],[prefix+n for n in expected])
        lib.chem_names_free(output,count)
        if identifier.startswith('system') and parameterized:
            equal(lib.chem_phase_names(C.byref(phase),0,inds([1,0]),C.byref(output),C.byref(count)),0)
            equal([output[i].decode() for i in range(count.value)],expected[::-1]); lib.chem_names_free(output,count)
    finally:
        lib.chem_phase_destroy(C.byref(phase))
        for s in values: lib.chem_species_destroy(C.byref(s))


@case(*['version:Version.'+name for name in ('FullVersion','VersionMajor','VersionMinor','VersionPatch','VersionTweak')])
def version(identifier):
    name = identifier.split('.')[-1]
    if name == 'FullVersion': assert lib.chem_version()
    else: lib.chem_version_component(['VersionMajor','VersionMinor','VersionPatch','VersionTweak'].index(name))
