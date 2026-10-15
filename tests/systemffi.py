"""Owned species/phase metadata bindings."""
from chemffi import *


class PropertyValue(C.Union):
    _fields_ = [('real',D),('integer',I),('boolean',I),('text',C.c_char_p)]


class Property(C.Structure):
    _fields_ = [('name',C.c_char_p),('type',I),('value',PropertyValue)]


class PhaseSpecies(C.Structure):
    _fields_ = [('species',P),('diffusion',D),('has_diffusion',I)]


bind('chem_parameterization_evaluate',D,C.POINTER(Parameterization),D,D,D)
bind('chem_species_create',I,C.c_char_p,Z,C.POINTER(Property),C.POINTER(P))
bind('chem_species_free',None,P)
bind('chem_species_clone',I,P,C.POINTER(P))
bind('chem_species_assign',I,C.POINTER(P),P)
bind('chem_species_name',C.c_char_p,P)
bind('chem_species_property_count',Z,P,I)
bind('chem_species_set_property',I,P,Property)
bind('chem_species_get_property',I,P,C.c_char_p,I,C.POINTER(PropertyValue))
bind('chem_species_parameterize',None,P,Parameterization)
bind('chem_species_parameterization',Parameterization,P)
bind('chem_species_third_body',None,P)
bind('chem_phase_create',I,C.c_char_p,Z,C.POINTER(PhaseSpecies),C.POINTER(P))
bind('chem_phase_free',None,P)
bind('chem_phase_name',C.c_char_p,P)
for suffix in ['count','state_size']:bind('chem_phase_'+suffix,Z,P)
bind('chem_phase_species',C.POINTER(PhaseSpecies),P,Z)
bind('chem_phase_variable_name',Z,P,Z,I,C.POINTER(Z),P,Z)
bind('chem_surface_rate',I,C.POINTER(PhaseSpecies),D,I,C.POINTER(Rate))
PROPERTY_FIELDS = ['real','integer','boolean','text']


def property_(name,type,value):
    union = PropertyValue()
    setattr(union,PROPERTY_FIELDS[type],value.encode() if type==3 else value)
    return Property(name.encode(),type,union)


class Species:
    def __init__(self,name,properties=()):
        self.handle = P()
        data = (Property*len(properties))(*[property_(*p) for p in properties])
        assert lib.chem_species_create(name.encode(),len(data),data,C.byref(self.handle)) == 0
    def __enter__(self):return self
    def __exit__(self,*args):lib.chem_species_free(self.handle)
    def get(self,name,type=0,status=0):
        out = PropertyValue()
        assert lib.chem_species_get_property(self.handle,name.encode(),type,C.byref(out)) == status
        if status:return
        value = getattr(out,PROPERTY_FIELDS[type])
        return value.decode() if type==3 else value
    def set(self,name,type,value):assert lib.chem_species_set_property(self.handle,property_(name,type,value)) == 0
    @property
    def parameterization(self):return lib.chem_species_parameterization(self.handle)
    @parameterization.setter
    def parameterization(self,p):lib.chem_species_parameterize(self.handle,p)
    @property
    def name(self):return lib.chem_species_name(self.handle).decode()


class Phase:
    def __init__(self,name,species):
        self.handle = P()
        data = (PhaseSpecies*len(species))(*species)
        assert lib.chem_phase_create(name.encode(),len(species),data,C.byref(self.handle)) == 0
    def __enter__(self):return self
    def __exit__(self,*args):lib.chem_phase_free(self.handle)
    @property
    def count(self):return lib.chem_phase_count(self.handle)
    @property
    def size(self):return lib.chem_phase_state_size(self.handle)
    @property
    def name(self):return lib.chem_phase_name(self.handle).decode()
    def member(self,i):return lib.chem_phase_species(self.handle,i).contents
    def names(self,prefix=True,permutation=None):
        perm = (Z*len(permutation))(*permutation) if permutation is not None else None
        result = []
        for i in range(self.size):
            length = lib.chem_phase_variable_name(self.handle,i,prefix,perm,None,0)
            assert length != Z(-1).value
            out = C.create_string_buffer(length+1)
            assert lib.chem_phase_variable_name(self.handle,i,prefix,perm,out,len(out)) == length
            result.append(out.value.decode())
        return result
