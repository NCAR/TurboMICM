"""C builder and external model fixture bindings (test tooling only)."""
from chemffi import *

MODEL = C.CFUNCTYPE(None,P,DP,DP,P)
UPDATE = C.CFUNCTYPE(None,P,P)
FINALIZE = C.CFUNCTYPE(None,P,P)
CONFIGURE = C.CFUNCTYPE(I,C.POINTER(Data),C.POINTER(ExternalModel),C.POINTER(Constraint),P)


class ModelDefinition(C.Structure):
    _fields_ = [('variable_count',Z),('parameter_count',Z),('constraint_count',Z),
                ('variable_names',C.POINTER(C.c_char_p)),('parameter_names',C.POINTER(C.c_char_p)),
                ('initial_variables',DP),('initial_parameters',DP),('absolute_tolerances',DP),
                ('used_species_count',Z),('used_species_names',C.POINTER(C.c_char_p)),('configure',CONFIGURE),('context',P)]


class BuilderSpec(C.Structure):
    _fields_ = [('base',Data),('model_count',Z),('models',C.POINTER(ModelDefinition)),('ignore_unused_species',I)]


bind('chem_builder_build',I,C.POINTER(BuilderSpec),C.POINTER(P))
for kind in ['species','parameter']:bind('chem_data_'+kind+'_index',Z,C.POINTER(Data),C.c_char_p)
bind('chem_maximum_cells',Z,P)
bind('chem_mechanism_phase_name',C.c_char_p,P)
bind('chem_cell_count',Z,P)


def lookup(data,name,parameter=False):
    fn = lib.chem_data_parameter_index if parameter else lib.chem_data_species_index
    index = fn(data,name.encode())
    assert index != Z(-1).value,name
    return index


class Definition:
    def __init__(self,variables=(),parameters=(),constraints=0,used=()):
        self.keep=[];self.errors=[];self.finalized=[]
        def array(type,values):
            out=(type*len(values))(*values);self.keep.append(out);return out
        self.array=array
        def configure(data,process,cs,context):
            try:
                process.contents.context=None
                self.configure(data,process.contents,cs)
                return 0
            except Exception as error:self.errors.append(error);return 1
        self.callback=CONFIGURE(configure)
        self.definition=ModelDefinition(variable_count=len(variables),parameter_count=len(parameters),constraint_count=constraints,
            variable_names=array(C.c_char_p,[n.encode() for n in variables]),
            parameter_names=array(C.c_char_p,[n.encode() for n in parameters]),
            used_species_count=len(used),used_species_names=array(C.c_char_p,[n.encode() for n in used]),configure=self.callback)
    def attach(self,process,forcing=None,jacobian=None,update=None,finalize=None):
        def guarded(fn):
            def run(*args):
                try:fn(*args)
                except Exception as error:self.errors.append(error)
            return run
        for field,signature,fn in [('forcing',MODEL,forcing),('negative_jacobian',MODEL,jacobian),('update',UPDATE,update)]:
            if fn is not None:
                cb=signature(guarded(fn));self.keep.append(cb);setattr(process,field,C.cast(cb,P))
        def finish(m,ctx):
            self.finalized.append((lib.chem_mechanism_data(m).contents.species_count,lib.chem_maximum_cells(m)))
            if finalize:finalize(m,ctx)
        cb=FINALIZE(guarded(finish));self.keep.append(cb);process.finalize=C.cast(cb,P)
    def check(self):
        if self.errors:raise self.errors.pop(0)


class Builder:
    def __init__(self,names,reactions=(),constraints=(),parameters=(),models=(),ignore_unused=True,maximum_cells=0):
        self.base=Mechanism(names,reactions,constraints,parameters,create=False)
        self.base.data.maximum_cells=maximum_cells
        self.base.data.phase_name=b'gas'
        self.models=list(models)
        self.definitions=(ModelDefinition*len(models))(*[m.definition for m in models])
        self.spec=BuilderSpec(self.base.data,len(models),self.definitions,int(ignore_unused))
    def build(self,cells=1,status=0):
        m=object.__new__(Mechanism);m.handle,m.state=P(),P();m.keep=[self]
        actual=lib.chem_builder_build(C.byref(self.spec),C.byref(m.handle))
        for model in self.models:model.check()
        assert actual==status,(actual,status)
        if actual:return None
        m.data=lib.chem_mechanism_data(m.handle).contents;m.n=m.data.species_count;m.nc=cells
        assert lib.chem_state_create(m.handle,cells,C.byref(m.state))==0
        m.y=lib.chem_concentrations(m.state);m.custom=lib.chem_parameters(m.state)
        return m
    def check(self):
        for model in self.models:model.check()
