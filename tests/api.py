"""ctypes bindings used only by the test tooling, never by the benchmark."""
import ctypes as C
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
lib = C.CDLL(str(ROOT / 'bin' / ('libchem.dylib' if sys.platform == 'darwin' else 'libchem.so')))
Z = C.c_size_t
D = C.c_double
I = C.c_int
P = C.c_void_p
DP = C.POINTER(D)
ZP = C.POINTER(Z)
MAX = Z(-1).value


class Conditions(C.Structure):
    _fields_ = [(name, D) for name in ('temperature','pressure','air_density')]


class Rate(C.Structure):
    _fields_ = [('type', I)] + [(name, D) for name in (
        'A B C D E k0_A k0_B k0_C k0_D kinf_A kinf_B kinf_C kinf_D kint_A kint_B kint_C kint_D Fc N X Y a0 n k0 z'.split())] + [
        ('branch', I), ('diffusion', D), ('molecular_weight', D), ('probability', D), ('coefficients', D*10),
        ('ncoefficients', Z), ('label', C.c_char_p), ('custom_index', Z), ('function', P), ('context', P),('mean_free_speed_factor',D)]


class Species(C.Structure):
    _fields_ = [('name', C.c_char_p)] + [(name,D) for name in ('initial','absolute_tolerance','molecular_weight','diffusion')] + [
        ('parameterized',I), ('parameterization',P), ('context',P), ('variable',Z), ('properties',P), ('nproperties',Z)]


class Reaction(C.Structure):
    _fields_ = [('nreactants',Z), ('nproducts',Z), ('reactants',ZP), ('products',ZP), ('yields',DP), ('rate',Rate)]


class Constraint(C.Structure):
    _fields_ = [('type',I), ('algebraic',Z), ('nterms',Z), ('nreactants',Z), ('species',ZP), ('coefficients',DP),
                ('constant',D), ('residual',P), ('jacobian',P), ('context',P),('name',C.c_char_p),
                ('delta_enthalpy',D),('reference_temperature',D),('custom_index',Z)]


class Range(C.Structure):
    _fields_ = [('first',Z), ('count',Z)]


class Triple(C.Structure):
    _fields_ = [('target',Z), ('left',Z), ('right',Z)]


class JacJob(C.Structure):
    _fields_ = [(name,Z) for name in ('reaction','column','omitted','nreactants','nproducts')]


class ModelData(C.Structure):
    _fields_ = [(name,Z) for name in ('nspecies','nvariables','nreactions','ncustom','nconstraints')] + [
        ('species',C.POINTER(Species)), ('reactions',C.POINTER(Reaction)), ('constraints',C.POINTER(Constraint)),
        ('custom_names',C.POINTER(C.c_char_p)), ('custom_defaults',DP), ('default_temperature',D), ('default_pressure',D),
        ('nnz',Z), ('variable_species',ZP), ('diagonal',ZP), ('index',C.POINTER(I)), ('algebraic',C.POINTER(C.c_ubyte))] + [
        (name,C.POINTER(Range)) for name in ('lower','upper','scale','updates')] + [
        (name,ZP) for name in ('lower_ids','upper_ids','scale_ids','lower_columns','upper_columns')] + [
        ('update_ops',C.POINTER(Triple)), ('njobs',Z), ('jobs',C.POINTER(JacJob))] + [
        (name,P) for name in ('force_kernel','jacobian_kernel','external_force','external_jacobian','external_context')] + [('finalized',I),('declared_sparsity',C.POINTER(C.c_ubyte)),('ignore_unused',I),('setup',P),('setup_context',P),('external_models',P),('nexternal_models',Z)]


class StateData(C.Structure):
    _fields_ = [('model',C.POINTER(ModelData)), ('cells',Z), ('stride',Z)] + [
        (name,DP) for name in ('variables','custom','rates')] + [('conditions',C.POINTER(Conditions)), ('absolute_tolerance',DP),
        ('relative_tolerance',D)] + [(name,DP) for name in ('jacobian','work','new_y','forcing','error')] + [('stages',DP*6)]


class Parameters(C.Structure):
    _fields_ = [('stages',I), ('new_function',I*6), ('a',D*15), ('c',D*15)] + [
        (name,D*6) for name in ('m','e','alpha','gamma')] + [('order',D), ('max_steps',Z)] + [
        (name,D) for name in ('roundoff','factor_min','factor_max','rejection_factor','safety','h_min','h_max','h_start')] + [
        ('init_iterations',Z), ('init_backtracks',Z)] + [(name,D) for name in ('init_tolerance','init_backtrack_factor','init_sufficient_decrease')]


class Solver(C.Structure):
    _fields_ = [('parameters',Parameters), ('controller',P), ('controller_context',P), ('factor',P), ('solve',P)]


class Stats(C.Structure):
    _fields_ = [(name,C.c_uint64) for name in ('internal_steps','accepted','rejected','function_calls','jacobian_updates','decompositions','solves','constraint_init_iterations')] + [('final_time',D)]


class Euler(C.Structure):
    _fields_ = [('small',D), ('h_start',D), ('max_steps',Z), ('nreductions',Z), ('reductions',D*5), ('factor',P), ('solve',P)]


def bind(name, result, *arguments):
    fn = getattr(lib,name)
    fn.restype, fn.argtypes = result, arguments
    return fn


MP, SP = C.POINTER(ModelData), C.POINTER(StateData)
class External(C.Structure):
    _fields_ = [(name,P) for name in ('forcing','jacobian','update','setup','context')]

ExternalKernel = C.CFUNCTYPE(None,MP,SP,DP,DP,P)
ExternalUpdate = C.CFUNCTYPE(None,SP,P)
ExternalSetup = C.CFUNCTYPE(None,MP,P)
bind('chem_model_external',I,MP,External)
bind('chem_model_custom',I,MP,C.c_char_p,D,ZP)
bind('chem_rate',Rate,I)
bind('chem_rate_value',D,C.POINTER(Rate),C.POINTER(Conditions),DP)
bind('chem_model_create',MP)
bind('chem_model_free',None,MP)
bind('chem_model_species',I,MP,C.c_char_p,D,I,ZP)
bind('chem_model_reaction',I,MP,ZP,Z,ZP,DP,Z,Rate)
bind('chem_model_constraint',I,MP,C.POINTER(Constraint))
bind('chem_constraint_sparsity',I,MP,C.POINTER(C.c_ubyte),Z)
bind('chem_model_reorder',I,MP)
bind('chem_model_finalize',I,MP)
bind('chem_model_sparsity',I,MP,C.POINTER(C.c_ubyte),Z)
bind('chem_model_lambda',C.POINTER(Rate),MP,C.c_char_p)
bind('chem_rate_offset',Z,MP,I)
bind('chem_status_string',C.c_char_p,I)
for _name in ('chem_species_id','chem_variable_id','chem_custom_id'):
    bind(_name,Z,MP,C.c_char_p)
bind('chem_state_create',SP,MP,Z)
bind('chem_state_copy',SP,SP)
bind('chem_state_free',None,SP)
for _name in ('chem_state_set','chem_state_set_custom'):
    bind(_name,I,SP,C.c_char_p,DP,Z)
bind('chem_state_tolerances',I,SP,D,DP,Z)
bind('chem_state_custom_array',I,SP,DP,Z,Z)
bind('chem_calculate_rates',None,SP)
for _name in ('chem_forcing','chem_jacobian'):
    bind(_name,None,SP,DP,DP)
bind('chem_constraints',None,SP,DP,DP,DP)
for _name in ('chem_generic_forcing','chem_generic_jacobian'):
    bind(_name,None,MP,Z,Z,DP,DP,DP)
for _name in ('chem_lu_mozart','chem_lu_doolittle'):
    bind(_name,None,MP,Z,Z,DP)
bind('chem_linear_solve',None,MP,Z,Z,DP,DP)
bind('chem_sparse_multiply',None,MP,Z,Z,DP,DP,DP)
bind('chem_lu_separate',I,MP,Z,Z,DP,DP,DP,I)
bind('chem_linear_solve_separate',None,MP,Z,Z,DP,DP,DP)
bind('chem_solver',Solver,I)
bind('chem_solve',I,C.POINTER(Solver),D,SP,C.POINTER(Stats))
bind('chem_initialize_constraints',I,C.POINTER(Solver),SP,C.POINTER(Stats))
bind('chem_normalized_error',D,SP,DP,DP,DP)
bind('chem_alpha_minus_jacobian',None,SP,D)
bind('chem_markowitz',I,Z,C.POINTER(C.c_ubyte),ZP)
bind('chem_euler_parameters',Euler)
bind('chem_euler_converged',I,SP,D,DP,DP)
bind('chem_backward_euler',I,C.POINTER(Euler),D,SP,C.POINTER(Stats))
bind('chem_advance',I,P,P,D,Z,SP,C.POINTER(Stats))
bind('chem_advance_many',I,P,P,D,Z,C.POINTER(SP),Z,I,C.POINTER(Stats))

ARRHENIUS,TROE,TERNARY,JPL19,BRANCHED,TUNNELING,TAYLOR,REVERSIBLE,USER,SURFACE,LAMBDA = range(11)


def arr(values):
    return (D*len(values))(*values)


def inds(values):
    return (Z*len(values))(*values)


def rate(kind=ARRHENIUS, **values):
    r = lib.chem_rate(kind)
    for key,value in values.items():
        if key == 'coefficients':
            for i,v in enumerate(value):
                r.coefficients[i] = v
            r.ncoefficients = len(value)
        else:
            setattr(r,key,value.encode() if isinstance(value,str) else value)
    return r


def value(r, temperature=301.24, pressure=0.0, air_density=42.2, custom=()):
    return lib.chem_rate_value(C.byref(r),C.byref(Conditions(temperature,pressure,air_density)),arr(custom))


class Model:
    def __init__(self, names=(), pointer=None):
        self.p = pointer if pointer is not None else lib.chem_model_create()
        assert self.p
        self.callbacks = []
        for name in names:
            self.species(name)

    @property
    def data(self):
        return self.p.contents

    def close(self):
        if self.p:
            lib.chem_model_free(self.p)
            self.p = None

    def species(self,name,initial=0.0,parameterized=0,**values):
        index = Z()
        assert lib.chem_model_species(self.p,name.encode(),initial,parameterized,C.byref(index)) == 0
        for key,val in values.items():
            setattr(self.data.species[index.value],key,val)
        return index.value

    def reaction(self,reactants,products,r=None):
        if r is None:
            r = rate()
        rs = [lib.chem_species_id(self.p,s.encode()) if isinstance(s,str) else s for s in reactants]
        ps = [lib.chem_species_id(self.p,s.encode()) if isinstance(s,str) else s for s,_ in products]
        code = lib.chem_model_reaction(self.p,inds(rs),len(rs),inds(ps),arr([v for _,v in products]),len(ps),r)
        assert code == 0, code

    def constraint(self,algebraic,terms,constant=1.0,equilibrium_reactants=None,name=None,delta_enthalpy=0.0,reference_temperature=298.15):
        algebraic = self.var(algebraic) if isinstance(algebraic,str) else algebraic
        species = inds([self.var(s) if isinstance(s,str) else s for s,_ in terms])
        coeffs = arr([v for _,v in terms])
        q = Constraint(int(equilibrium_reactants is not None),algebraic,len(terms),equilibrium_reactants or 0,species,coeffs,constant,None,None,None)
        q.name = name.encode() if name else None
        q.delta_enthalpy = delta_enthalpy; q.reference_temperature = reference_temperature
        assert lib.chem_model_constraint(self.p,C.byref(q)) == 0

    def var(self,name):
        return lib.chem_variable_id(self.p,name.encode())

    def finalize(self,reorder=False):
        if reorder:
            assert lib.chem_model_reorder(self.p) == 0
        assert lib.chem_model_finalize(self.p) == 0
        return self


class State:
    def __init__(self,model,cells=1,pointer=None):
        self.model = model
        self.p = pointer if pointer is not None else lib.chem_state_create(model.p,cells)
        assert self.p

    @property
    def data(self):
        return self.p.contents

    def close(self):
        if self.p:
            lib.chem_state_free(self.p)
            self.p = None

    def set(self,name,values,custom=False):
        if isinstance(values,(int,float)):
            values = [values]*self.data.cells
        code = (lib.chem_state_set_custom if custom else lib.chem_state_set)(self.p,name.encode(),arr(values),len(values))
        assert code == 0, code

    def get(self,name,cell=0):
        return self.data.variables[self.model.var(name)*self.data.stride+cell]

    def rates(self):
        lib.chem_calculate_rates(self.p)

    def tolerances(self,relative=1e-6,absolute=(1e-3,)):
        if isinstance(absolute,(int,float)):
            absolute = [absolute]
        assert lib.chem_state_tolerances(self.p,relative,arr(absolute),len(absolute)) == 0

    def solve(self,dt,method=3,solver=None):
        solver = solver if solver is not None else lib.chem_solver(method)
        stats = Stats()
        code = lib.chem_solve(C.byref(solver),dt,self.p,C.byref(stats))
        assert code == 0, (code,stats.final_time,stats.internal_steps)
        return stats

    def advance(self,dt,method=3,parameters=None,substeps=1):
        parameters = parameters if parameters is not None else (lib.chem_euler_parameters() if method == 0 else lib.chem_solver(method))
        step = lib.chem_step_euler if isinstance(parameters,Euler) else lib.chem_step_rosenbrock
        stats = Stats()
        code = lib.chem_advance(C.cast(step,P),C.byref(parameters),dt,substeps,self.p,C.byref(stats))
        assert code == 0,(code,stats.final_time,stats.internal_steps)
        return stats
