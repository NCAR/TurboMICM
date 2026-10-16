"""Small ctypes binding used exclusively by the build and tests."""
import ctypes as C
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
lib = C.CDLL(str(ROOT/'lib/libchem.so'))
P = C.c_void_p
I, Z, D = C.c_int, C.c_size_t, C.c_double
IP, DP, BP = C.POINTER(I), C.POINTER(D), C.POINTER(C.c_ubyte)


class Conditions(C.Structure):
    _fields_ = [('temperature',D),('pressure',D),('air_density',D)]


class Rate(C.Structure):
    _fields_ = [('type',I),('p',D*32),('custom_index',I),('coefficient_count',I),('function',P),('context',P)]


class Stats(C.Structure):
    _fields_ = [(name,Z) for name in ('internal_steps','accepted','rejected','function_calls',
                'jacobian_updates','decompositions','solves','constraint_init_iterations')] + [('final_time',D)]


class Rosenbrock(C.Structure):
    _fields_ = [('stages',I),('max_steps',I),('new_function',I*6)] + [(name,D) for name in (
        'round_off','factor_min','factor_max','rejection_factor','safety_factor','h_min','h_max','h_start','estimator_order')] + [
        ('a',D*15),('c',D*15),('m',D*6),('e',D*6),('alpha',D*6),('gamma',D*6),
        ('init_max_iterations',I),('init_max_backtracks',I),('init_tolerance',D),('init_backtrack_factor',D),
        ('init_sufficient_decrease',D),('controller',P),('controller_context',P)]


class BackwardEuler(C.Structure):
    _fields_ = [('small',D),('h_start',D),('max_iterations',I),('reduction_count',I),('reductions',D*16)]


class External(C.Structure):
    _fields_ = [('elements',Z),('algebraic_count',Z),('rows',IP),('columns',IP),('algebraic_species',IP)] + [
        (name,P) for name in ('bind','update','initialize_constraints','force','jacobian',
                              'constraint_force','constraint_jacobian','context')]


class ReactionInfo(C.Structure):
    _fields_ = [('reactants',Z),('products',Z),('parameterized_reactants',Z),
                ('reactant_ids',IP),('product_ids',IP),('parameterized_ids',IP),('yields',DP),('rate',C.POINTER(Rate))]


class JacobianComparison(C.Structure):
    _fields_ = [('passed',I),('max_absolute_error',D),('max_relative_error',D),
                ('worst_cell',Z),('worst_row',Z),('worst_column',Z),('worst_analytical',D),('worst_finite_difference',D)]


ForcingFunction=C.CFUNCTYPE(None,Z,Z,DP,DP,P)


def bind(name, result, *args):
    function = getattr(lib,'ch_'+name)
    function.restype, function.argtypes = result,list(args)
    globals()[name] = function
    return function


bind('status_string',C.c_char_p,I)
bind('rate_default',Rate,I)
bind('rate_prepare',None,C.POINTER(Rate))
bind('rate_calculate',D,C.POINTER(Rate),C.POINTER(Conditions),DP)
bind('model_create',P)
bind('model_free',None,P)
bind('model_add_species',I,P,C.c_char_p,D,I)
bind('model_set_species_property',I,P,I,C.c_char_p,D)
bind('model_add_parameter',I,P,C.c_char_p,D)
bind('model_add_reaction',I,P,Z,IP,Z,IP,DP,C.POINTER(Rate))
bind('model_add_linear',I,P,C.c_char_p,I,Z,IP,DP,D)
bind('model_add_equilibrium',I,P,C.c_char_p,I,Z,IP,DP,Z,IP,DP,D,D)
bind('model_add_external',I,P,C.POINTER(External))
bind('model_finalize',I,P,I,I)
for name in ('species_count','reaction_count','parameter_count','jacobian_size','constraint_count'):
    bind('model_'+name,Z,P)
for name in ('species_name','parameter_name'):
    bind('model_'+name,C.c_char_p,P,Z)
for name in ('species_index','parameter_index'):
    bind('model_'+name,I,P,C.c_char_p)
bind('model_is_algebraic',I,P,Z)
bind('model_rate',C.POINTER(Rate),P,Z)
bind('model_jacobian_index',I,P,Z,Z)
bind('state_create',P,P,Z)
bind('state_clone',P,P)
bind('state_free',None,P)
bind('state_model',P,P)
bind('state_cells',Z,P)
for name in ('values','parameters','rates'):
    bind('state_'+name,DP,P)
bind('state_conditions',C.POINTER(Conditions),P)
for name in ('set','set_parameter'):
    bind('state_'+name,I,P,C.c_char_p,DP,Z)
bind('state_tolerances',I,P,D,DP,Z)
bind('update_rates',I,P)
for name in ('force','jacobian','constraint_force','constraint_jacobian'):
    bind(name,None,P,DP,DP)
bind('rosenbrock_default',Rosenbrock,I)
bind('backward_euler_default',BackwardEuler)
bind('solve',I,P,D,C.POINTER(Rosenbrock),C.POINTER(Stats))
bind('solve_many',I,C.POINTER(P),Z,D,C.POINTER(Rosenbrock),I,IP,C.POINTER(Stats))
bind('initialize_constraints',I,P,C.POINTER(Rosenbrock),C.POINTER(Stats))
bind('solve_backward_euler',I,P,D,C.POINTER(BackwardEuler),C.POINTER(Stats))
bind('backward_euler_converged',I,P,DP,DP,C.POINTER(BackwardEuler))
bind('normalized_error',D,P,DP,DP,DP)
bind('markowitz',I,Z,BP,IP)
bind('sparse_create',P,Z,BP)
bind('sparse_free',None,P)
bind('sparse_size',Z,P)
bind('sparse_index',I,P,Z,Z)
bind('lu_factor',I,P,Z,DP,I)
bind('lu_factor_copy',I,P,Z,DP,DP,I)
bind('lu_solve',None,P,Z,DP,DP)
bind('sparse_multiply',None,P,Z,DP,DP,DP)
bind('alpha_minus_jacobian',None,P,D,DP)
bind('state_set_lu',I,P,I)
bind('model_clone_builder',P,P)
bind('species_create',P,C.c_char_p)
bind('species_clone',P,P)
bind('species_free',None,P)
bind('species_name',C.c_char_p,P)
bind('species_property_count',Z,P,I)
bind('species_set_property',I,P,C.c_char_p,I,P)
bind('species_get_property',I,P,C.c_char_p,I,P)
bind('species_parameterize',I,P,I,D,D,D,D)
bind('species_is_parameterized',I,P)
bind('species_parameter_value',D,P,C.POINTER(Conditions))
bind('phase_create',P,C.c_char_p)
bind('phase_free',None,P)
bind('phase_name',C.c_char_p,P)
bind('phase_add',I,P,P)
bind('phase_species_count',Z,P)
bind('phase_state_size',Z,P)
bind('phase_species',P,P,Z)
bind('phase_state_name',C.c_char_p,P,Z,I)
bind('model_add_species_spec',I,P,P)
bind('model_add_phase',I,P,P,I)
bind('model_species_info',P,P,Z)
bind('model_builder_species_count',Z,P)
bind('state_copy',I,P,P)
bind('state_move',None,C.POINTER(P),C.POINTER(P))
bind('state_relative_tolerance',D,P)
bind('state_absolute_tolerances',DP,P)
bind('state_set_exact',I,P,C.c_char_p,DP,Z)
bind('state_set_parameter_exact',I,P,C.c_char_p,DP,Z)
bind('state_set_parameter_matrix',I,P,DP,Z,Z)
bind('version',C.c_char_p)
bind('version_component',I,I)
bind('model_reaction_info',I,P,Z,C.POINTER(ReactionInfo))
bind('model_rate_group',I,P,I,C.POINTER(Z),C.POINTER(Z))
bind('model_rate_in_group',C.POINTER(Rate),P,I,Z)
bind('model_species_used',I,P,Z)
bind('model_kinetic_pattern',I,P,BP,Z)
bind('add_force',I,P,DP,DP)
bind('add_negative_jacobian',I,P,DP,DP)
bind('finite_difference',I,Z,Z,DP,D,ForcingFunction,P,DP)
bind('compare_jacobian',JacobianComparison,Z,Z,BP,DP,DP,D,D)
bind('check_sparsity',JacobianComparison,Z,Z,BP,DP,D)


def doubles(values):
    return (D*len(values))(*values)


def ints(values):
    return (I*len(values))(*values)


def ptr(array):
    return array.ctypes.data_as(DP)


def require(status):
    assert status >= 0, status_string(status).decode()
    return status


def near(actual, expected, rtol=1e-13, atol=0):
    actual,expected = np.asarray(actual),np.asarray(expected)
    assert actual.shape == expected.shape, (actual.shape,expected.shape)
    delta = np.abs(actual-expected)
    assert np.all(np.isfinite(actual)) and np.all(delta <= atol+rtol*np.abs(expected)), (actual,expected,delta)


class Model:
    def __init__(self, names=()):
        self.handle = model_create()
        assert self.handle
        self.ids = {}
        for name in names:
            self.species(name)

    def close(self):
        if self.handle:
            model_free(self.handle)
            self.handle = None

    def __enter__(self):
        return self

    def __exit__(self,*_):
        self.close()

    def species(self, name, atol=1e-3, third_body=0, **props):
        self.ids[name] = require(model_add_species(self.handle,name.encode(),atol,third_body))
        for key,value in props.items():
            require(model_set_species_property(self.handle,self.ids[name],key.encode(),value))
        return self.ids[name]

    def parameter(self, name, initial=0):
        return require(model_add_parameter(self.handle,name.encode(),initial))

    def reaction(self, reactants, products, rate):
        ids = ints([self.ids[r] for r in reactants])
        products = list(products.items()) if isinstance(products,dict) else list(products)
        prod = ints([self.ids[p] for p,_ in products])
        return require(model_add_reaction(self.handle,len(ids),ids,len(prod),prod,
                       doubles([y for _,y in products]),C.byref(rate)))

    def linear(self, row, terms, constant, name='linear'):
        terms = list(terms.items())
        return require(model_add_linear(self.handle,name.encode(),self.ids[row],len(terms),
            ints([self.ids[s] for s,_ in terms]),doubles([v for _,v in terms]),constant))

    def equilibrium(self, row, reactants, products, constant, enthalpy=0, name='equilibrium'):
        r,p = list(reactants.items()),list(products.items())
        return require(model_add_equilibrium(self.handle,name.encode(),self.ids[row],len(r),
            ints([self.ids[s] for s,_ in r]),doubles([v for _,v in r]),len(p),
            ints([self.ids[s] for s,_ in p]),doubles([v for _,v in p]),constant,enthalpy))

    def finalize(self, reorder=False, ignore_unused=True):
        require(model_finalize(self.handle,reorder,ignore_unused))
        self.names = [model_species_name(self.handle,i).decode() for i in range(model_species_count(self.handle))]
        self.indices = {name:i for i,name in enumerate(self.names)}
        return self


class State:
    def __init__(self,model,cells=1):
        self.model,self.cells = model,cells
        self.handle = state_create(model.handle,cells)
        assert self.handle

    def close(self):
        if self.handle:
            state_free(self.handle)
            self.handle = None

    def __enter__(self):
        return self

    def __exit__(self,*_):
        self.close()

    @property
    def y(self):
        return np.ctypeslib.as_array(state_values(self.handle),shape=(len(self.model.names),self.cells))

    @property
    def rates(self):
        return np.ctypeslib.as_array(state_rates(self.handle),shape=(model_reaction_count(self.model.handle),self.cells))

    @property
    def custom(self):
        return np.ctypeslib.as_array(state_parameters(self.handle),shape=(model_parameter_count(self.model.handle),self.cells))

    def set(self,name,value):
        value = np.atleast_1d(value).astype(float)
        require(state_set(self.handle,name.encode(),ptr(value),value.size))

    def parameter(self,name,value):
        value = np.atleast_1d(value).astype(float)
        require(state_set_parameter(self.handle,name.encode(),ptr(value),value.size))

    def conditions(self,temperature=301.24,pressure=101325,air=None):
        cond = state_conditions(self.handle)
        for c in range(self.cells):
            t = np.broadcast_to(temperature,(self.cells,))[c]
            p = np.broadcast_to(pressure,(self.cells,))[c]
            a = p/(1.380649e-23*6.02214076e23*t) if air is None else np.broadcast_to(air,(self.cells,))[c]
            cond[c] = Conditions(t,p,a)

    def update(self):
        require(update_rates(self.handle))

    def tolerance(self,rtol=1e-6,atol=1e-3):
        atol = np.atleast_1d(atol).astype(float)
        require(state_tolerances(self.handle,rtol,ptr(atol),atol.size))

    def solve(self,dt,method=3,params=None,allow_failure=False):
        stats = Stats()
        if method == 'be':
            params = params or backward_euler_default()
            status = solve_backward_euler(self.handle,dt,C.byref(params),C.byref(stats))
        else:
            params = params or rosenbrock_default(method)
            status = solve(self.handle,dt,C.byref(params),C.byref(stats))
        if not allow_failure:
            require(status)
        return status,stats

    def forcing(self, constraints=False):
        result = np.zeros_like(self.y)
        force(self.handle,ptr(self.y),ptr(result))
        if constraints:
            constraint_force(self.handle,ptr(self.y),ptr(result))
        return result

    def jacobian(self,constraints=False):
        size = model_jacobian_size(self.model.handle)
        values = np.zeros((size,self.cells))
        jacobian(self.handle,ptr(self.y),ptr(values))
        if constraints:
            constraint_jacobian(self.handle,ptr(self.y),ptr(values))
        n = len(self.model.names)
        result = np.zeros((self.cells,n,n))
        for r in range(n):
            for c in range(n):
                i = model_jacobian_index(self.model.handle,r,c)
                if i >= 0:
                    result[:,r,c] = -values[i,:]
        return result
