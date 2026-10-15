"""ctypes test adapter. Only build/test processes import Python."""
import ctypes as C
from pathlib import Path

D = C.c_double
I = C.c_int
Z = C.c_size_t
P = C.c_void_p
DP = C.POINTER(D)
IP = C.POINTER(I)


class Rate(C.Structure):
    _fields_ = [('type', I), ('p', D * 24), ('parameter', I), ('coefficient_count', I),
                ('third_body_count', I), ('function', P), ('context', P), ('name', C.c_char_p)]


class Product(C.Structure):
    _fields_ = [('species', I), ('yield_', D)]


class Reaction(C.Structure):
    _fields_ = [('reactant_count', Z), ('product_count', Z), ('reactants', IP),
                ('products', C.POINTER(Product)), ('rate', Rate)]


class Constraint(C.Structure):
    _fields_ = [('type', I), ('algebraic_species', I), ('reactant_count', Z), ('product_count', Z),
                ('reactants', C.POINTER(Product)), ('products', C.POINTER(Product)),
                ('constant', D), ('equilibrium_constant', D), ('enthalpy', D), ('gas_constant', D),
                ('reference_temperature', D), ('parameter', I), ('function', P), ('context', P), ('name', C.c_char_p)]


class MatrixIndex(C.Structure):
    _fields_ = [('row', Z), ('column', Z)]


class ExternalModel(C.Structure):
    _fields_ = [('nonzero_count', Z), ('nonzeros', C.POINTER(MatrixIndex)),
                ('forcing', P), ('negative_jacobian', P), ('update', P), ('context', P)]


class Data(C.Structure):
    _fields_ = [('species_count', Z), ('reaction_count', Z), ('parameter_count', Z),
                ('species_names', C.POINTER(C.c_char_p)), ('parameter_names', C.POINTER(C.c_char_p)),
                ('initial_concentrations', DP), ('initial_parameters', DP), ('absolute_tolerances', DP),
                ('temperature', D), ('pressure', D), ('reactions', C.POINTER(Reaction)),
                ('constraint_count', Z), ('constraints', C.POINTER(Constraint)),
                ('external_model_count', Z), ('external_models', C.POINTER(ExternalModel))]


class Stats(C.Structure):
    _fields_ = [(n, Z) for n in ('internal_steps', 'accepted', 'rejected', 'function_calls',
                'jacobian_updates', 'decompositions', 'solves')] + [('final_time', D), ('constraint_init_iterations', Z)]


class Rosenbrock(C.Structure):
    _fields_ = [('stages', I), ('new_function', I * 6), ('max_steps', Z)] + [
        (n, D * count) for n, count in [('a', 15), ('c', 15), ('m', 6), ('e', 6), ('alpha', 6), ('gamma', 6)]]
    _fields_ += [(n, D) for n in ('order', 'factor_min', 'factor_max', 'rejection_factor', 'safety', 'h_min', 'h_max', 'h_start')]
    _fields_ += [('controller', P), ('controller_context', P), ('constraint_init_max_iterations', Z),
                 ('constraint_init_max_backtracks', Z)]
    _fields_ += [(n, D) for n in ('constraint_init_tolerance', 'constraint_init_backtrack_factor', 'constraint_init_sufficient_decrease')]


class BackwardEuler(C.Structure):
    _fields_ = [('small', D), ('h_start', D), ('max_iterations', Z), ('reduction_count', Z), ('reductions', D * 5)]


class Linear(C.Structure):
    _fields_ = [('factor', P), ('solve', P)]


lib = C.CDLL(str(Path(__file__).resolve().parent.parent / 'build/libchem.so'))

def bind(name, result, *args):
    f = getattr(lib, name)
    f.restype, f.argtypes = result, args
    return f


bind('chem_rate_default', Rate, I)
bind('chem_rate_evaluate', D, C.POINTER(Rate), D, D, D, DP)
bind('chem_constraint_default', Constraint, I)
bind('chem_mechanism_create', I, C.POINTER(Data), C.POINTER(P))
bind('chem_mechanism_free', None, P)
bind('chem_rate_by_name', C.POINTER(Rate), P, C.c_char_p)
bind('chem_mechanism_data', C.POINTER(Data), P)
bind('chem_state_create', I, P, Z, C.POINTER(P))
bind('chem_state_free', None, P)
bind('chem_is_algebraic', I, P, Z)
bind('chem_constraint_nonzero', I, P, Z, Z)
bind('chem_jacobian_size', Z, P)
bind('chem_jacobian_index', Z, P, Z, Z)
for name in ('concentrations', 'parameters', 'temperatures', 'pressures', 'air_densities', 'absolute_tolerances'):
    bind('chem_' + name, DP, P)
bind('chem_rate_constants', DP, P)
bind('chem_relative_tolerance', None, P, D)
bind('chem_update_rates', I, P)
bind('chem_calculate_air_density', I, P)
bind('chem_status_string', C.c_char_p, I)
bind('chem_alpha_minus_jacobian', None, P, Z, DP, D)
for name in ('forcing', 'negative_jacobian', 'constraint_forcing', 'constraint_negative_jacobian'):
    bind('chem_' + name, None, P, DP, DP)
bind('chem_rosenbrock_parameters', Rosenbrock, I)
bind('chem_solve', I, P, C.POINTER(Rosenbrock), D, C.POINTER(Stats))
bind('chem_initialize_constraints', I, P, C.POINTER(Rosenbrock), C.POINTER(Stats))
bind('chem_backward_euler_parameters', BackwardEuler)
bind('chem_solve_backward_euler', I, P, C.POINTER(BackwardEuler), D, C.POINTER(Stats))
bind('chem_backward_euler_converged', I, P, C.POINTER(BackwardEuler), DP, DP)
bind('chem_normalized_error', D, P, DP, DP, DP)
for name in ('mozart', 'doolittle'):
    bind('chem_linear_' + name, Linear)
bind('chem_set_linear_solver', I, P, Linear)


class Mechanism:
    """Own all borrowed descriptor arrays for the lifetime of the C mechanism."""
    def __init__(self, names, reactions=(), constraints=(), parameters=(), cells=1, expected_status=0, atols=None):
        self.keep = []
        self.n, self.nc = len(names), cells
        def array(typ, values):
            result = (typ * len(values))(*values)
            self.keep.append(result)
            return result
        self.array = array
        rx = []
        for reactants, products, rate in reactions:
            rr, pp = array(I, reactants), array(Product, [Product(*v) for v in products])
            rx.append(Reaction(len(rr), len(pp), rr, pp, rate))
        cs = []
        for constraint, reactants, products in constraints:
            rr, pp = array(Product, [Product(*v) for v in reactants]), array(Product, [Product(*v) for v in products])
            constraint.reactant_count, constraint.product_count = len(rr), len(pp)
            constraint.reactants, constraint.products = rr, pp
            cs.append(constraint)
        self.data = Data(species_count=len(names), reaction_count=len(rx), parameter_count=len(parameters),
            species_names=array(C.c_char_p, [n.encode() for n in names]),
            parameter_names=array(C.c_char_p, [n.encode() for n in parameters]),
            temperature=298.15, pressure=101325, reactions=array(Reaction, rx),
            absolute_tolerances=array(D,atols) if atols is not None else None,
            constraint_count=len(cs), constraints=array(Constraint, cs))
        self.handle, self.state = P(), P()
        status = lib.chem_mechanism_create(C.byref(self.data), C.byref(self.handle))
        assert status == expected_status, (status, expected_status)
        if status:
            return
        assert lib.chem_state_create(self.handle, cells, C.byref(self.state)) == 0
        self.y = lib.chem_concentrations(self.state)
        self.custom = lib.chem_parameters(self.state)

    def close(self):
        lib.chem_state_free(self.state)
        lib.chem_mechanism_free(self.handle)
        self.state, self.handle = P(), P()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def set_y(self, rows):
        assert len(rows) == self.nc
        for c, row in enumerate(rows):
            assert len(row) == self.n
            for i, v in enumerate(row):
                self.y[i * self.nc + c] = v

    def rows(self):
        # Solve can swap the concentration and candidate allocations.
        self.y = lib.chem_concentrations(self.state)
        return [[self.y[i * self.nc + c] for i in range(self.n)] for c in range(self.nc)]

    def force(self, constraints=False):
        out = (D * (self.n * self.nc))()
        f = lib.chem_constraint_forcing if constraints else lib.chem_forcing
        f(self.state, self.y, out)
        return [[out[i * self.nc + c] for i in range(self.n)] for c in range(self.nc)]

    def jacobian(self, constraints=False):
        count = lib.chem_jacobian_size(self.handle)
        out = (D * (count * self.nc))()
        f = lib.chem_constraint_negative_jacobian if constraints else lib.chem_negative_jacobian
        f(self.state, self.y, out)
        result = []
        missing = Z(-1).value
        for c in range(self.nc):
            matrix = []
            for row in range(self.n):
                values = []
                for col in range(self.n):
                    index = lib.chem_jacobian_index(self.handle, row, col)
                    values.append(0 if index == missing else out[index * self.nc + c])
                matrix.append(values)
            result.append(matrix)
        return result

    def solve(self, dt, method=3, parameters=None):
        p = parameters or lib.chem_rosenbrock_parameters(method)
        stats = Stats()
        status = lib.chem_solve(self.state, C.byref(p), dt, C.byref(stats))
        self.rows()
        return status, stats


def rate(kind, values=None):
    r = lib.chem_rate_default(kind)
    if values is not None:
        for i, v in enumerate(values):
            r.p[i] = v
    return r


def equilibrium(row, reactants, products, k=.033, enthalpy=0):
    c = lib.chem_constraint_default(1)
    c.algebraic_species, c.equilibrium_constant, c.enthalpy = row, k, enthalpy
    return c, reactants, products


def linear(row, terms, constant):
    c = lib.chem_constraint_default(0)
    c.algebraic_species, c.constant = row, constant
    return c, terms, []
